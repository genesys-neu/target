#!/usr/bin/env python

#*********************************************************************************
# DISTRIBUTION STATEMENT A. Approved for public release. Distribution is unlimited.
#
# This material is based upon work supported under Air Force Contract No. FA8702-15-D-0001. Any opinions, findings, conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the U.S. Air Force.
#
# (c) 2023 Massachusetts Institute of Technology.
#
# Subject to FAR52.227-11 Patent Rights - Ownership by the contractor (May 2014)
#
# The software/firmware is provided to you on an As-Is basis
#
# Delivered to the U.S. Government with Unlimited Rights, as defined in DFARS Part 252.227-7013 or 7014 (Feb 2014). Notwithstanding any copyright notice, U.S. Government rights in this work are defined by DFARS 252.227-7013 or DFARS 252.227-7014 as detailed above. Use of this work other than as specifically authorized by the U.S. Government may violate any copyrights that exist in this work.
#********************************************************************************/

# This is the main program for demonstration of distributed retrodirective
# beamforming. It must be run with a node ID argument (see main function
# below). In order for data to be logged to file, run this via the shell
# wrapper script (test_node.sh).

import sys
import signal
import threading
import time
import os
import datetime

import numpy as np
import pyfftw
from scipy.io import savemat, loadmat
from scipy import signal as sp
from gnuradio import gr, blocks, uhd
from gnuradio.blocks import parse_file_metadata
from gnuradio.filter import firdes
import pmt

def hdr_to_mat(file_in):

    # Convert .hdr to .mat file for post-processing

    # Initialize the output dict
    data_out = {'rx_freq': [], 'rx_rate': [], 'rx_time': [], 'n_items': [], 'n_bytes': []}

    fid = open(file_in, "rb")
    n_headers = 0
    n_read = 0
    while(True):

        # Read in next header
        hdr_start = fid.tell()
        header_str = fid.read(parse_file_metadata.HEADER_LENGTH)
        if len(header_str) == 0:
            break

        # Convert from string to PMT
        try:
            header = pmt.deserialize_str(header_str)
        except RuntimeError:
            sys.stderr.write('Could not deserialize header: invalid or corrupt data file.\n')
            sys.exit(1)

        # Parse the header
        info = parse_file_metadata.parse_header(header, False)

        # Parse any extra header
        if info['extra_len'] > 0:

            # Read in
            extra_str = fid.read(info['extra_len'])
            if len(extra_str) == 0:
                break

            # Convert from string to PMT
            try:
                extra = pmt.deserialize_str(extra_str)
            except RuntimeError:
                sys.stderr.write('Could not deserialize extra header: invalid or corrupt data file.\n')
                sys.exit(1)

            # Parse the extra header and save back into info
            info = parse_file_metadata.parse_extra_dict(extra, info, False)

        # NOTE: may want to change this to save any extra tags below. Right now
        # we only save the tags from the USRP.

        # Add to the output dict
        if 'rx_freq' in info:
            data_out['rx_freq'].append(pmt.to_double(info['rx_freq']))
        if 'rx_rate' in info:
            data_out['rx_rate'].append(info['rx_rate'])
        if 'rx_time' in info:
            data_out['rx_time'].append(info['rx_time'])
        if 'nitems' in info:
            data_out['n_items'].append(info['nitems'])
        if 'nbytes' in info:
            data_out['n_bytes'].append(info['nbytes'])

        # Move on to the next header
        n_headers += 1
        n_read += parse_file_metadata.HEADER_LENGTH + info["extra_len"]
        fid.seek(n_read, 0)

    # Save output file
    file_out = file_in + '.mat'
    savemat(file_out, data_out)

def make_const_signal(n_samps, n_pad):

    # Makes a constant signal and pads with zeros
    # Output power is 90% of maximum
    x = [1+0j,]*n_samps
    x = 0.9*np.concatenate([np.array(x), np.zeros(n_pad)])
    return x

def make_tone_signal(n_samps, n_pad):

    # Makes a tone signal and pads with zeros
    # Frequency is 1% of the sample rate
    # Output power is 90% of maximum
    f = 0.01
    n = np.arange(n_samps)
    x = np.exp(1j*2*np.pi*f*n)
    x = 0.9*np.concatenate([x, np.zeros(n_pad)])
    return x

def get_gc(node_id, n_bits):

    # Get Gold code from externally generated .mat file
    data = loadmat('../../matlab/gc.mat')
    b_mat = data['b_mat']

    # Basic error checking
    if np.shape(b_mat)[1] < n_bits:
        raise RuntimeError('Not enough bits in reference matrix')
    elif np.shape(b_mat)[0] < (node_id+1):
        raise RuntimeError('Not enough codes in reference matrix')

    # Return bits as bools
    bits = b_mat[node_id, :n_bits]
    return bits.astype(bool)

def make_bpsk_signal(node_id, n_samps, n_pad, payload_bits=None):

    # Makes a BPSK signal with bits based on the node's ID
    # Output power is 90% of maximum
    # NOTE: it could be useful to make this a class, so we can use the modulator
    # without having to create the preamble, filter, etc. every time.

    # Samples per symbol
    samps_per_sym = 4

    # Number of symbols (same as number of bits for BPSK)
    n_syms = int(n_samps/samps_per_sym)

    # Generate bits from Gold code or at random
    gold = 1
    if gold:
        bits = get_gc(node_id, n_syms)
    else:
        np.random.seed(node_id)
        bits = np.random.rand(n_syms) > 0.5

    # Add payload bits if we have any
    if payload_bits != None:
        bits = np.concatenate([bits, payload_bits])
        n_samps += (samps_per_sym*len(payload_bits))

    # Map bits to symbols
    syms = 1 - 2*bits

    # Upsample and zero pad
    syms_up = np.zeros(n_samps)
    syms_up[::samps_per_sym] = syms
    syms_up = np.concatenate([syms_up, np.zeros(n_pad)])

    # RRC filter
    alpha = 0.35
    n_taps = 11*samps_per_sym
    rrc_filt = firdes.root_raised_cosine(1, samps_per_sym, 1, alpha, n_taps)
    x = np.convolve(syms_up, rrc_filt)
    x = x[:len(syms_up)]

    # Scale and return
    x = 0.9*(x/np.max(np.abs(x)))
    return x

class usrp_ctrl(gr.basic_block):

    """
    Block that controls a node's transmit/receive operations using a USRP.
    """

    def __init__(self, node_id, rx_usrp):
        gr.basic_block.__init__(
            self,
            name='usrp_ctrl',
            in_sig=None,
            out_sig=None
        )

        # Set node ID
        self.node_id = node_id

        # The USRP receiver. We need a copy of this here since we are stopping
        # the receiver every time we need to transmit.
        self.rx_usrp = rx_usrp

        # Input and output ports
        self.message_port_register_in(pmt.intern('msg_in'))
        self.message_port_register_out(pmt.intern('msg_out'))

        # Message handler
        self.set_msg_handler(pmt.intern('msg_in'), self.process_msg)

        # Main thread
        self.thread = threading.Thread(target=self.ctrl_loop)
        self.thread.daemon = True

        # Sounding signal (1 ms of data)

        # NOTE: padding with 100 usec of zeros here to flush the USRP's buffers
        # during transmit. May want to find a different way of doing that, or
        # put the padding elsewhere for efficiency (avoid processing zeros).
        n_samps = 1000
        n_pad = 100

        # Either a constant, a tone, or BPSK. If we're the leader, send 16
        # payload bits after the preamble. These are initially zero, but are
        # XOR'd with a randomizer so we don't send all zeros:
        #   0x1D53 = 0001 1101 0101 0011

        # NOTE: when using the constant, the RX DC offset correction should be
        # disabled, otherwise it will be filtered out.
        #tx_snd_sig = make_const_signal(n_samps, n_pad)
        #tx_snd_sig = make_tone_signal(n_samps, n_pad)

        if self.node_id == 1:
            phase_est = 0
            phase_int_rand = round(phase_est*10)^0x1D53
            payload_bits = np.binary_repr(phase_int_rand, 16)
            payload_bits = [bool(int(i)) for i in payload_bits]
            tx_snd_sig = make_bpsk_signal(node_id, n_samps, n_pad, payload_bits)
            self.tx_snd_sig_ready = True
        else:
            tx_snd_sig = make_bpsk_signal(node_id, n_samps, n_pad)
        self.tx_snd_sig = pmt.init_c32vector(len(tx_snd_sig), tx_snd_sig)

        # Beamforming signal (1 ms of data)
        # See comments above about padding. All nodes use a random seed of 9
        # to generate this signal, so they transmit the same data.
        # NOTE: this means the maximum number of nodes is 10 (target, leader,
        # and eight followers)
        tx_bf_sig = make_bpsk_signal(9, n_samps, n_pad)
        self.tx_bf_sig = pmt.init_c32vector(len(tx_bf_sig), tx_bf_sig)

        # Initialize the TX metadata
        self.tx_freq_shift = 0
        self.tx_phase_shift = 0

        # Follower needs three phase terms for beamforming. The order is:
        #   [leader phase, self phase, target phase]
        if self.node_id >= 2:
            self.phases = np.array([None,]*3)

        # Followers must time sync to the leader
        # (Can be disabled for debugging, make sure all have the same PPS)
        if self.node_id >= 2:
            self.have_epoch_time = False

        # Epoch time and period (sec)
        self.epoch_time = 0
        self.epoch_period = 0.1

    def start(self):

        # Start the main thread
        self.thread.start()
        return gr.basic_block.start(self)

    def process_msg(self, msg):

        # If we're the leader
        if self.node_id == 1:

            # Get node detected
            node_det = pmt.dict_ref(msg, pmt.intern('node_det'), pmt.PMT_NIL)
            if node_det == pmt.PMT_NIL:
                print('WARNING: missing node_det key')
                node_det = 'none'
            else:
                node_det = pmt.symbol_to_string(node_det)

            # Get phase estimate
            phase_est = pmt.dict_ref(msg, pmt.intern('phase_est'), pmt.PMT_NIL)
            if phase_est == pmt.PMT_NIL:
                print('WARNING: missing phase_est key')
                phase_est = 0
            else:
                phase_est = pmt.to_double(phase_est)

            # Set TX phase to beamform to the target
            if node_det == 'target':
                self.tx_phase_shift = -phase_est
                #print('phase update at %s' % self.rx_usrp.get_time_now().get_real_secs())

            # Update sounding waveform with estimate of follower's phase
            # NOTE: this needs updating if we ever have more than one follower
            if node_det == 'follower':
                # Phase is encoded into 16 bits, integer part with one decimal.
                # See above for explanation about the randomizer.
                phase_int_rand = round(phase_est*10)^0x1D53
                payload_bits = np.binary_repr(phase_int_rand, 16)
                payload_bits = [bool(int(i)) for i in payload_bits]
                n_samps = 1000
                n_pad = 100
                tx_snd_sig = make_bpsk_signal(self.node_id, n_samps, n_pad, payload_bits)
                self.tx_snd_sig = pmt.init_c32vector(len(tx_snd_sig), tx_snd_sig)
                self.tx_snd_sig_ready = True

        # Otherwise if we're a follower
        elif self.node_id >= 2:

            # Get node detected
            node_det = pmt.dict_ref(msg, pmt.intern('node_det'), pmt.PMT_NIL)
            if node_det == pmt.PMT_NIL:
                print('WARNING: missing node_det key')
                node_det = 'none'
            else:
                node_det = pmt.symbol_to_string(node_det)

            # Get phase estimate
            phase_est = pmt.dict_ref(msg, pmt.intern('phase_est'), pmt.PMT_NIL)
            if phase_est == pmt.PMT_NIL:
                print('WARNING: missing phase_est key')
                phase_est = 0
            else:
                phase_est = pmt.to_double(phase_est)

            # Get the packet detection time. The first one is based on the RX
            # time tag from the USRP and the first sample that crosses the
            # detection threshold. The second one is relative to the start of
            # the reference waveform. So the full time is the sum of these.
            pkt_time = pmt.dict_ref(msg, pmt.intern('pkt_time'), pmt.PMT_NIL)
            if pkt_time == pmt.PMT_NIL:
                print('WARNING: missing pkt_time key')
                pkt_time = 0
            else:
                pkt_time = pmt.to_double(pkt_time)
            time_est = pmt.dict_ref(msg, pmt.intern('time_est'), pmt.PMT_NIL)
            if time_est == pmt.PMT_NIL:
                print('WARNING: missing time_est key')
                time_est = 0
            else:
                time_est = pmt.to_double(time_est)
            det_time = pkt_time + time_est

            # If the message was from the leader
            if node_det == 'leader':

                # Save the phase estimate
                self.phases[0] = phase_est

                # Set TX frequency to match the leader
                freq_est = pmt.dict_ref(msg, pmt.intern('freq_est'), pmt.PMT_NIL)
                if freq_est == pmt.PMT_NIL:
                    print('WARNING: missing freq_est key')
                    freq_est = 0
                else:
                    freq_est = pmt.to_double(freq_est)
                self.tx_freq_shift = freq_est

                # NOTE: ideally we would update the epoch time every packet to
                # handle any distance changes between the nodes. But right now
                # we only set the epoch time based on the first detection and
                # let the main loop handle updating the rest of the time.

                # NOTE: the extra delay of 47.1 usec is due to the USRP's DSP.
                # This was measured with two devices on the same clock and a
                # very short cable between the two TX/RX.

                # Use the time estimate to adjust our clock. We set the time to
                # enter into the network at the start of the next epoch, which
                # is in 75 ms. If that's too soon it will be caught by the main
                # loop and adjusted until we have enough time to enter.
                if not self.have_epoch_time:
                    self.epoch_time = det_time + 0.075 - 47.1e-6

                # If this is the first we've obtained time, adjust the RX to
                # deliver 1.25 ms of data from now on
                if not self.have_epoch_time:
                    self.rx_stream_cmd.num_samps = 0.00125*1e6
                    self.have_epoch_time = True
                    print('Follower obtained epoch time (%s)' % self.epoch_time)

                # Get the phase feedback
                phase_fb = pmt.dict_ref(msg, pmt.intern('phase_fb'), pmt.PMT_NIL)
                if phase_fb == pmt.PMT_NIL:
                    print('WARNING: missing phase_fb key')
                    phase_fb = 0
                else:
                    phase_fb = pmt.to_double(phase_fb)
                if self.phases[1] == phase_fb:
                    print('WARNING: received the same phase feedback as before')
                self.phases[1] = phase_fb

            # Otherwise, if the message was from the target
            elif node_det == 'target':

                # Save the phase
                self.phases[2] = phase_est

                # Compute the delay between the expected and actual detection
                # times for the target. This is needed to adjust the phase based
                # on our frequency offset. We also need to compensate for the
                # USRP's DSP delay when doing this.
                target_delay = (self.epoch_time + 0.050) - det_time + 47.1e-6

                # If we have all three phases, set the beamforming phase. Note
                # that we must adjust this by any frequency offset relative to
                # the leader, to account for the leader's phase at various times
                # during the protocol. The phases are:
                #
                #   p0 / phases[0] : leader-to-follower @ 25 ms. Received by the
                #   followers without downconversion at the right phase, so
                #   adjust by -2*pi*f*0.025.
                #
                #   p1 / phases[1] : follower-to-leader @ 0 ms. Sent by the
                #   followers without upconversion at the right phase, so adjust
                #   by 2*pi*f*0.000 (at time 0, so the phase is 0).
                #
                #   p2 / phases[2] : target-to-follower @ 50 ms. Received by the
                #   followers without downconversion at the right phase, and
                #   potentially with an extra delay from the target, so adjust
                #   by -2*pi*f*(0.050 + d) where d is the extra delay.
                #
                # Also, for the beamforming transmission @ 75 ms, the followers
                # must upconvert with the right phase: 2*pi*f*0.075. So putting
                # this together, along with the signs in the below, the phase
                # factor is:
                #   (p0 - 2*pi*f*0.025) - (p1 + 2*pi*f*0.00) - (p2 - 2*pi*f*(0.050 + d)) + 2*pi*f*0.075
                #   = (p0 - p1 - p2) + 2*pi*f*(-0.025 - 0.000 + (0.050 + d) + 0.075)
                #   = (p0 - p1 - p2) + 2*pi*f*(0.100 + d)
                if not np.any(self.phases == None):
                    self.tx_phase_shift = self.phases[0] - self.phases[1] - self.phases[2]
                    self.tx_phase_shift += np.degrees(2*np.pi*self.tx_freq_shift*(0.1 + target_delay))
                    #print('phase update at %s' % self.rx_usrp.get_time_now().get_real_secs())

    def ctrl_loop(self):

        # Wait for a PPS event
        t = self.rx_usrp.get_time_last_pps()
        while(t == self.rx_usrp.get_time_last_pps()):
            time.sleep(0.1)
        time.sleep(0.1)

        # Get the start of the next second (epoch time)
        self.epoch_time = self.rx_usrp.get_time_last_pps().get_full_secs() + 1.0

        # Initialize the RX stream command
        self.rx_stream_cmd = uhd.stream_cmd(uhd.stream_cmd_t.STREAM_MODE_NUM_SAMPS_AND_DONE)
        self.rx_stream_cmd.stream_now = False

        # NOTE: do an initial streaming operation to get things started (1.25 ms
        # of data). Without this the first two packets are returned at once and
        # only with one tag. I'm not sure why, but this prevents that, so we'll
        # have a time tag on every packet.
        self.rx_stream_cmd.num_samps = 0.00125*1e6
        rx_time = self.rx_usrp.get_time_now().get_real_secs() + 0.01
        self.rx_stream_cmd.time_spec = uhd.time_spec_t(rx_time)
        self.rx_usrp.issue_stream_cmd(self.rx_stream_cmd)

        # If we're a follower, adjust the stream command so we receive an entire
        # epoch's worth of samples. This is needed so we can detect the leader
        # and schedule our TX/RX operations. Note that the second part of this
        # if statement should only be true if we've been told not to do over the
        # air timing synchronization during initialization (so for debugging).
        if (self.node_id >= 2) and (not self.have_epoch_time):
            # NOTE: this must still be a multiple of 1.25 ms packets, which
            # we've set the USRP up to deliver
            self.rx_stream_cmd.num_samps = self.epoch_period*1e6

        # Run a loop forever, depending on what type of node we are

        # Target
        if self.node_id == 0:
            print('Node type: target')
            while True:

                # Receive sounding from followers @ 0 ms
                rx_time = self.epoch_time + 0.000
                self.wait_until_before(rx_time)
                self.do_rx(rx_time)

                # Receive sounding from leader @ 25  ms
                rx_time = self.epoch_time + 0.025
                self.wait_until_before(rx_time)
                self.do_rx(rx_time)

                # Transmit sounding to leader and followers @ 50 ms
                tx_time = self.epoch_time + 0.050
                self.wait_until_before(tx_time)
                self.do_tx(tx_time, self.tx_snd_sig)

                # Receive beamforming payload @ 75 ms
                rx_time = self.epoch_time + 0.075
                self.wait_until_before(rx_time)
                self.do_rx(rx_time)

                # Update for the next epoch
                self.epoch_time += self.epoch_period

        # Leader
        elif self.node_id == 1:
            print('Node type: leader')
            while True:

                # Receive sounding from followers @ 0 ms
                rx_time = self.epoch_time + 0.000
                self.wait_until_before(rx_time)
                self.do_rx(rx_time)

                # Transmit sounding to followers @ 25 ms
                tx_time = self.epoch_time + 0.025
                self.wait_until_before(tx_time)
                if not self.tx_snd_sig_ready:
                    print('WARNING: sending the same phase feedback as before')
                self.do_tx(tx_time, self.tx_snd_sig)
                self.tx_snd_sig_ready = False

                # Receive from the target @ 50 ms
                rx_time = self.epoch_time + 0.050
                self.wait_until_before(rx_time)
                self.do_rx(rx_time)

                # Transmit beamforming payload to the target @ 75 ms
                # NOTE: leader frequency shift is always 0
                # NOTE: we wait here to make sure we've got the phase update
                tx_time = self.epoch_time + 0.075
                self.wait_until_before(tx_time)
                self.do_tx(tx_time, self.tx_bf_sig, self.tx_freq_shift, self.tx_phase_shift)

                # Update for the next epoch
                self.epoch_time += self.epoch_period

        # Follower
        else:
            print('Node type: follower')
            while True:

                # If we don't have timing from the leader yet
                if not self.have_epoch_time:

                    # Receive an entire epoch to try and detect the leader
                    # This is scheduled two epochs from now so we have a duty
                    # cycle of 50% (otherwise, get error code / late command)
                    rx_time = self.rx_usrp.get_time_now().get_real_secs() + 2*self.epoch_period
                    self.wait_until_before(rx_time)
                    self.do_rx(rx_time)

                # Otherwise, if we have timing from the leader, operate normally
                else:

                    # Since we get epoch time from the other thread, we may be
                    # blocked in the above RX call, which takes place in two
                    # epochs. By the time we make it here we may already be
                    # behind, so jump ahead if necessary.
                    time_margin = self.epoch_time - self.rx_usrp.get_time_now().get_real_secs()
                    while time_margin < 0.010:
                        print('WARNING: next epoch in %.3f ms, jumping ahead by an epoch' % (time_margin*1e3))
                        self.epoch_time += self.epoch_period
                        time_margin = self.epoch_time - self.rx_usrp.get_time_now().get_real_secs()

                    # Transmit sounding to the leader @ 0 ms
                    # NOTE: wait here to make sure we've got the freq update
                    tx_time = self.epoch_time + 0.000
                    self.wait_until_before(tx_time)
                    self.do_tx(tx_time, self.tx_snd_sig, self.tx_freq_shift)

                    # Receive sounding from the leader @ 25  ms
                    rx_time = self.epoch_time + 0.025
                    self.wait_until_before(rx_time)
                    self.do_rx(rx_time)

                    # Receive from the target @ 50 ms
                    rx_time = self.epoch_time + 0.050
                    self.wait_until_before(rx_time)
                    self.do_rx(rx_time)

                    # Transmit beamforming payload to the target @ 75 ms
                    # NOTE: wait here to make sure we've got the freq/phase updates
                    tx_time = self.epoch_time + 0.075
                    self.wait_until_before(tx_time)
                    self.do_tx(tx_time, self.tx_bf_sig, self.tx_freq_shift, self.tx_phase_shift)

                    # NOTE: updating the epoch time here instead of in the other
                    # thread based on leader detection times.
                    self.epoch_time += self.epoch_period

    def do_tx(self, tx_time, tx_sig, tx_freq_shift=0, tx_phase_shift=0):

        # Creates time tag in the format required by the USRP
        int_time = int(tx_time)
        frac_time = tx_time - int_time
        time_tuple = pmt.make_tuple(pmt.from_uint64(int_time), pmt.from_double(frac_time))
        tx_meta = pmt.dict_add(pmt.make_dict(), pmt.intern('tx_time'), time_tuple)

        # Add our frequency and phase shift tags (Hz and degrees)
        tx_meta = pmt.dict_add(tx_meta, pmt.intern('freq_shift'), pmt.from_double(tx_freq_shift))
        tx_meta = pmt.dict_add(tx_meta, pmt.intern('phase_shift'), pmt.from_double(tx_phase_shift))

        # Form the PDU and send it to the USRP
        tx_pdu = pmt.cons(tx_meta, tx_sig)
        #print('sending for %s at %s' % (tx_time, self.rx_usrp.get_time_now().get_real_secs()))
        self.message_port_pub(pmt.intern('msg_out'), tx_pdu)

    def do_rx(self, rx_time):

        # Update the stream command and send it to the USRP
        self.rx_stream_cmd.time_spec = uhd.time_spec_t(rx_time)
        self.rx_usrp.issue_stream_cmd(self.rx_stream_cmd)

    def wait_until_before(self, t, t_margin=0.0075):

        # Waits until t_margin seconds before t, based on the USRP's time
        t_wait = t - self.rx_usrp.get_time_now().get_real_secs() - t_margin
        if t_wait <= 0:
            # NOTE: in this case we may want to skip ahead to the next epoch
            print('WARNING: not enough time margin')
        else:
            time.sleep(t_wait)

        # All done
        return

class freq_phase_shift(gr.basic_block):

    """
    Block that performs a frequency and phase shift on the input signal.
    """

    def __init__(self, samp_rate, node_id):
        gr.basic_block.__init__(
            self,
            name='freq_phase_shift',
            in_sig=None,
            out_sig=None
        )

        # Save sample rate and node ID
        self.samp_rate = samp_rate
        self.node_id = node_id

        # Input and output ports
        self.message_port_register_in(pmt.intern('msg_in'))
        self.message_port_register_out(pmt.intern('msg_out'))

        # Message handler (main processing function)
        self.set_msg_handler(pmt.intern('msg_in'), self.process_msg)

        # Time vector for frequency shift (maximum packet size is 1.25 ms),
        # only needed if we're a follower
        if self.node_id >= 2:
            self.n = np.arange(1250)

    def process_msg(self, msg):

        # Get metadata
        # NOTE: only followers do a frequency shift
        meta = pmt.car(msg)
        if self.node_id >= 2:
            freq_shift = pmt.dict_ref(meta, pmt.intern('freq_shift'), pmt.PMT_NIL)
            if freq_shift == pmt.PMT_NIL:
                print('WARNING: missing freq_shift key')
            else:
                freq_shift = pmt.to_double(freq_shift)
        phase_shift = pmt.dict_ref(meta, pmt.intern('phase_shift'), pmt.PMT_NIL)
        if phase_shift == pmt.PMT_NIL:
            print('WARNING: missing phase_shift key')
        else:
            phase_shift = pmt.to_double(phase_shift)

        # Get data
        data = pmt.cdr(msg)
        data = pmt.c32vector_elements(data)

        # Followers do frequency and phase shift, leader only does phase shift
        if self.node_id >= 2:
            shift_sig = np.exp(1j*(2*np.pi*(freq_shift/self.samp_rate)*self.n[:len(data)] + np.radians(phase_shift)))
        elif self.node_id == 1:
            shift_sig = np.exp(1j*np.radians(phase_shift))

        # Apply the shift
        data_shift = shift_sig*np.array(data)

        # Output data and message
        data_shift = pmt.init_c32vector(len(data_shift), data_shift)
        msg_out = pmt.cons(meta, data_shift)

        # Publish message
        self.message_port_pub(pmt.intern('msg_out'), msg_out)

class pkt_detector(gr.basic_block):

    """
    Block that performs power detection on a stream of input samples and returns
    packets of data for further processing.
    """

    def __init__(self, samp_rate, pow_thresh, t_min1, t_min2, t_max):
        gr.sync_block.__init__(
            self,
            name='pkt_detector',
            in_sig=[np.complex64],
            out_sig=None
        )

        # Save inputs
        self.samp_period = 1/samp_rate
        self.pow_thresh = pow_thresh
        self.n_min1 = int(t_min1*samp_rate)
        self.n_min2 = int(t_min2*samp_rate)
        self.n_max = int(t_max*samp_rate)

        # Initialize
        # NOTE: packet buffer has to be complex128 for conversion to c32vector
        self.rx_time = 0
        self.meta = None
        self.buff = np.empty(self.n_max, dtype=np.complex128)
        self.n_pkt = 0
        self.n_samps = 0
        self.n_high = 0
        self.n_low = 0

        # Output port
        self.message_port_register_out(pmt.intern('msg_out'))

    def work(self, input_items, output_items):

        # Look for rx_time tags from the USRP. If one is found, update the time
        # and clear the sample counter
        n_in = len(input_items[0])
        tags = self.get_tags_in_window(0, 0, n_in, pmt.intern('rx_time'))
        n_tags = len(tags)
        if n_tags > 1:
            print('WARNING: found more than one rx_time tag')
        elif n_tags == 1:
            tag = tags[0]
            offset = tag.offset - self.nitems_read(0)
            if offset != 0:
                print('WARNING: rx_time tag not on first sample')
            int_time = pmt.tuple_ref(tag.value, 0)
            frac_time = pmt.tuple_ref(tag.value, 1)
            self.rx_time = pmt.to_uint64(int_time) + pmt.to_double(frac_time)
            self.n_samps = 0

        # Find input samples that exceed the power threshold
        x = input_items[0]
        p = np.abs(x)**2 > self.pow_thresh

        # Process samples one at a time
        # NOTE: this could be vectorized to make it more efficient.
        for i in range(n_in):

            # If power has been detected
            if p[i]:

                # Update the counters
                self.n_high += 1
                self.n_low = 0

                # If this is the first sample of the packet, make our metadata
                if self.n_pkt == 0:
                    pkt_time = self.rx_time + self.n_samps*self.samp_period 
                    self.meta = pmt.dict_add(pmt.make_dict(), pmt.intern('pkt_time'), pmt.from_double(pkt_time))

                # Add to the packet if there is space, otherwise send it
                if self.n_pkt < self.n_max:
                    self.buff[self.n_pkt] = x[i]
                    self.n_pkt += 1
                else:
                    print('WARNING: packet exceeded maximum size (%i), sending now' % self.n_max)
                    data = pmt.init_c32vector(self.n_pkt, self.buff[:self.n_pkt])
                    msg_out = pmt.cons(self.meta, data)
                    self.message_port_pub(pmt.intern('msg_out'), msg_out)
                    #print('Published packet with %i samples' % self.n_pkt)
                    self.n_pkt = 0
                    self.meta = None

            # Otherwise, if no power was detected and there is a packet buffered
            elif self.n_pkt > 0:

                # Update the counters
                self.n_high = 0
                self.n_low += 1

                # If we've gone enough samples without power being detected
                if self.n_low >= self.n_min1:

                    # Send the packet if it's large enough, otherwise drop it
                    if self.n_pkt >= self.n_min2:
                        data = pmt.init_c32vector(self.n_pkt, self.buff[:self.n_pkt])
                        msg_out = pmt.cons(self.meta, data)
                        self.message_port_pub(pmt.intern('msg_out'), msg_out)
                        #print('Published packet with %i samples' % self.n_pkt)
                    #else:
                        #print('WARNING: packet less than minimum size (%i, %i), dropping' % (self.n_pkt, self.n_min2))
                    self.n_pkt = 0
                    self.meta = None

                # Otherwise, continue to add to the packet if there's space
                else:
                    if self.n_pkt < self.n_max:
                        self.buff[self.n_pkt] = x[i]
                        self.n_pkt += 1
                    else:
                        print('WARNING: packet exceeded maximum size (%i), sending now' % self.n_max)
                        data = pmt.init_c32vector(self.n_pkt, self.buff[:self.n_pkt])
                        msg_out = pmt.cons(self.meta, data)
                        self.message_port_pub(pmt.intern('msg_out'), msg_out)
                        #print('Published packet with %i samples' % self.n_pkt)
                        self.n_pkt = 0
                        self.meta = None

            # Otherwise, no power detected and no packet buffered, nothing to do

            # Update sample count
            self.n_samps += 1

        # Consume the input buffer and return (output is PDUs, so no samples
        # into the output buffer)
        self.consume(0, n_in)
        return 0

class pkt_estimator(gr.basic_block):

    """
    Block that performs parameter estimation on the received packet (frequency,
    time, and phase of arrival) for synchronization and beamforming.
    """

    def __init__(self, samp_rate, node_id):
        gr.sync_block.__init__(
            self,
            name='pkt_estimator',
            in_sig=None,
            out_sig=None
        )

        # Save inputs
        self.node_id = node_id
        self.samp_rate = samp_rate

        # Input and output ports
        self.message_port_register_in(pmt.intern('msg_in'))
        self.message_port_register_out(pmt.intern('msg_out'))

        # Message handler (main processing function)
        self.set_msg_handler(pmt.intern('msg_in'), self.process_msg)

        # Create reference signals
        n_samps = 1000
        n_pad = 100

        # If we're the leader, look for the followers and target
        if node_id == 1:
            # NOTE: if more UAVs are added, this will need to be
            # modified to search for multiple followers.
            self.x_follower = make_bpsk_signal(2, n_samps, n_pad)
            self.x_follower /= np.linalg.norm(self.x_follower)
            self.x_target = make_bpsk_signal(0, n_samps, n_pad)
            self.x_target /= np.linalg.norm(self.x_target)

        # If we're a follower, look for the leader and target
        elif node_id >= 2:
            self.x_leader = make_bpsk_signal(1, n_samps, n_pad)
            self.x_leader /= np.linalg.norm(self.x_leader)
            self.x_target = make_bpsk_signal(0, n_samps, n_pad)
            self.x_target /= np.linalg.norm(self.x_target)

        # If we're the target, look for the leader, follower, and beamforming
        else:
            self.x_leader = make_bpsk_signal(1, n_samps, n_pad)
            self.x_leader /= np.linalg.norm(self.x_leader)
            self.x_follower = make_bpsk_signal(2, n_samps, n_pad)
            self.x_follower /= np.linalg.norm(self.x_follower)
            self.x_bf = make_bpsk_signal(9, n_samps, n_pad)
            self.x_bf /= np.linalg.norm(self.x_bf)

        # Frequency axis for FFT-based frequency offset estimator
        # This runs at a sample rate 100x lower (10 KSPS), 0.15 Hz resolution
        self.n_dec = 100
        self.n_fft = 2**16
        self.f = ((samp_rate/self.n_dec)/2)*np.arange(-1, 1, 2/self.n_fft)
        # NOTE: do the FFT shift here so we don't need to do it on the data
        # every time. We use the inverse shift to align with the FFT's output.
        self.f = np.fft.ifftshift(self.f)
        self.fft_in = pyfftw.empty_aligned(self.n_fft, dtype='complex128')
        self.fft_out = pyfftw.empty_aligned(self.n_fft, dtype='complex128')
        self.fft_obj = pyfftw.FFTW(self.fft_in, self.fft_out)

        # Time vector for when we do frequency correction, based on the maximum
        # possible packet size of 1.25 ms
        self.n = np.arange(1250)

        # We keep a running sum to average over many estimates. The target only
        # estimates this for the beamforming signal, the leader for the target, 
        # and the followers for both the leader and the target.
        if self.node_id == 0:
            self.bf_freq_est_sum = 0
            self.bf_n_freq_ests = 0
        elif self.node_id == 1:
            self.target_freq_est_sum = 0
            self.target_n_freq_ests = 0
        elif self.node_id >= 2:
            self.leader_freq_est_sum = 0
            self.leader_n_freq_ests = 0
            self.target_freq_est_sum = 0
            self.target_n_freq_ests = 0

        # Detection threshold (normalized between 0 and 1)
        self.det_thresh = 0.1

        # If we're the target, keep track of individual and beamformed powers
        if self.node_id == 0:
            self.leader_pow_est = 0
            self.follower_pow_est = 0
            self.bf_pow_est = 0

    def process_msg(self, msg):

        # Get metadata
        meta = pmt.car(msg)
        pkt_time = pmt.dict_ref(meta, pmt.intern('pkt_time'), pmt.PMT_NIL)
        if pkt_time == pmt.PMT_NIL:
            print('WARNING: missing pkt_time key')
        else:
            pkt_time = pmt.to_double(pkt_time)

        # Get data
        data = pmt.cdr(msg)
        data = pmt.c32vector_elements(data)
        y = np.array(data)

        # Remove any residual mean
        y -= np.mean(y)

        # Estimate power
        pow_est = np.var(y)

        # Normalize
        y /= np.linalg.norm(y)

        # Detect which waveform is present
        # NOTE: this may not work if the frequency offset is too large
        det = 'none'

        # If we're the target
        if self.node_id == 0:

            # Possible signals: leader, follower, beamforming
            r1, _ = xcorr_via_fft(y, self.x_leader)
            r2, _ = xcorr_via_fft(y, self.x_follower)
            r3, _ = xcorr_via_fft(y, self.x_bf)
            r_max = [np.max(np.abs(r1)**2), np.max(np.abs(r2)**2), np.max(np.abs(r3)**2)]
            l = ['leader', 'follower', 'beamforming']
            max_idx = np.argmax(r_max)
            if r_max[max_idx] > self.det_thresh:
                det = l[max_idx]

        # If we're the leader
        elif self.node_id == 1:

            # Possible signals: follower, target
            # NOTE: if more UAVs are added, this will need to be
            # modified to search for multiple followers.
            r1, _ = xcorr_via_fft(y, self.x_follower)
            r2, _ = xcorr_via_fft(y, self.x_target)
            r_max = [np.max(np.abs(r1)**2), np.max(np.abs(r2)**2)]
            l = ['follower', 'target']
            max_idx = np.argmax(r_max)
            if r_max[max_idx] > self.det_thresh:
                det = l[max_idx]

        # If we're a follower
        elif self.node_id >= 2:

            # Possible signals: leader, target
            r1, _ = xcorr_via_fft(y, self.x_leader)
            r2, _ = xcorr_via_fft(y, self.x_target)
            r_max = [np.max(np.abs(r1)**2), np.max(np.abs(r2)**2)]
            l = ['leader', 'target']
            max_idx = np.argmax(r_max)
            if r_max[max_idx] > self.det_thresh:
                det = l[max_idx]

        # If no signal was detected there's nothing to do
        if det == 'none':
            print('WARNING: no signal detected')
            return

        # Target only continues processing if this is a beamformed signal
        if (self.node_id == 0) and (det != 'beamforming'):

            # Save the power estimate for the individual signals
            if det == 'leader':
                self.leader_pow_est = pow_est
            elif det == 'follower':
                self.follower_pow_est = pow_est

            # For debugging: time, power level, and node detected
            print('%.9f, %.4f, %s' % (pkt_time, pow_est, det))
            sys.stdout.flush()

            return

        # For debugging: expected and actual beamformed powers at the target
        if (self.node_id == 0) and (det == 'beamforming'):
            exp_bf_pow_est = (np.sqrt(self.leader_pow_est) + np.sqrt(self.follower_pow_est))**2
            percent_max = pow_est/exp_bf_pow_est
            if percent_max < 1:
                graph = '.'*int(percent_max*100) + ' '*int((1 - percent_max)*100) + '|'
            else:
                graph = '.'*100 + '|' + '.'*int((percent_max - 1)*100)
            print('exp %.4f, act %.4f, %.2f\t%s' % (exp_bf_pow_est, pow_est, percent_max, graph))
            sys.stdout.flush()

        # Estimate the frequency offset and keep a running sum to average over
        # many estimates. Keep one sum per signal we're interested in. Note that
        # the leader does not do any frequency correction for the followers,
        # as we assume they are compensating for that when transmitting.
        freq_est = float('NaN')

        # Target estimates offset for beamforming
        if (self.node_id == 0) and (det == 'beamforming'):
            this_freq_est = self.do_freq_est(y)
            self.bf_freq_est_sum += this_freq_est
            self.bf_n_freq_ests += 1
            freq_est = self.bf_freq_est_sum/self.bf_n_freq_ests

        # Leader estimates offset for target
        elif (self.node_id == 1) and (det == 'target'):
            this_freq_est = self.do_freq_est(y)
            self.target_freq_est_sum += this_freq_est
            self.target_n_freq_ests += 1
            freq_est = self.target_freq_est_sum/self.target_n_freq_ests

        # Followers estimate offset for leader and target
        elif self.node_id >= 2:
            this_freq_est = self.do_freq_est(y)
            if det == 'leader':
                self.leader_freq_est_sum += this_freq_est
                self.leader_n_freq_ests += 1
                freq_est = self.leader_freq_est_sum/self.leader_n_freq_ests
            elif det == 'target':
                self.target_freq_est_sum += this_freq_est
                self.target_n_freq_ests += 1
                freq_est = self.target_freq_est_sum/self.target_n_freq_ests

        # If we have a frequency estimate, correct for it
        if not np.isnan(freq_est):
            y = y*np.exp(1j*2*np.pi*(-freq_est/self.samp_rate)*self.n[:len(y)])

        # Estimate time and phase by correlating with the reference signal
        time_est = float('NaN')
        phase_est = float('NaN')

        # Target estimates offsets for beamforming
        if (self.node_id == 0) and (det == 'beamforming'):
            time_est, phase_est = self.do_time_phase_est(y, self.x_bf)

        # Leader estimates offsets for followers and target
        elif self.node_id == 1:
            if det == 'follower':
                # NOTE: if more UAVs are added, the leader will need to be
                # modified to search for multiple followers.
                time_est, phase_est = self.do_time_phase_est(y, self.x_follower)
            elif det == 'target':
                time_est, phase_est = self.do_time_phase_est(y, self.x_target)

        # Followers estimate offsets for leader and target
        elif self.node_id >= 2:
            if det == 'leader':
                time_est, phase_est = self.do_time_phase_est(y, self.x_leader)
            elif det == 'target':
                time_est, phase_est = self.do_time_phase_est(y, self.x_target)

        # NOTE: could align, estimate SNR, and demodulate the signals here.

        # If we're a follower that detected the leader, decode the last 16 bits
        # to get our phase feedback
        # NOTE: this needs updating if more follower UAVs are added.
        phase_fb = float('NaN')
        if (self.node_id >= 2) and (det == 'leader'):

            # NOTE: skip the first 1000 samples, then keep 16 symbols (4 samples
            # per symbol). The 22 is from the filter delay and -1 is since we're
            # zero-based indexing.
            start_idx = int(time_est*self.samp_rate) + 22 + 1000 - 1
            end_idx = start_idx + 4*16
            y_demod = y[start_idx:end_idx]*np.exp(-1j*np.radians(phase_est))
            """
            # For debugging (copy/paste to MATLAB)
            print('y=[', end='')
            for i in range(len(y_demod)):
                print('%.3f+1j*%.3f,' % (np.real(y_demod[i]), np.imag(y_demod[i])), end='')
            print('];')
            """
            payload_bits = np.real(y_demod[1::4]) < 0
            payload_bits = [str(int(i)) for i in payload_bits]
            payload_bits = ''.join(payload_bits)

            # Need to see if the phase was negative, then derandomize and scale
            if len(payload_bits) == 16:
                if payload_bits[0] == '1':
                    phase_int_rand = int(payload_bits, 2) - 2**16
                else:
                    phase_int_rand = int(payload_bits, 2)
                phase_int = phase_int_rand^0x1D53
                phase_fb = phase_int/10.0
            else:
                print('WARNING: error decoding phase (length = %s), setting to 0' % len(payload_bits))
                phase_int_rand = 0
                phase_int = phase_int_rand^0x1D53
                phase_fb = phase_int/10.0

        # Form output metadata
        meta_out = meta
        #meta_out = pmt.dict_add(meta_out, pmt.intern('pow_est'), pmt.from_double(pow_est))
        meta_out = pmt.dict_add(meta_out, pmt.intern('node_det'), pmt.intern(det))
        meta_out = pmt.dict_add(meta_out, pmt.intern('freq_est'), pmt.from_double(freq_est))
        meta_out = pmt.dict_add(meta_out, pmt.intern('time_est'), pmt.from_double(time_est))
        meta_out = pmt.dict_add(meta_out, pmt.intern('phase_est'), pmt.from_double(phase_est))
        meta_out = pmt.dict_add(meta_out, pmt.intern('phase_fb'), pmt.from_double(phase_fb))

        # Publish message
        self.message_port_pub(pmt.intern('msg_out'), meta_out)

        # For debugging: estimates from this packet (only if we're not the target)
        if self.node_id != 0:
            print('%.9f, %.4f, %s, %.3f, %.9f, %.2f, %.2f' % (pkt_time, pow_est, det, freq_est, time_est, phase_est, phase_fb))
            sys.stdout.flush()

    def do_freq_est(self, y):

        # Frequency estimator for BPSK. Works by squaring and looking for a tone
        # at twice the carrier frequency.

        y_sq = y**2
        y_sq = sp.resample_poly(y_sq, 1, self.n_dec)
        self.fft_in[:] = 0
        self.fft_in[:len(y_sq)] = y_sq
        y_fft = self.fft_obj()
        # NOTE: we've already done the inverse FFT shift on the vector of
        # frequencies, so don't need to shift the data here.
        max_idx = np.argmax(np.abs(y_fft))
        freq_est = self.f[max_idx]/2
        return freq_est

    def do_time_phase_est(self, y, x):


        # Correlate with reference, compute power, and get the timing and phase
        # offset estimates from the peak

        # NOTE: we interpolate for higher resolution, for when the peak is
        # between samples. The efficiency of this could be improved.
        n_interp = 10
        y = sp.resample_poly(y, n_interp, 1)
        x = sp.resample_poly(x, n_interp, 1)
        r, l = xcorr_via_fft(y, x)
        max_idx = np.argmax(np.abs(r)**2)
        time_est = l[max_idx]*(1/(n_interp*self.samp_rate))
        phase_est = np.degrees(np.angle(r[max_idx]))
        return time_est, phase_est

def xcorr(a, b):

    # Cross-correlation - works like xcorr() in MATLAB

    # Zero pad one of the inputs so they are the same length
    n_pad = len(a) - len(b)
    if n_pad > 0:
        b = np.concatenate([b, np.zeros(n_pad)])
    elif n_pad < 0:
        a = np.concatenate([a, np.zeros(-n_pad)])

    # Perform the cross correlation
    r = np.correlate(a, b, 'full')

    # Lags
    max_lag = len(a) - 1
    l = np.arange(-max_lag, max_lag+1)

    # Done
    return r, l

def xcorr_via_fft(a, b):

    # Cross-correlation - works like xcorr() in MATLAB
    # Uses FFT for convolution

    # Zero pad one of the inputs so they are the same length
    n_pad = len(a) - len(b)
    if n_pad > 0:
        b = np.concatenate([b, np.zeros(n_pad)])
    elif n_pad < 0:
        a = np.concatenate([a, np.zeros(-n_pad)])

    # Perform forward FFT, multiply, and inverse FFT
    # Note that we need to zero pad since multiplication of the FFTs is circular
    # convolution but we want linear convolution
    n_pad = len(a) - 1
    a_f = np.fft.fft(np.concatenate([a, np.zeros(n_pad)]))
    b_f = np.fft.fft(np.concatenate([np.conj(b[::-1]), np.zeros(n_pad)]))
    r_f = a_f*b_f
    r = np.fft.ifft(r_f)

    # Lags
    max_lag = len(a) - 1
    l = np.arange(-max_lag, max_lag+1)

    # Done
    return r, l

def sync_usrp_to_external(usrp):

    # Synchronizes the USRP to external time/frequency references

    # Set clock and time sources to the external reference
    usrp.set_clock_source('external', 0)
    usrp.set_time_source('external', 0)

    # Make sure we have 10 MHz lock
    max_n_tries = 3
    for i in range(max_n_tries):
        time.sleep(1)
        usrp_locked = usrp.get_mboard_sensor('ref_locked').to_bool()
        if not usrp_locked:
            if i == (max_n_tries-1):
                raise RuntimeError('Failed to lock USRP to 10 MHz reference')
            else:
                print('USRP is not locked to 10 MHz reference, trying again...')
        else:
            print('USRP is locked to 10 MHz reference')
            break

    # Set the time to 0 on the next PPS edge
    t = usrp.get_time_last_pps()
    while(t == usrp.get_time_last_pps()):
        time.sleep(0.1)
    time.sleep(0.1)
    usrp.set_time_next_pps(uhd.time_spec(0))
    t = usrp.get_time_last_pps()
    while(t == usrp.get_time_last_pps()):
        time.sleep(0.1)
    time.sleep(0.1)
    t = usrp.get_time_last_pps().get_real_secs()
    if t != 0:
        raise RuntimeError('Failed to set USRP time to 0 on the PPS edge')
    else:
        print('Set USRP time to 0 on the PPS edge')

class node(gr.top_block):

    """
    Flowgraph that instantiates the hardware and signal processing for each type
    of beamforming node.
    """

    def __init__(self, node_id):

        gr.top_block.__init__(self, "node")

        # ----------------------------------------------------------------------
        # General setup
        # ----------------------------------------------------------------------

        # Sample rate (SPS)
        samp_rate = 1e6

        # Center frequency (Hz)
        freq = 915e6

        # TX/RX gains (dB)
        # NOTE: make sure to use attenuators in a wired setup!
        tx_gain = 57
        rx_gain = 53
        
        """
        # Can set serial numbers for each node if needed (required if all nodes
        # are attached to the same computer)
        # Example: dev_args += 'serial=1234567'
        dev_args = ''
        if node_id == 0:
            dev_args += ''
        elif node_id == 1:
            dev_args += ''
        elif node_id == 2:
            dev_args += ''
        else:
            print('WARNING: no serial number configured for this node number')
        """
        
        # Blocks for debugging / saving data
        enable_tx_tag_debug = 0
        enable_rx_iq_file_sink = 0
        enable_rx_msg_debug = 0
        enable_rx_msg_file_sink = 0

        # USRP stream arguments
        stream_args = uhd.stream_args(cpu_format="fc32",
                                      args='',
                                      channels=[0])

        # USRP tune request - no digital tuning (to prevent frequency offset)
        tune_req = uhd.tune_request_t(freq)
        tune_req.dsp_freq_policy = uhd.tune_request_t.POLICY_MANUAL
        tune_req.dsp_freq = 0

        # ----------------------------------------------------------------------
        # Transmitter blocks
        # ----------------------------------------------------------------------

        # USRP transmitter
        self.tx_usrp = uhd.usrp_sink(dev_args,
                                     stream_args,
                                     'packet_len')
        # Specify clock rate during TX initialization only
        self.tx_usrp.set_clock_rate(32e6, 0)
        self.tx_usrp.set_subdev_spec('A:A', 0)
        self.tx_usrp.set_center_freq(tune_req, 0)
        self.tx_usrp.set_gain(tx_gain, 0)
        self.tx_usrp.set_antenna('TX/RX', 0)
        self.tx_usrp.set_samp_rate(samp_rate)

        # Frequency and phase offset adjustment (only needed for the leader and
        # followers)
        if node_id != 0:
            self.tx_freq_phase_shift = freq_phase_shift(samp_rate, node_id)

        # PDU to tagged stream
        self.tx_pdu_to_ts = blocks.pdu_to_tagged_stream(blocks.complex_t, 'packet_len')

        # Tag debugger
        if enable_tx_tag_debug:
            self.tx_tag_debug = blocks.tag_debug(gr.sizeof_gr_complex*1, '', "")
            self.tx_tag_debug.set_display(True)

        # Connections
        if node_id != 0:
            self.msg_connect((self.tx_freq_phase_shift, 'msg_out'), (self.tx_pdu_to_ts, 'pdus'))
        self.connect((self.tx_pdu_to_ts, 0), (self.tx_usrp, 0))
        if enable_tx_tag_debug:
             self.connect((self.tx_pdu_to_ts, 0), (self.tx_tag_debug, 0))

        # ----------------------------------------------------------------------
        # Receiver blocks
        # ----------------------------------------------------------------------

        # USRP receiver (doesn't start streaming immediately)
        rx_stream_on_start = False
        self.rx_usrp = uhd.usrp_source(dev_args,
                                       stream_args,
                                       rx_stream_on_start)
        self.rx_usrp.set_subdev_spec('A:A', 0)
        self.rx_usrp.set_center_freq(tune_req, 0)
        self.rx_usrp.set_gain(rx_gain, 0)
        self.rx_usrp.set_antenna('TX/RX', 0)
        self.rx_usrp.set_samp_rate(samp_rate)

        # Disable AGC
        self.rx_usrp.set_rx_agc(False, 0)

        # NOTE: DC offset correction and IQ imbalance correction are left on.
        # They don't seem to affect phase or frequency synchronization. It looks
        # like they actually help quite a bit for beamforming.
        self.rx_usrp.set_auto_dc_offset(False, 0)
        #self.rx_usrp.set_auto_iq_balance(False, 0)

        # NOTE: for some reason if auto IQ balance is turned off at the RX, a
        # frequency offset appears between the TX and RX (even when sync'd to
        # the same reference). The offset will only go away if the TX tune
        # request is sent again. Not sure why this happens, but need to keep
        # it in mind if we ever disable the IQ imbalance correction. Based on
        # separate testing in GRC, this can also happen with the RX tune request
        # (so to be safe, resend both TX/RX tunes afterwards).
        #self.tx_usrp.set_center_freq(tune_req, 0)
        #self.rx_usrp.set_center_freq(tune_req, 0)

        # Set the number of output items to be the receive packet size (1.25 ms of data)
        self.rx_usrp.set_min_noutput_items(0.00125*samp_rate)
        self.rx_usrp.set_max_noutput_items(0.00125*samp_rate)

        # Make the receive timeout large so we block here instead of between
        # calls to work by the scheduler. This way when data arrives we can
        # return immediately with a tag, and the next time data arrives it will
        # also be tagged (since only one tag per call to work).
        self.rx_usrp.set_recv_timeout(1.0, True)

        # Synchronize to our external reference after the device is done being
        # initialized (otherwise time between TX/RX will not be in sync)
        sync_usrp_to_external(self.rx_usrp)

        # File sink with metadata for IQ samples
        if enable_rx_iq_file_sink:
            ts = datetime.datetime.now().strftime('%m-%d-%Y_%H.%M.%S')
            self.rx_iq_file_sink = blocks.file_meta_sink(gr.sizeof_gr_complex*1,
                                                         'data/rx_iq_%s_%s.dat' % (node_id, ts),
                                                         samp_rate,
                                                         1,
                                                         blocks.GR_FILE_FLOAT,
                                                         True,
                                                         2**32,
                                                         pmt.make_dict(),
                                                         True)
            self.rx_iq_file_sink.set_unbuffered(False)
        else:
            self.rx_null_sink = blocks.null_sink(gr.sizeof_gr_complex*1)

        # Packet detector
        # Power threshold set based on noise with RX gain of 70 dB @ 915 MHz
        # Minimum/maximum packet durations are 0.25 and 1.25 ms
        # If the signal is low for more than 0.025 ms, it's dropped
        pow_thresh = 10e-3
        t_min1 = 0.025e-3
        t_min2 = 0.25e-3
        t_max = 1.25e-3
        self.rx_pkt_detector = pkt_detector(samp_rate, pow_thresh, t_min1, t_min2, t_max)

        # Packet estimator
        self.rx_pkt_estimator = pkt_estimator(samp_rate, node_id)

        # Message debugger
        if enable_rx_msg_debug:
            self.rx_msg_debug = blocks.message_debug()

        # File sink with metadata for detected packets
        # This requires conversion from PDU to a tagged stream
        if enable_rx_msg_file_sink:
            ts = datetime.datetime.now().strftime('%m-%d-%Y_%H.%M.%S')
            self.rx_pdu_to_ts = blocks.pdu_to_tagged_stream(blocks.complex_t,
                                                            'packet_len')
            self.rx_msg_file_sink = blocks.file_meta_sink(gr.sizeof_gr_complex*1,
                                                          'data/rx_msg_%s_%s.dat' % (node_id, ts),
                                                          samp_rate,
                                                          1,
                                                          blocks.GR_FILE_FLOAT,
                                                          True,
                                                          2**32,
                                                          pmt.make_dict(),
                                                          True)
            self.rx_msg_file_sink.set_unbuffered(False)

        # Connections
        self.connect((self.rx_usrp, 0), (self.rx_pkt_detector, 0))
        self.msg_connect((self.rx_pkt_detector, 'msg_out'), (self.rx_pkt_estimator, 'msg_in'))
        if enable_rx_iq_file_sink:
            self.connect((self.rx_usrp, 0), (self.rx_iq_file_sink, 0))
        else:
            self.connect((self.rx_usrp, 0), (self.rx_null_sink, 0))
        if enable_rx_msg_file_sink:
            self.msg_connect((self.rx_pkt_detector, 'msg_out'), (self.rx_pdu_to_ts, 'pdus'))
            self.connect((self.rx_pdu_to_ts, 0), (self.rx_msg_file_sink, 0))
        if enable_rx_msg_debug:
            self.msg_connect((self.rx_pkt_estimator, 'msg_out'), (self.rx_msg_debug, 'print'))

        # ----------------------------------------------------------------------
        # Control blocks
        # ----------------------------------------------------------------------

        # Controller block
        self.usrp_ctrl = usrp_ctrl(node_id, self.rx_usrp)

        # Connections
        if node_id != 0:
            self.msg_connect((self.usrp_ctrl, 'msg_out'), (self.tx_freq_phase_shift, 'msg_in'))
        else:
            self.msg_connect((self.usrp_ctrl, 'msg_out'), (self.tx_pdu_to_ts, 'pdus'))
        self.msg_connect((self.rx_pkt_estimator, 'msg_out'), (self.usrp_ctrl, 'msg_in'))

def main():

    # Main function - creates flowgraph and runs it

    # Get node ID from command line
    # Target: 0
    # Leader: 1
    # Follower: >= 2
    n_arg = len(sys.argv)
    if n_arg == 2:
        node_id = int(sys.argv[1])
        print('Node ID: %i' % node_id)
    else:
        raise RuntimeError('Program requires one argument (node ID)')

    # Make output data directory if needed
    if not os.path.exists('data'):
        os.mkdir('data')

    # Create the flowgraph
    tb = node(node_id)

    # Setup handler to stop on CTRL+C
    def sig_handler(sig=None, frame=None):
        tb.stop()
        tb.wait()
        sys.exit(0)
    signal.signal(signal.SIGINT, sig_handler)
    signal.signal(signal.SIGTERM, sig_handler)

    # Start the flowgraph
    tb.start()

    # Prompt to quit
    try:
        input('Press Enter to quit:\n')
    except EOFError:
        pass

    # Stop the flowgraph
    tb.stop()
    tb.wait()

    # Convert RX metadata to MATLAB format
    if os.path.exists('data/rx_iq_%s.dat.hdr' % node_id):
        hdr_to_mat('data/rx_iq_%s.dat.hdr' % node_id)
    if os.path.exists('data/rx_msg_%s.dat.hdr' % node_id):
        hdr_to_mat('data/rx_msg_%s.dat.hdr' % node_id)

if __name__ == '__main__':
    main()

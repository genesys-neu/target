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

# Some common utilities

import time
import datetime
import calendar
import sys

from gnuradio import uhd

import octoclock

# Custom style sheet for QT applications
qss = """
QwtPlotCanvas
{
    background: black;
}

DisplayPlot {
    qproperty-zoomer_color: white;
    qproperty-line_color1: yellow;
    qproperty-line_color2: #00FF00;
}

TimeDomainDisplayPlot {
    qproperty-tag_text_color: black;
    qproperty-tag_background_color: yellow;
    qproperty-tag_background_style: SolidPattern;
}

FrequencyDisplayPlot {
    qproperty-marker_lower_intensity_visible: false;
    qproperty-marker_upper_intensity_visible: false;
    qproperty-marker_noise_floor_amplitude_visible: false;
}
"""

def sync_to_internal(usrp):

    # Synchronizes a USRP to its internal references

    # Set USRP clock and time sources
    usrp.set_clock_source('internal', 0)
    usrp.set_time_source('internal', 0)
    
    # Set time from PC on the next PPS edge
    set_time_from_pc(usrp)

def sync_to_gpsdo(usrp, sync_time_to_gps):

    # Synchronizes a USRP to an onboard GPSDO

    # If GPS is required, make sure we have lock
    gps_locked = usrp.get_mboard_sensor('gps_locked').to_bool()
    if gps_locked:
        print('GPS is locked')
    else:
        print('GPS is NOT locked')
        if sync_time_to_gps:
            raise RuntimeError('GPS is NOT locked and GPS time is required')

    # Set USRP clock and time sources
    usrp.set_clock_source('gpsdo', 0)
    usrp.set_time_source('gpsdo', 0)

    # Make sure the USRP is locked to the 10 MHz reference
    time.sleep(1)
    while(not usrp.get_mboard_sensor('ref_locked').to_bool()):
        print('USRP 10 MHz reference is NOT locked, trying again in 1 second...')
        time.sleep(1)
    print('USRP 10 MHz reference is locked')

    # Set time from GPS on the next PPS edge
    if sync_time_to_gps:

        # Wait for the next PPS, then wait 200 ms for the NMEA string to
        # propagate, then get GPS time, and finally set the USRP time to GPS
        # time +1 second on the next PPS.
        t = usrp.get_time_last_pps()
        while(t == usrp.get_time_last_pps()):
            time.sleep(0.1)
        time.sleep(0.2)
        gps_time = usrp.get_mboard_sensor('gps_time').to_int()
        usrp.set_time_next_pps(uhd.time_spec(gps_time + 1)) 

        # Wait a bit and make sure the USRP time matches GPS on the next PPS
        time.sleep(1.1)
        t = usrp.get_time_last_pps()
        while(t == usrp.get_time_last_pps()):
            time.sleep(0.1)
        time.sleep(0.2)
        gps_time = usrp.get_mboard_sensor('gps_time').to_int()
        usrp_time = usrp.get_time_last_pps().get_real_secs()
        print('USRP time: %s' % usrp_time)
        print('GPS time : %s' % gps_time)
        if usrp_time != gps_time:
            raise RuntimeError('USRP time is NOT locked to GPS')
        print('USRP time is locked to GPS')

    # Set time from PC on the next PPS edge
    else:
        set_time_from_pc(usrp)

def sync_to_octoclock_g(usrp, sync_time_to_gps):

    # Synchronizes a USRP to the OctoClock-G

    # If we need to synchronize to GPS time
    if sync_time_to_gps:

        # Create an OctoClock
        clk = octoclock.octoclock('')

        # Make sure we have GPS lock
        gps_locked = clk.get_gps_locked()
        if gps_locked:
            print('OctoClock-G GPS is locked')
        else:
            print('OctoClock-G GPS is NOT locked')
            raise RuntimeError('OctoClock-G GPS is NOT locked and GPS time is required')

    # Set USRP clock and time sources
    usrp.set_clock_source('external', 0)
    usrp.set_time_source('external', 0)

    # Make sure the USRP is locked to the 10 MHz reference
    time.sleep(1)
    while(not usrp.get_mboard_sensor('ref_locked').to_bool()):
        print('USRP 10 MHz reference is NOT locked, trying again in 1 second...')
        time.sleep(1)
    print('USRP 10 MHz reference is locked')

    # Set time from GPS on the next PPS edge
    if sync_time_to_gps:

        # Wait for the next PPS, then wait 200 ms for the NMEA string to
        # propagate, then get GPS time, and finally set the USRP time to GPS
        # time +1 second on the next PPS.
        t = usrp.get_time_last_pps()
        while(t == usrp.get_time_last_pps()):
            time.sleep(0.1)
        time.sleep(0.2)
        gps_time = clk.get_gps_time()
        usrp.set_time_next_pps(uhd.time_spec(gps_time + 1)) 

        # Wait a bit and make sure the USRP time matches GPS on the next PPS
        time.sleep(1.1)
        t = usrp.get_time_last_pps()
        while(t == usrp.get_time_last_pps()):
            time.sleep(0.1)
        time.sleep(0.2)
        gps_time = clk.get_gps_time()
        usrp_time = usrp.get_time_last_pps().get_real_secs()
        print('USRP time: %s' % usrp_time)
        print('GPS time : %s' % gps_time)
        if usrp_time != gps_time:
            raise RuntimeError('USRP time is NOT locked to GPS')
        print('USRP time is locked to GPS')

    # Set time from PC on the next PPS edge
    else:
        set_time_from_pc(usrp)

def set_time_from_pc(usrp):

    # Wait for the PC time to change seconds, then get the PC time and set the
    # USRP time to PC time + 1 on the next PPS.
    t = int(time.time())
    while(t == int(time.time())):
        time.sleep(0.1)
    time.sleep(0.2)
    pc_time = datetime.datetime.utcnow().replace(microsecond=0)
    pc_time = float(calendar.timegm(pc_time.timetuple()))
    usrp.set_time_next_pps(uhd.time_spec(pc_time + 1))

    # Wait a bit and make sure the USRP time roughly matches the PC time on the
    # next PPS
    time.sleep(1.1)
    t = int(time.time())
    while(t == int(time.time())):
        time.sleep(0.1)
    time.sleep(0.2)
    pc_time = datetime.datetime.utcnow().replace(microsecond=0)
    pc_time = float(calendar.timegm(pc_time.timetuple()))
    usrp_time = usrp.get_time_last_pps().get_real_secs()
    print('USRP time: %s' % usrp_time)
    print('PC time  : %s' % pc_time)
    if usrp_time != pc_time:
        raise RuntimeError('USRP time is NOT locked to PC time')
    print('USRP time is locked to PC time')

    # Make sure console is up to date
    sys.stdout.flush()

def tune_without_digital(usrp, freq):

    # Tunes without using the CORDIC

    tune_req = uhd.tune_request_t(freq)
    tune_req.dsp_freq_policy = uhd.tune_request_t.POLICY_MANUAL
    tune_req.dsp_freq = 0
    usrp.set_center_freq(tune_req, 0)

def get_gps_start_t(usrp):

    # Gets start time from GPS (start of the next minute)

    t_now = usrp.get_time_now().get_real_secs()
    dt_now = datetime.datetime.utcfromtimestamp(t_now)
    start_dt = dt_now + datetime.timedelta(minutes=1)
    start_dt = start_dt.replace(second=0, microsecond=0)
    start_t = calendar.timegm(start_dt.timetuple())

    return start_t

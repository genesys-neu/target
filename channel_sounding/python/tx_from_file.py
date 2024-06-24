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

# ------------------------------------------------------------------------------
# Imports and setup
# ------------------------------------------------------------------------------

import sys
import ctypes
import signal
import time
import datetime
import math

from gnuradio import gr, eng_notation, qtgui, blocks, uhd
from gnuradio.filter import firdes
import pmt
import sip

import config
import util

# Get versions of Python and GNU Radio
py_version = int(sys.version[0])
gr_version = float(gr.version()[:3])

if gr_version >= 3.8:
    from PyQt5 import Qt
else:
    from PyQt4 import Qt

if py_version >= 3:
    from packaging import version
else:
    from distutils.version import StrictVersion

if sys.platform.startswith('linux'):
    try:
        x11 = ctypes.cdll.LoadLibrary('libX11.so')
        x11.XInitThreads()
    except:
        print("Warning: failed to XInitThreads()")

if py_version >= 3:
    if version.Version("4.5.0") <= version.Version(Qt.qVersion()) < version.Version("5.0.0"):
        style = gr.prefs().get_string('qtgui', 'style', 'raster')
        Qt.QApplication.setGraphicsSystem(style)
else:
    if StrictVersion("4.5.0") <= StrictVersion(Qt.qVersion()) < StrictVersion("5.0.0"):
        style = gr.prefs().get_string('qtgui', 'style', 'raster')
        Qt.QApplication.setGraphicsSystem(style)

# ------------------------------------------------------------------------------
# Parse arguments
# ------------------------------------------------------------------------------

# Not in debug mode by default
# Debug mode allows changing the gain and frequency
DEBUG = False
n_arg = len(sys.argv)
if (n_arg == 2) and (int(sys.argv[1]) == 1):
    DEBUG = True
    print('Running in debug mode')

# ------------------------------------------------------------------------------
# Flowgraph
# ------------------------------------------------------------------------------

class tx_from_file(gr.top_block, Qt.QWidget):

    def __init__(self):

        # Set title
        self.win_title = "TX From File"
        if DEBUG:
            self.win_title = "TX From File (DEBUG)"

        # Initialize the flowgraph and QT GUI
        gr.top_block.__init__(self, self.win_title)
        Qt.QWidget.__init__(self)
        self.setWindowTitle(self.win_title)
        qtgui.util.check_set_qss()
        try:
            self.setWindowIcon(Qt.QIcon.fromTheme('gnuradio-grc'))
        except:
            pass
        self.top_scroll_layout = Qt.QVBoxLayout()
        self.setLayout(self.top_scroll_layout)
        self.top_scroll = Qt.QScrollArea()
        self.top_scroll.setFrameStyle(Qt.QFrame.NoFrame)
        self.top_scroll_layout.addWidget(self.top_scroll)
        self.top_scroll.setWidgetResizable(True)
        self.top_widget = Qt.QWidget()
        self.top_scroll.setWidget(self.top_widget)
        self.top_layout = Qt.QVBoxLayout(self.top_widget)
        self.top_grid_layout = Qt.QGridLayout()
        self.top_layout.addLayout(self.top_grid_layout)
        self.settings = Qt.QSettings("GNU Radio", "tx_from_file")
        try:
            if py_version >= 3:
                if version.Version(Qt.qVersion()) < version.Version("5.0.0"):
                    self.restoreGeometry(self.settings.value("geometry").toByteArray())
                else:
                    self.restoreGeometry(self.settings.value("geometry"))
            else:
                if StrictVersion(Qt.qVersion()) < StrictVersion("5.0.0"):
                    self.restoreGeometry(self.settings.value("geometry").toByteArray())
                else:
                    self.restoreGeometry(self.settings.value("geometry"))
        except:
            pass

        # ----------------------------------------------------------------------
        # Variables
        # ----------------------------------------------------------------------

        self.dev_args = config.tx_dev_args
        self.samp_rate = config.samp_rate
        self.gain = config.tx_gain
        self.freq = config.freq
        self.sync_src = config.sync_src
        self.sync_time_to_gps = config.sync_time_to_gps
        self.tx_file = config.tx_file
        self.plots_on = config.tx_plots_on

        print('Configuration:')
        print('Device args: %s' % self.dev_args)
        print('Sample rate: %s SPS' % self.samp_rate)
        print('Gain: %s dB' % self.gain)
        print('Frequency: %s Hz' % self.freq)
        print('Sync source: %s' % self.sync_src)
        print('Sync time to GPS: %s' % self.sync_time_to_gps)
        print('TX file: %s' % self.tx_file)
        print('Plots on: %s' % self.plots_on)

        # ----------------------------------------------------------------------
        # Blocks
        # ----------------------------------------------------------------------

        # Gain entry
        self.gain_entry = Qt.QToolBar(self)
        self.gain_entry.addWidget(Qt.QLabel('Gain (dB)' + ": "))
        self.gain_edit = Qt.QLineEdit(str(self.gain))
        self.gain_entry.addWidget(self.gain_edit)
        self.gain_edit.returnPressed.connect(
            lambda: self.set_gain(eng_notation.str_to_num(str(self.gain_edit.text()))))
        self.top_layout.addWidget(self.gain_entry)

        # Frequency entry
        self.freq_entry = Qt.QToolBar(self)
        self.freq_entry.addWidget(Qt.QLabel('Frequency (Hz)' + ": "))
        self.freq_edit = Qt.QLineEdit(str(self.freq))
        self.freq_entry.addWidget(self.freq_edit)
        self.freq_edit.returnPressed.connect(
            lambda: self.set_freq(eng_notation.str_to_num(str(self.freq_edit.text()))))
        self.top_layout.addWidget(self.freq_entry)

        # Disable changing gain/frequency if we're not in debug mode
        if not DEBUG:
            self.gain_edit.setEnabled(False)
            self.freq_edit.setEnabled(False)

        # USRP transmitter
        self.usrp_sink = uhd.usrp_sink(
            self.dev_args,
            uhd.stream_args(
                cpu_format="sc16",
                args='',
                channels=list(range(0,1)),
            ),
            '',
        )
        self.usrp_sink.set_subdev_spec('A:A', 0)
        util.tune_without_digital(self.usrp_sink, self.freq)
        self.usrp_sink.set_gain(self.gain, 0)
        self.usrp_sink.set_antenna('TX/RX', 0)
        self.usrp_sink.set_samp_rate(self.samp_rate)

        # Synchronize to 10 MHz and PPS references, and set time
        if self.sync_src == 'internal':
            if self.sync_time_to_gps:
                raise RuntimeError('Cannot sync time to GPS when using internal references')
            util.sync_to_internal(self.usrp_sink)
        elif self.sync_src == 'gpsdo':
            util.sync_to_gpsdo(self.usrp_sink, self.sync_time_to_gps)
        elif self.sync_src == 'octoclock-g':
            util.sync_to_octoclock_g(self.usrp_sink, self.sync_time_to_gps)
        else:
            raise RuntimeError('Unknown sync source: %s', self.sync_src)

        # Only show the frequency sink if plots are on
        if self.plots_on:

            # QT frequency sink GUI
            self.qt_freq_sink = qtgui.freq_sink_c(
                1024, #size
                firdes.WIN_BLACKMAN_hARRIS, #wintype
                0, #fc
                self.samp_rate, #bw
                "", #name
                1
            )
            self.qt_freq_sink.set_update_time(0.10)
            self.qt_freq_sink.set_y_axis(-140, 10)
            self.qt_freq_sink.set_y_label('Relative Gain', 'dB')
            self.qt_freq_sink.set_trigger_mode(qtgui.TRIG_MODE_FREE, 0.0, 0, "")
            self.qt_freq_sink.enable_autoscale(False)
            self.qt_freq_sink.enable_grid(True)
            self.qt_freq_sink.set_fft_average(1.0)
            self.qt_freq_sink.enable_axis_labels(True)
            self.qt_freq_sink.enable_control_panel(False)
            self.qt_freq_sink.disable_legend()
            self.qt_freq_sink_win = sip.wrapinstance(self.qt_freq_sink.pyqwidget(), Qt.QWidget)
            self.top_layout.addWidget(self.qt_freq_sink_win)

            # Convert full scale to [-1, +1]
            self.multiply_const = blocks.multiply_const_cc(1.0/(2**15))

            # Complex int16 (sc16) to complex floating point
            self.sc16_to_complex = blocks.interleaved_short_to_complex(True, False)

        # File source playing out sc16 samples
        if gr_version >= 3.8:
            self.file_source = blocks.file_source(gr.sizeof_short*2, self.tx_file, True, 0, 0)
        else:
            self.file_source = blocks.file_source(gr.sizeof_short*2, self.tx_file, True)
        self.file_source.set_begin_tag(pmt.PMT_NIL)

        # ----------------------------------------------------------------------
        # Connections
        # ----------------------------------------------------------------------

        # File input to USRP TX
        self.connect((self.file_source, 0), (self.usrp_sink, 0))

        # View on frequency sink for debugging
        if self.plots_on:
            self.connect((self.file_source, 0), (self.sc16_to_complex, 0))
            self.connect((self.sc16_to_complex, 0), (self.multiply_const, 0))
            self.connect((self.multiply_const, 0), (self.qt_freq_sink, 0))

        # --------------------------------------------------------------------------
        # Wait to start
        # --------------------------------------------------------------------------

        # Make sure console is up to date
        print('Prompting for input (press enter when ready to start)')
        sys.stdout.flush()

        if py_version == 3:
            input('Press enter when ready to start')
        else:
            raw_input('Press enter when ready to start')

        # Get the start time. When using GPS, this is the start of the next
        # minute. Otherwise, it is 5 seconds from now.
        if self.sync_time_to_gps:
            start_t = util.get_gps_start_t(self.usrp_sink)
        else:
            start_t = math.ceil(self.usrp_sink.get_time_now().get_real_secs()) + 5

        # Set the start time and print status
        self.usrp_sink.set_start_time(uhd.time_spec(start_t))
        t_now = self.usrp_sink.get_time_now().get_real_secs()
        t_wait = int(start_t - t_now)
        start_dt = datetime.datetime.utcfromtimestamp(start_t)
        print('USRP will start in %s seconds at %s (UTC)' % (t_wait, start_dt))

        # Make sure console is up to date
        sys.stdout.flush()

    # --------------------------------------------------------------------------
    # Methods
    # --------------------------------------------------------------------------

    def closeEvent(self, event):
        self.settings = Qt.QSettings("GNU Radio", "tx_from_file")
        self.settings.setValue("geometry", self.saveGeometry())
        event.accept()

    def set_gain(self, gain):
        self.gain = gain
        Qt.QMetaObject.invokeMethod(self.gain_edit, "setText", Qt.Q_ARG("QString", eng_notation.num_to_str(self.gain)))
        self.usrp_sink.set_gain(self.gain, 0)

    def set_freq(self, freq):
        self.freq = freq
        Qt.QMetaObject.invokeMethod(self.freq_edit, "setText", Qt.Q_ARG("QString", eng_notation.num_to_str(self.freq)))
        util.tune_without_digital(self.usrp_sink, self.freq)

# ------------------------------------------------------------------------------
# Main program
# ------------------------------------------------------------------------------

if __name__ == "__main__":

    # Starting message
    print('Started at %s (UTC)' % datetime.datetime.utcnow())

    # Create the QT application
    qapp = Qt.QApplication(sys.argv)
    qapp.setStyleSheet(util.qss)

    # Create and start the flowgraph
    tx = tx_from_file()
    tx.start()
    tx.show()

    # Function to quit the QT application
    def sig_handler(sig=None, frame=None):
        Qt.QApplication.quit()
    signal.signal(signal.SIGINT, sig_handler)
    signal.signal(signal.SIGTERM, sig_handler)

    # Function to quit the flowgraph
    def quitting():
        tx.stop()
        tx.wait()
    qapp.aboutToQuit.connect(quitting)

    # Start the QT application
    qapp.exec_()

    # Stopping message
    print('Stopped at %s (UTC)' % datetime.datetime.utcnow())

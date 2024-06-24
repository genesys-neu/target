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

import time
import sys
import signal

from gnuradio import gr, uhd
import numpy as np
from scipy.io import savemat

from tx_from_file import tx_from_file
from rx_to_file import rx_to_file
import util
import config

# Get version of GNU Radio
gr_version = float(gr.version()[:3])

if gr_version >= 3.8:
    from PyQt5 import Qt
else:
    from PyQt4 import Qt

# Create the QT application
qapp = Qt.QApplication(sys.argv)
qapp.setStyleSheet(util.qss)

# Function to quit the QT application
def sig_handler(sig=None, frame=None):
	Qt.QApplication.quit()
signal.signal(signal.SIGINT, sig_handler)
signal.signal(signal.SIGTERM, sig_handler)

# Create the two USRPs. By default these use the sync source as specified in the
# config file for both time and frequency references.
tx = tx_from_file()
rx = rx_to_file()

# Change both USRPs to use an external time reference. This will allow us to
# trigger a measurement at the same time but have each USRP use their own 10 MHz
# frequency references for counting clock ticks.
tx.usrp_sink.set_time_source('external', 0)
rx.usrp_source.set_time_source('external', 0)

# Wait a few seconds for setup
time.sleep(3)

# Helper function for resetting the time on each USRP
def reset_times(tx, rx):

	# Sync the USRPs to the same time on the next PPS edge
	t = tx.usrp_sink.get_time_last_pps()
	while(t == tx.usrp_sink.get_time_last_pps()):
		time.sleep(0.1)
	time.sleep(0.2)
	tx.usrp_sink.set_time_next_pps(uhd.time_spec(0))
	rx.usrp_source.set_time_next_pps(uhd.time_spec(0))

	# Wait for the next PPS
	t = tx.usrp_sink.get_time_last_pps()
	while(t == tx.usrp_sink.get_time_last_pps()):
		time.sleep(0.1)
	time.sleep(0.2)

# Number of measurements
n = 1000

# Collect data
tx_ticks = np.empty(n)
rx_ticks = np.empty(n)
tx_times = np.empty(n)
rx_times = np.empty(n)
for i in range(n):

	# Reset the time
	reset_times(tx, rx)

	# Wait for the next PPS and record the tick counts and times
	t = tx.usrp_sink.get_time_last_pps()
	while(t == tx.usrp_sink.get_time_last_pps()):
		time.sleep(0.1)
	time.sleep(0.2)
	tx_time = tx.usrp_sink.get_time_last_pps()
	rx_time = rx.usrp_source.get_time_last_pps()
	tx_ticks[i] = tx_time.to_ticks(61.44e6)
	rx_ticks[i] = rx_time.to_ticks(61.44e6)
	tx_times[i] = tx_time.get_real_secs()
	rx_times[i] = rx_time.get_real_secs()

	# Print status
	print('Measurement %i of %i' % (i, n))
	print('TX USRP ticks: %s' % tx_ticks[i])
	print('RX USRP ticks: %s' % rx_ticks[i])
	print('TX USRP time: %s' % tx_times[i])
	print('RX USRP time: %s\n' % rx_times[i])

# Save to .mat file
data_out = {'tx_ticks': [], 'rx_ticks': [], 'tx_times': [], 'rx_times': []}
data_out['tx_ticks'] = list(tx_ticks)
data_out['rx_ticks'] = list(rx_ticks)
data_out['tx_times'] = list(tx_times)
data_out['rx_times'] = list(rx_times)
file_out = '../data/' + config.sync_src + '.mat'
savemat(file_out, data_out)
print('Saved file: %s' % file_out)

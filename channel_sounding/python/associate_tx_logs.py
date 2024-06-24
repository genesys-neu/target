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

# This goes through the TX log files and puts them in the directories with the
# corresponding RX data. Note that the timestamp in the TX log file names cannot
# be used, since the TX2 clock resets to some time around February 2016 each time
# the system was rebooted. But the file contents have timestamps from GPS, which
# is what we look for here.

from os import listdir
from os.path import isfile, isdir, join
import datetime
import calendar

# Directories with the original TX logs (text files) and RX data (folders)
tx_log_dir = '../tx_logs'
rx_data_dir = '../data'

# All TX log files
tx_logs = [f for f in listdir(tx_log_dir) if isfile(join(tx_log_dir, f))]

# All RX data folders
rx_folders = [f for f in listdir(rx_data_dir) if isdir(join(rx_data_dir, f))]

# Associate each of the log files
for i in tx_logs:

	# Make sure it's a TX log file
	if i[:6] != 'tx_log':
		continue

	# Read in the data
	file_in = tx_log_dir + '/' + i
	print(file_in)
	fid = open(file_in, 'r')
	file_data = fid.read()
	fid.close()

	# Look for the timestamp when the USRP will start
	idx = file_data.find('seconds at')
	if idx < 0:
		print('Could not find timestamp in %s, skipping' % i)
		continue
	start_str = file_data[idx+11:idx+30]
	start_dt = datetime.datetime.strptime(start_str, '%Y-%m-%d %H:%M:%S')
	start_t = calendar.timegm(start_dt.timetuple())

	# See if we have a corresponding RX data folder
	has_rx = sum([j == str(start_t) for j in rx_folders])
	if has_rx == 0:
		print('No corresponding RX folder for %s, skipping' % i)
		continue
	if has_rx > 1:
		print('More than one matching RX folder for %s, skipping' % i)
		continue

	# Make sure the folder doesn't already have a TX log file in it
	rx_files = listdir(rx_data_dir + '/' + str(start_t))
	has_tx_log = sum([(j.find('tx_log') >= 0) for j in rx_files])
	if has_tx_log:
		print('RX folder for %s already has a TX log file, skipping' % start_t)
		continue

	# Write the file contents to the RX data directory
	print('Writing %s to RX directory %s' % (i, start_t))
	file_out = rx_data_dir + '/' + str(start_t) + '/' + i
	fid = open(file_out, 'w')
	fid.write(file_data)
	fid.close()

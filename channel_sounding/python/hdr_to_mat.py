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

# Convert .hdr metadata file to .mat for processing

import sys

from scipy.io import savemat
from gnuradio.blocks import parse_file_metadata
import pmt

# Input file is given as an argument to the program
n_arg = len(sys.argv)
if n_arg < 2:
	raise RuntimeError('Missing input file')
file_in = sys.argv[1]

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

	# Add to the output dict
	data_out['rx_freq'].append(pmt.to_double(info['rx_freq']))
	data_out['rx_rate'].append(info['rx_rate'])
	data_out['rx_time'].append(info['rx_time'])
	data_out['n_items'].append(info['nitems'])
	data_out['n_bytes'].append(info['nbytes'])

	# Move on to the next header
	n_headers += 1
	n_read += parse_file_metadata.HEADER_LENGTH + info["extra_len"]
	fid.seek(n_read, 0)

# Save output file
file_out = file_in + '.mat'
savemat(file_out, data_out)

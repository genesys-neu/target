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

# Configuration for the TX/RX programs

# Device arguments
tx_dev_args = ''
rx_dev_args = ''

# Sample rate (SPS)
samp_rate = 1e6

# Carrier frequency (Hz)
freq = 915e6

# RF gain (dB)
tx_gain = 0
rx_gain = 0

# Sync source (for 10 MHz and PPS)
# Options are internal, gpsdo (onboard), or octoclock-g
#sync_src = 'internal'
#sync_src = 'gpsdo'
sync_src = 'octoclock-g'

# Sync time to GPS
# If this is off, we sync to the PC clock and start in 5 seconds
# If it's on, we sync to GPS time and start at the next minute
# (GPS time sync requires a gpsdo or octoclock-g)
sync_time_to_gps = 1

# TX playout file
tx_file = '../matlab/tone.dat'
#tx_file = '../matlab/zc.dat'

# Plots on/off
tx_plots_on = 0
rx_plots_on = 1
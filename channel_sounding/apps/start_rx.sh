#!/bin/bash

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

# Runs the RX program

# Get timestamp
TIME=$(date -u +%s)

# Start the RX program with logging
python ../python/rx_to_file.py 2>&1 | tee ../data/rx_log_$TIME.txt

# Move the log to the most recently created folder
DIR=$(ls -td ../data/*/ | head -1)
mv ../data/rx_log_$TIME.txt ${DIR%?}

# Convert metadata in .hdr to .mat
python ../python/hdr_to_mat.py ../data/${DIR%?}/rx.dat.hdr

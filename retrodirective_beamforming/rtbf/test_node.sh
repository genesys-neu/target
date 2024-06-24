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

# Shell wrapper script for the main program. This sets the node ID, runs the
# program, and logs information to file for post-processing.

# Node ID is the only argument
NODE=$1

# Description of the node
if [ "$NODE" == "0" ]; then
	DESC="target"
elif [ "$NODE" == "1" ]; then
	DESC="leader"
elif [ "$NODE" == "2" ]; then
	DESC="follower"
else
	echo "Bad node ID"
	exit
fi

# Timestamp
TS=$(date +%m-%d-%Y_%H.%M.%S)

# Log file name
FILE=${DESC}_log_${TS}.txt

# Start the program, printing everything to stdout + file
python test_node.py $NODE 2>&1 | tee data/${FILE}

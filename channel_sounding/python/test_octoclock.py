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

import octoclock

args = ''
clock = octoclock.octoclock(args)

print('External reference detected: %s' % clock.get_ext_ref_detected())
print('GPS detected: %s' % clock.get_gps_detected())
print(clock.get_using_ref())
print(clock.get_switch_pos())
print(clock.get_gps_gpgga())
print(clock.get_gps_gprmc())
print('GPS time: %s' % clock.get_gps_time())
print('GPS locked: %s' % clock.get_gps_locked())
print(clock.get_gps_servo())
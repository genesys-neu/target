/*********************************************************************************
* DISTRIBUTION STATEMENT A. Approved for public release. Distribution is unlimited.
*
* This material is based upon work supported under Air Force Contract No. FA8702-15-D-0001. Any opinions, findings, conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the U.S. Air Force.
*
* (c) 2023 Massachusetts Institute of Technology.
*
* Subject to FAR52.227-11 Patent Rights - Ownership by the contractor (May 2014)
*
* The software/firmware is provided to you on an As-Is basis
*
* Delivered to the U.S. Government with Unlimited Rights, as defined in DFARS Part 252.227-7013 or 7014 (Feb 2014). Notwithstanding any copyright notice, U.S. Government rights in this work are defined by DFARS 252.227-7013 or DFARS 252.227-7014 as detailed above. Use of this work other than as specifically authorized by the U.S. Government may violate any copyrights that exist in this work.
*********************************************************************************/

#include <iostream>

#include "octoclock.h"

octoclock::octoclock(const std::string &args)
{
	std::cout << "Creating an OctoClock with args: " << args << std::endl;
	clock = uhd::usrp_clock::multi_usrp_clock::make(args);
	gps_detected = clock->get_sensor("gps_detected").to_bool();
}

octoclock::~octoclock()
{
}

bool octoclock::get_ext_ref_detected()
{
	return clock->get_sensor("ext_ref_detected").to_bool();
}

bool octoclock::get_gps_detected()
{
	return gps_detected;
}

std::string octoclock::get_using_ref()
{
	return clock->get_sensor("using_ref").to_pp_string();
}

std::string octoclock::get_switch_pos()
{
	return clock->get_sensor("switch_pos").to_pp_string();
}

std::string octoclock::get_gps_gpgga()
{
	if (gps_detected)
		return clock->get_sensor("gps_gpgga").to_pp_string();
	else
		return "GPS_GPGGA: No GPS detected.";
}

std::string octoclock::get_gps_gprmc()
{
	if (gps_detected)
		return clock->get_sensor("gps_gprmc").to_pp_string();
	else
		return "GPS_GPRMC: No GPS detected.";
}

double octoclock::get_gps_time()
{
	if (gps_detected)
		return clock->get_sensor("gps_time").to_real();
	else
		return -1;
}

bool octoclock::get_gps_locked()
{
	if (gps_detected)
		return clock->get_sensor("gps_locked").to_bool();
	else
		return false;
}

std::string octoclock::get_gps_servo()
{
	if (gps_detected)
		return clock->get_sensor("gps_servo").to_pp_string();
	else
		return "GPS_SERVO: No GPS detected.";
}
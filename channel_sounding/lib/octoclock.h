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

#include <uhd/usrp_clock/multi_usrp_clock.hpp>

class octoclock
{
public:
	octoclock(const std::string &args);
	~octoclock();

	uhd::usrp_clock::multi_usrp_clock::sptr clock;
	bool gps_detected;

	bool get_ext_ref_detected();
	bool get_gps_detected();
	std::string get_using_ref();
	std::string get_switch_pos();
	std::string get_gps_gpgga();
	std::string get_gps_gprmc();
	double get_gps_time();
	bool get_gps_locked();
	std::string get_gps_servo();
};
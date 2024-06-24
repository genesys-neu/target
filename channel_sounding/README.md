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

# Description

Software for channel sounding measurements between a UAV transmitter (with TX2 processor and B210 SDR) and ground receiver (with laptop processor and B210 SDR). There are separate programs for transmitting the test signals and receiving data to file for post-processing.

# Setup

First make the OctoClock module:

```
cd lib
make
```

Then create the test signals using the following scripts in the `matlab` folder:

```
gen_tone.m
gen_zc_seq.m
```

## Use

Modify the run-time configuration using `python/config.py` and run the main programs from the `apps` folder. Note that the programs with `_debug` in the name do not save any data. A quick look at the data can be obtained using the script `view_rx_data.m` in the `matlab` folder. The scripts `run_rx_tone.m` and `run_rx_zc_seq.m` will perform the post-processing and plot the channel magnitude and phase over time. These results can be saved for further processing and development of channel models.

There are some other utilities within the `python` folder, such as associating the TX log files with the RX data folders, testing clock drift between two SDRs, and testing to make sure the OctoClock module was properly built.

## Software Versions

**Laptop**: Python 3.9.9, GNU Radio 3.8.4.0, UHD 3.15.0.0-MacPorts-Release  
**TX2**: Python 2.7.12, GNU Radio 3.7.13.4, UHD 3.13.1.HEAD-0-gbbce3e45
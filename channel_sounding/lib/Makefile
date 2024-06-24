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

OS := $(shell uname -s)
ifeq ($(OS), Linux)
	UHD_INC=/usr/local/include
	UHD_LIB=/usr/local/lib
	PY_VER=2.7
	PY_INC=/usr/include/python$(PY_VER)
	PY_LIB=/usr/lib/python$(PY_VER)
	LIB_FLAG=-shared
endif
ifeq ($(OS), Darwin)
	UHD_INC=/opt/local/include
	UHD_LIB=/opt/local/lib
	PY_VER=3.9
	PY_INC=/opt/local/Library/Frameworks/Python.framework/Versions/$(PY_VER)/include/python$(PY_VER)
	PY_LIB=/opt/local/Library/Frameworks/Python.framework/Versions/$(PY_VER)/lib
	LIB_FLAG=-dynamiclib
endif

all:
	swig -c++ -python octoclock.i
	g++ -fPIC -c octoclock.cc -I $(UHD_INC)
	g++ -fPIC -c octoclock_wrap.cxx -I $(UHD_INC) -I $(PY_INC)
	g++ $(LIB_FLAG) octoclock.o octoclock_wrap.o -o _octoclock.so -L $(PY_LIB) -lpython$(PY_VER) -L $(UHD_LIB) -luhd
	mv _octoclock.so ../python/
	mv octoclock.py ../python/
clean:
	rm -f *.o
	rm -f octoclock_wrap.cxx
	rm -f ../python/*.so
	rm -f ../python/octoclock.py
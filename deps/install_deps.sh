#!/bin/bash
make clean
cd deps
wget http://sourceforge.net/projects/boost/files/boost/1.60.0/boost_1_60_0.tar.gz
if [ -f boost_1_60_0.tar.gz ]; then
tar xzf boost_1_60_0.tar.gz
cd boost_1_60_0
./bootstrap.sh --with-libraries=serialization,system,thread,program_options --prefix=../boost
./b2 install
cd ..
rm -rf boost_1_60_0
rm -rf boost_1_60_0.tar.gz
cd ..
else
echo "Boost download error!"
fi

#!/bin/bash

apt install -y libboost-all-dev
apt install -y libopencv-dev
apt install -y libeigen3-dev
apt install -y libsndfile1-dev
apt install -y cmake
cd 3rdParty && tar -xvf arma*xz
cp config.hpp armadillo-10.1.2/include/armadillo_bits
cd pugixml && mkdir -p build && cd build && cmake .. && make

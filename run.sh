#!/bin/bash
cd -- "$(dirname "$0")"
mkdir -p build
cd build
cmake -DEXEC_TIME=ON -DCMAKE_BUILD_TYPE=Release ..
make
chmod +x main
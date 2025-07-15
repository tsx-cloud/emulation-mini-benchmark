#!/bin/bash
clang++ --target=x86_64-w64-windows-gnu -std=c++20 -march=znver2 -O2 -fuse-ld=lld \
  -static-libgcc -static-libstdc++ \
  -I/usr/x86_64-w64-mingw32/include \
  -L/usr/x86_64-w64-mingw32/lib \
  -o benchmark_x64.exe benchmark.cpp
  
  
 clang++ --target=x86_64-linux-gnu -std=c++20 -march=znver2 -O2 \
  -I/usr/x86_64-linux-gnu/include \
  -L/usr/x86_64-linux-gnu/lib \
  -o benchmark_x64 benchmark.cpp
  
  
clang++ --target=aarch64-linux-gnu -std=c++20 -march=armv8.2-a -O2 \
  -o benchmark_arm64 benchmark.cpp
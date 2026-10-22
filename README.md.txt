# C++ Boris Pusher

This repository tracks and documents the development of a C++ particle simulation for my master's project.

I’m using it to record progress and checkpoints while learning C++ over the next six months.

## Progress

# 21st October 2026 ( Converting from prior python code ) 
- Implemented a C++ Boris pusher for test particles in a dipole magnetic field.
- Added Earth-radius impact detection and impact-time recording.
- Added OpenMP parallelisation across particles.
- Tested with 1,000,000 particles; 86,000 impacted by `T = 300` in about 120 seconds
- Python --> C++ results in speedup of around 200x 

# 22nd October 2026 ( Converting functions to 3D grid ready for data input ) 
- Implemented a 3D grid from the initial dipole magnetic field function. 
- Implemented linear interpolation for points in between ijk. 
- Compares trajectory timings, errors and mismatches to function version.
- At resolution of 0.0125, relative field error is 6.14159e-05
- At resolution of 0.0125, 8619 particles impact same, 13 mismatch , 3 impacts on analytical and 10 on grid. 
- Hit-time difference: Mean ( 0.0401456 ), Median ( 0.0161152 ), 95th ( 0.172847 ), Max ( 1.13578 ) 
- Increasing duration from 300.0 to 500.0 causes greater mismatch, aswell as 310.764 second difference

## Build

Requires a C++17 compiler with OpenMP support:

```bash
g++ -O3 -fopenmp -std=c++17 -Wall -Wextra -pedantic BorisPusher.cpp -o BorisPusher
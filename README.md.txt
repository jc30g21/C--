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

## Build

Requires a C++17 compiler with OpenMP support:

```bash
g++ -O3 -fopenmp -std=c++17 -Wall -Wextra -pedantic BorisPusher.cpp -o BorisPusher
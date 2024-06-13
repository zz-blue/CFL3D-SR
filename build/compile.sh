#!/bin/bash
# -j on cluster AEROLAB, and -j 8 recommended on cluster TIANHE
# In makefile, add -DFASTIO if cfl3d_mpi, but remove it if cfl3d_seq

make scruball -f makefile_intel -j
make linkall -f makefile_intel -j
make cfl3d_libs -f makefile_intel -j
make cfl3d_mpi  -f makefile_intel -j

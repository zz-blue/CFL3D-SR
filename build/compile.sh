#!/bin/bash
# -j on cluster AEROLAB, and -j 8 recommended on cluster TIANHE
# In makefile, add -DFASTIO if cfl3d_mpi, but remove it if cfl3d_seq

make scruball -f makefile_gnu
make linkall -f makefile_gnu
make cfl3d_libs -f makefile_gnu -j
make cfl3d_mpi  -f makefile_gnu -j

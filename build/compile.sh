#!/bin/bash
# -j on cluster AEROLAB, and -j 8 recommended on cluster TIANHE
# In makefile, add -DFASTIO if cfl3d_mpi, but remove -DFASTIO if cfl3d_seq
make scruball -f makefile_nvidia
make linkall -f makefile_nvidia
make cfl3d_libs -f makefile_nvidia -j
make cfl3d_mpi  -f makefile_nvidia -j

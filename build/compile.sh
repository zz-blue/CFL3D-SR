#!/bin/bash
# use -j on cluster AEROLAB, and use -j 8 on cluster TIANHE
# define -DFASTIO for cfl3d_mpi, and not for cfl3d_seq in makefile
make scruball -f makefile_nvidia
make linkall -f makefile_nvidia
make cfl3d_libs -f makefile_nvidia -j
make cfl3d_mpi  -f makefile_nvidia -j

#!/bin/bash
# cfl3d_mpi for parallel version, cfl3d_seq for sequential version
# makefile_gfortran for GNU compiler, and makefile_ifort for INTEL compiler
# use -j on cluster AEROLAB, and use -j 8 on cluster TIANHE
# define -DFASTIO in makefile for cfl3d_mpi, 
# and do not define if for cfl3d_seq
make scruball -f makefile_gfortran
make linkall -f makefile_gfortran
make cfl3d_libs -f makefile_gfortran -j
make cfl3d_mpi  -f makefile_gfortran -j

c  ---------------------------------------------------------------------------
c  CFL3D is a structured-grid, cell-centered, upwind-biased, Reynolds-averaged
c  Navier-Stokes (RANS) code. It can be run in parallel on multiple grid zones
c  with point-matched, patched, overset, or embedded connectivities. Both
c  multigrid and mesh sequencing are available in time-accurate or
c  steady-state modes.
c
c  Copyright 2001 United States Government as represented by the Administrator
c  of the National Aeronautics and Space Administration. All Rights Reserved.
c 
c  The CFL3D platform is licensed under the Apache License, Version 2.0 
c  (the "License"); you may not use this file except in compliance with the 
c  License. You may obtain a copy of the License at 
c  http://www.apache.org/licenses/LICENSE-2.0. 
c 
c  Unless required by applicable law or agreed to in writing, software 
c  distributed under the License is distributed on an "AS IS" BASIS, WITHOUT 
c  WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the 
c  License for the specific language governing permissions and limitations 
c  under the License.
c  ---------------------------------------------------------------------------
c
      subroutine varnam
c
      open(unit=121,file='field.nam',form='formatted')
      write(121,'(a)') 'RHO' 
      write(121,'(a)') 'U_X'
      write(121,'(a)') 'U_Y'
      write(121,'(a)') 'U_Z'
      write(121,'(a)') 'P'
      write(121,'(a)') 'T'
      write(121,'(a)') 'MUL'
      write(121,'(a)') 'MUT'
      write(121,'(a)') 'OMEGA'
      write(121,'(a)') 'K'
      write(121,'(a)') 'M'
      write(121,'(a)') 'VORT_X'
      write(121,'(a)') 'VORT_Y'
      write(121,'(a)') 'VORT_Z'
      write(121,'(a)') 'Q'
      write(121,'(a)') 'GRAD_RHO'
      write(121,'(a)') 'GRAD_P'
      write(121,'(a)') 'FD'
      write(121,'(a)') 'LRANS'
      write(121,'(a)') 'LLES'
      write(121,'(a)') 'XBLEND'
      close(121)
c
      open(unit=121,file='field_avg.nam',form='formatted')
      write(121,'(a)') 'RHO_AVG' 
      write(121,'(a)') 'RHO_RMS' 
      write(121,'(a)') 'U_AVG_X'
      write(121,'(a)') 'U_AVG_Y'
      write(121,'(a)') 'U_AVG_Z'
      write(121,'(a)') 'U_RMS_X'
      write(121,'(a)') 'U_RMS_Y'
      write(121,'(a)') 'U_RMS_Z'
      write(121,'(a)') 'U_REY_X'
      write(121,'(a)') 'U_REY_Y'
      write(121,'(a)') 'U_REY_Z'
      write(121,'(a)') 'P_AVG'
      write(121,'(a)') 'P_RMS'
      write(121,'(a)') 'T_AVG'
      write(121,'(a)') 'T_RMS'
      write(121,'(a)') 'MUL_AVG'
      write(121,'(a)') 'MUL_RMS'
      write(121,'(a)') 'MUT_AVG'
      write(121,'(a)') 'MUT_RMS'
      write(121,'(a)') 'OMEGA_AVG'
      write(121,'(a)') 'OMEGA_RMS'
      write(121,'(a)') 'K_AVG'
      write(121,'(a)') 'K_RMS'
      write(121,'(a)') 'M_AVG'
      write(121,'(a)') 'M_RMS'
      write(121,'(a)') 'VORT_AVG_X'
      write(121,'(a)') 'VORT_AVG_Y'
      write(121,'(a)') 'VORT_AVG_Z'
      write(121,'(a)') 'VORT_RMS_X'
      write(121,'(a)') 'VORT_RMS_Y'
      write(121,'(a)') 'VORT_RMS_Z'
      write(121,'(a)') 'Q_AVG'
      write(121,'(a)') 'Q_RMS'
      write(121,'(a)') 'GRAD_RHO_AVG'
      write(121,'(a)') 'GRAD_RHO_RMS'
      write(121,'(a)') 'GRAD_P_AVG'
      write(121,'(a)') 'GRAD_P_RMS'
      write(121,'(a)') 'FD_AVG'
      write(121,'(a)') 'FD_RMS'
      write(121,'(a)') 'LRANS_AVG'
      write(121,'(a)') 'LRANS_RMS'
      write(121,'(a)') 'LLES_AVG'
      write(121,'(a)') 'LLES_RMS'
      write(121,'(a)') 'XBLEND_AVG'
      write(121,'(a)') 'XBLEND_RMS'
      close(121)
c
      open(unit=121,file='surf.nam',form='formatted')
      write(121,'(a)') 'U_X' 
      write(121,'(a)') 'U_Y' 
      write(121,'(a)') 'U_Z'
      write(121,'(a)') 'P'
      write(121,'(a)') 'T'
      write(121,'(a)') 'M'
      write(121,'(a)') 'CP'
      write(121,'(a)') 'MUT'
      write(121,'(a)') 'DIS'
      write(121,'(a)') 'CH'
      write(121,'(a)') 'YPLUS'
      write(121,'(a)') 'CF_X'
      write(121,'(a)') 'CF_Y'
      write(121,'(a)') 'CF_Z'
      close(121)
c
      open(unit=121,file='surf_avg.nam',form='formatted')
      write(121,'(a)') 'U_AVG_X'
      write(121,'(a)') 'U_AVG_Y'
      write(121,'(a)') 'U_AVG_Z'
      write(121,'(a)') 'U_RMS_X'
      write(121,'(a)') 'U_RMS_Y'
      write(121,'(a)') 'U_RMS_Z'
      write(121,'(a)') 'P_AVG'
      write(121,'(a)') 'P_RMS'
      write(121,'(a)') 'T_AVG'
      write(121,'(a)') 'T_RMS'
      write(121,'(a)') 'M_AVG'
      write(121,'(a)') 'M_RMS'
      write(121,'(a)') 'CP_AVG'
      write(121,'(a)') 'CP_RMS'
      write(121,'(a)') 'MUT_AVG'
      write(121,'(a)') 'MUT_RMS'
      write(121,'(a)') 'DIS_AVG'
      write(121,'(a)') 'DIS_RMS'
      write(121,'(a)') 'CH_AVG'
      write(121,'(a)') 'CH_RMS'
      write(121,'(a)') 'YPLUS_AVG'
      write(121,'(a)') 'YPLUS_RMS'
      write(121,'(a)') 'CF_AVG_X'
      write(121,'(a)') 'CF_AVG_Y'
      write(121,'(a)') 'CF_AVG_Z'
      write(121,'(a)') 'CF_RMS_X'
      write(121,'(a)') 'CF_RMS_Y'
      write(121,'(a)') 'CF_RMS_Z'
      close(121)
c
      return
      end
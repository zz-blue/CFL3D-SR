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
c  --------------------------------------------------------------------------
c
      subroutine fcd(ax,ay,az,area,at,f,q1,q2,q3,q4,n,nvtq)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Compute central flux at the interface provided four stencil cells
c     Extension to Weiss-Smith's preconditioning: J.R. Edwards, July,1998
c       cprec = 0 ---> original code used
c             > 0 ---> modified code used
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
c       ax ay az unit direction vector, area, at velocity of grid motion
      dimension ax(n),ay(n),az(n),area(n),at(n),
     .          q1(nvtq,5),q2(nvtq,5),q3(nvtq,5),q4(nvtq,5),f(nvtq,5)
c
      common /fluid/ gamma,gm1
      
      x1   = gamma/gm1
c
      !if (real(cprec) .eq. 0.) then   !no preconditioning
cdir$ ivdep
         do i=1,n
c
c         pressure and enthalpy
c
         t6 = 1./12.*(-q1(i,5)+7.*q2(i,5)+7.*q3(i,5)-q4(i,5)) !pressure
         
         q1(i,5) = x1*q1(i,5)/q1(i,1)+0.5*(q1(i,2)*q1(i,2)
     .                   +q1(i,3)*q1(i,3)+q1(i,4)*q1(i,4))
         q2(i,5) = x1*q2(i,5)/q2(i,1)+0.5*(q2(i,2)*q2(i,2)
     .                   +q2(i,3)*q2(i,3)+q2(i,4)*q2(i,4)) 
         q3(i,5) = x1*q3(i,5)/q3(i,1)+0.5*(q3(i,2)*q3(i,2)
     .                   +q3(i,3)*q3(i,3)+q3(i,4)*q3(i,4)) 
         q4(i,5) = x1*q4(i,5)/q4(i,1)+0.5*(q4(i,2)*q4(i,2)
     .                   +q4(i,3)*q4(i,3)+q4(i,4)*q4(i,4))
         
         t1 = 1./12.*(-q1(i,1)+7.*q2(i,1)+7.*q3(i,1)-q4(i,1))
         t2 = 1./12.*(-q1(i,2)+7.*q2(i,2)+7.*q3(i,2)-q4(i,2))
         t3 = 1./12.*(-q1(i,3)+7.*q2(i,3)+7.*q3(i,3)-q4(i,3))
         t4 = 1./12.*(-q1(i,4)+7.*q2(i,4)+7.*q3(i,4)-q4(i,4))
         t5 = 1./12.*(-q1(i,5)+7.*q2(i,5)+7.*q3(i,5)-q4(i,5)) !enthalpy
         
c
         t0  = ax(i)*t2+ay(i)*t3+az(i)*t4+at(i)
         f1  = t0*t1
         f2  = f1*t2        
         f3  = f1*t3
         f4  = f1*t4
         f5  = f1*t5
       
         f2  = f2+ax(i)*t6
         f3  = f3+ay(i)*t6
         f4  = f4+az(i)*t6
         f5  = f5-at(i)*t6
c
c         include factor area
c
         f(i,1) = f1*area(i)
         f(i,2) = f2*area(i)
         f(i,3) = f3*area(i)
         f(i,4) = f4*area(i)
         f(i,5) = f5*area(i)
         
         end do     
      return
      end
      
      
      subroutine flxblndcof(jdim,kdim,idim,xlscale,xrans,xblend)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Compute flux blending coefficient
c
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension xlscale(jdim-1,kdim-1,idim-1),
     +          xrans(jdim-1,kdim-1,idim-1)
      dimension xblend(jdim,kdim,idim)
c
      common /twod/ i2d
      common /axisym/ iaxi2plane,iaxi2planeturb,istrongturbdis,iforcev0
      common /numerics/ iblend,lxblend,blendlim(2)
c
      if(iblend.eq.1) then
        do i=1,idim-1
          do j=1,jdim-1
            do k=1,kdim-1      
              xlhyb = xlscale(j,k,i)
              xrs   = xrans(j,k,i)
              xtmp  = cctanh(10.*ccmax(xlhyb/xrs-0.6,0.0))
              xtmp  = ccmax(xtmp,blendlim(1))
              xblend(j,k,i) = ccmin(xtmp,blendlim(2))
            end do
          end do
        end do
      else 
c     please add other flux blending functions
      end if
c     flux conservation, to be paralleled dim = 1 analogous to vi0,vj0,vk0 for vist3d
      xblend(1,1:kdim,1:idim)    = blendlim(2)
      xblend(1:jdim,1,1:idim)    = blendlim(2)
      xblend(jdim,1:kdim,1:idim) = blendlim(2)
      xblend(1:jdim,kdim,1:idim) = blendlim(2)
      if( i2d .ne. 1 .and. iaxi2planeturb .ne. 1 ) then
        xblend(1:jdim,1:kdim,1)    = blendlim(2)
        xblend(1:jdim,1:kdim,idim) = blendlim(2)
      end if
c
      return
      end

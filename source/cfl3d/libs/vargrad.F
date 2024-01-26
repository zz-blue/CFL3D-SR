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
c  -------------------------------------------------------------------------
c 
      subroutine vargrad(jdim,kdim,idim,q,sj,sk,si,vol,ux,wt,blank,
     .                   iover,qj0,qk0,qi0,bcj,bck,bci,nbl,volj0,volk0,
     .                   voli0,maxbl,vormax,ivmax,jvmax,kvmax,px,rx,
     .                   vdsp,nvdsp)
c
c     $Id$
c
c**********************************************************************
c     Purpose:  Evaluate velocity derivatives at cell centers
c               and compute gradient-based quantities vorticity, 
c               Q, divrho, divp
c**********************************************************************
c
c      input arrays   : 
c       q             : primitive variables at cell centers
c                       (rho,u,v,w,p)
c       vol           : cell volumes
c       sj,sk,si      : metrics
c                       (direction cosines,areas,speeds of cell faces)
c       ifj,ifk,ifi   : =0 don't compute contribution from that direction
c                       >0 do      "
c      output arrays  :
c       ux            : velocity derivatives    
c                           1 = ux    4 = vx    7 = wx
c                           2 = uy    5 = vy    8 = wy
c                           3 = uz    6 = vz    9 = wz
c      scratch arrays :   
c       wt(1-3)       : cross-flow temporaries
c
c**********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension q(jdim,kdim,idim,5),ux(jdim-1,kdim-1,idim-1,9)
      dimension sj(jdim,kdim,idim-1,5),sk(jdim,kdim,idim-1,5),
     .          si(jdim,kdim,idim,5)
      dimension vol(jdim,kdim,idim-1),blank(jdim,kdim,idim)
      dimension wt(jdim,kdim,9)
      dimension bcj(kdim,idim-1,2),bck(jdim,idim-1,2),bci(jdim,kdim,2)
      dimension qj0(kdim,idim-1,5,4),
     .          qk0(jdim,idim-1,5,4),qi0(jdim,kdim,5,4)
      dimension volj0(kdim,idim-1,4),
     .          volk0(jdim,idim-1,4),voli0(jdim,kdim,4)
      dimension vormax(maxbl),ivmax(maxbl),jvmax(maxbl),kvmax(maxbl)
      dimension px(jdim-1,kdim-1,idim-1,3),rx(jdim-1,kdim-1,idim-1,3),
     .          vdsp(jdim,kdim,idim,nvdsp)
c
      common /twod/ i2d
c
c
      call delv(jdim,kdim,idim,q,sj,sk,si,vol,ux,wt,
     .          blank,iover,qj0,qk0,qi0,bcj,bck,bci,nblk,
     .          volj0,volk0,voli0,maxbl,vormax,ivmax,jvmax,kvmax)
      ldim = 1
      call delq(jdim,kdim,idim,q,sj,sk,si,vol,rx,wt,
     .          blank,iover,qj0,qk0,qi0,bcj,bck,bci,nblk,
     .          volj0,volk0,voli0,maxbl,vormax,ivmax,jvmax,kvmax,
     .          1)
      ldim = 5
      call delq(jdim,kdim,idim,q,sj,sk,si,vol,px,wt,
     .          blank,iover,qj0,qk0,qi0,bcj,bck,bci,nblk,
     .          volj0,volk0,voli0,maxbl,vormax,ivmax,jvmax,kvmax,
     .          ldim)
c
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
c
      do i=1,idim1
      do k=1,kdim1
      do j=1,jdim1
c
        s11 = ux(j,k,i,1)
        s22 = ux(j,k,i,5)
        s33 = ux(j,k,i,9)
        s12 = 0.5*(ux(j,k,i,2) + ux(j,k,i,4))
        s13 = 0.5*(ux(j,k,i,3) + ux(j,k,i,7))
        s23 = 0.5*(ux(j,k,i,6) + ux(j,k,i,8))
        xis = s11*s11 + s22*s22 + s33*s33 +
     .        2.*s12*s12 + 2.*s13*s13 + 2.*s23*s23
        w12 = 0.5*(ux(j,k,i,2) - ux(j,k,i,4))
        w13 = 0.5*(ux(j,k,i,3) - ux(j,k,i,7))
        w23 = 0.5*(ux(j,k,i,6) - ux(j,k,i,8))
        wis = 2.*w12*w12 + 2.*w13*w13 + 2.*w23*w23
c
        gradr = rx(j,k,i,1)*rx(j,k,i,1)
     .        + rx(j,k,i,2)*rx(j,k,i,2)
     .        + rx(j,k,i,3)*rx(j,k,i,3)
        gradp = px(j,k,i,1)*px(j,k,i,1)
     .        + px(j,k,i,2)*px(j,k,i,2)
     .        + px(j,k,i,3)*px(j,k,i,3)
c
c       vortx,vorty,vortz
        vdsp(j,k,i,1) = -2.*w23
        vdsp(j,k,i,2) =  2.*w13
        vdsp(j,k,i,3) = -2.*w12
c       q0
        vdsp(j,k,i,4) = 0.5*(wis-xis)
c       magnitude of gradr,gradp
        vdsp(j,k,i,5) = sqrt(gradr)
        vdsp(j,k,i,6) = sqrt(gradp)
c
      end do
      end do
      end do
c
      return
      end

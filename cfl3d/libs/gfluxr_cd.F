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
      subroutine gfluxr_cd(i,npl,jdim,kdim,idim,q,qj0,sj,t,nvtq,bcj)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Compute central flux for the 
c     right-hand-side in the J-direction from the inviscid terms.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
      dimension sj(jdim*kdim,idim-1,5)
      dimension q(jdim,kdim,idim,5),qj0(kdim,idim-1,5,4)
      dimension bcj(kdim,idim-1,2)
      dimension t(nvtq,30)
      !1 ~20: 4 stencil cells
      !21~25: central flux
c
      idim1 = idim-1
      jdim1 = jdim-1
      kdim1 = kdim-1
c         
      n     = jdim*kdim*npl
      nr    = n-1
c
c      higher order
c
      do l=1,5
c
c     interior face 3~(nr=n-1)
c
cdir$ ivdep
      do izz=3,nr
      t(izz,0 +l) = q(izz-2,1,i,l)
      t(izz,5 +l) = q(izz-1,1,i,l)
      t(izz,10+l) = q(izz  ,1,i,l)
      t(izz,15+l) = q(izz+1,1,i,l)
      end do
c
c     left boundary face 1
c
      do ipl=1,npl
      ii = i+ipl-1
      jc = (ipl-1)*(kdim*jdim) + 1 - jdim
      do kk=1,kdim1
      jc = jc + jdim
      !cell-face type, t(1~4) = qj0(1)
      tc = 1.0-bcj(kk,ii,1)
      tf =     bcj(kk,ii,1)
      t(jc,0 +l) = tc*qj0(kk,ii,l,2)
     .             + tf*qj0(kk,ii,l,1)
      t(jc,5 +l) = tc*qj0(kk,ii,l,1)
     .             + tf*qj0(kk,ii,l,1)
      t(jc,10+l) = tc*q(1,kk,ii,l)
     .             + tf*qj0(kk,ii,l,1)
      t(jc,15+l) = tc*q(2,kk,ii,l)  
     .             + tf*qj0(kk,ii,l,1)
      end do
      end do
c
c     left boundary face 2
c
      if(jdim1.ne.2) then
      do ipl=1,npl
      ii = i+ipl-1
      jc = (ipl-1)*(kdim*jdim) + 2 - jdim
      do kk=1,kdim1
      jc = jc + jdim
      !cell-face type, t(1) = 2.*qj0(1)-q(1)
      tc = 1.0-bcj(kk,ii,1)
      tf =     bcj(kk,ii,1)
      t(jc,0 +l) = tc*qj0(kk,ii,l,1)
     .             + tf*(2.*qj0(kk,ii,l,1)-q(1,kk,ii,l))
      t(jc,5 +l) = tc*q(1,kk,ii,l)
     .             + tf*q(1,kk,ii,l)
      t(jc,10+l) = tc*q(2,kk,ii,l)
     .             + tf*q(2,kk,ii,l)
      t(jc,15+l) = tc*q(3,kk,ii,l)
     .             + tf*q(3,kk,ii,l)
      end do
      end do
      end if
c
c     right boundary face jdim
c
      do ipl=1,npl
      ii = i+ipl-1
      jc = (ipl-1)*(kdim*jdim)
      do kk=1,kdim1
      jc = jc + jdim
      !cell-face type, t(1~4) = qj0(3)
      tc = 1.0-bcj(kk,ii,2)
      tf =     bcj(kk,ii,2)
      t(jc,0 +l) = tc*q(jdim1-1,kk,ii,l)
     .             + tf*qj0(kk,ii,l,3)
      t(jc,5 +l) = tc*q(jdim1,kk,ii,l)
     .             + tf*qj0(kk,ii,l,3)
      t(jc,10+l) = tc*qj0(kk,ii,l,3)
     .             + tf*qj0(kk,ii,l,3)
      t(jc,15+l) = tc*qj0(kk,ii,l,4)
     .             + tf*qj0(kk,ii,l,3)
      end do
      end do
c
c     right boundary face jdim1=jdim-1
c     
      if(jdim1.ne.2) then
      do ipl=1,npl
      ii = i+ipl-1
      jc = (ipl-1)*(kdim*jdim) - 1
      do kk=1,kdim1
      jc = jc + jdim
      !cell-face type, t(4) = 2.*qj0(3)-q(jdim1)
      tc = 1.0-bcj(kk,ii,2)
      tf =     bcj(kk,ii,2)
      t(jc,0 +l) = tc*q(jdim1-2,kk,ii,l)
     .             + tf*q(jdim1-2,kk,ii,l)
      t(jc,5 +l) = tc*q(jdim1-1,kk,ii,l)
     .             + tf*q(jdim1-1,kk,ii,l)
      t(jc,10+l) = tc*q(jdim1  ,kk,ii,l)
     .             + tf*q(jdim1  ,kk,ii,l)
      t(jc,15+l) = tc*qj0(kk,ii,l,3)
     .             + tf*(2.*qj0(kk,ii,l,3)-q(jdim1,kk,ii,l))
      end do
      end do
      end if
c
c     special treatment for face(2)/face(jdim1) if jdim1=2
c          
      if(jdim1.eq.2) then
      do ipl=1,npl
      ii = i+ipl-1
      jc = (ipl-1)*(kdim*jdim) + 2 - jdim
      do kk=1,kdim1
      jc = jc + jdim
      !q(0)/q(3) replaced by ghost cells qj0(1)/qj0(3)
      !q(1)/q(2) from above
      tc = 1.0-bcj(kk,ii,1)
      tf =     bcj(kk,ii,1)
      t(jc,0 +l) = tc*qj0(kk,ii,l,1)
     .             + tf*(2.*qj0(kk,ii,l,1)-q(1,kk,ii,l))
      t(jc,5 +l) = tc*q(1,kk,ii,l)
     .             + tf*q(1,kk,ii,l)
      tc = 1.0-bcj(kk,ii,2)
      tf =     bcj(kk,ii,2)
      t(jc,10+l) = tc*q(2,kk,ii,l)
     .             + tf*q(2,kk,ii,l)
      t(jc,15+l) = tc*qj0(kk,ii,l,3)
     .             + tf*(2.*qj0(kk,ii,l,3)-q(2,kk,ii,l))
      end do
      end do
      endif
c
c     fill end point for safety (just for completeness)
      t(n,10+l) = t(nr,10+l)
      t(n,15+l) = t(nr,15+l)
c
      end do
c
      jkpro = jdim*kdim*npl-jdim
c
      call fcd(sj(1,i,1),sj(1,i,2),sj(1,i,3),sj(1,i,4),sj(1,i,5),
     .         t(1,21),t(1,1),t(1,6),t(1,11),t(1,16),
     .         jkpro,nvtq)
      
      return
      end

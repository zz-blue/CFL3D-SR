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
      subroutine hfluxr_cd(i,npl,jdim,kdim,idim,q,qk0,sk,t,nvtq,bck)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Compute central flux for the 
c     right-hand-side in the K-direction from the inviscid terms
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
      dimension sk(jdim*kdim,idim-1,5)
      dimension q(jdim,kdim,idim,5),qk0(jdim,idim-1,5,4)
      dimension bck(jdim,idim-1,2)
      dimension t(nvtq,30)
      !1 ~20: 4 stencil cells
      !21~25: central flux
c
      idim1 = idim-1
      jdim1 = jdim-1
      kdim1 = kdim-1
c         
      n     = jdim*kdim*npl
      nr    = n-jdim
      nr2   = n-2*jdim
c
c      higher order
c
      do l=1,5
c
c     interior face (1+2*jdim)~(nr2=n-2*jdim)
c
cdir$ ivdep
      do izz=1+2*jdim,nr2
      t(izz,0 +l) = q(izz-2*jdim,1,i,l)
      t(izz,5 +l) = q(izz-1*jdim,1,i,l)
      t(izz,10+l) = q(izz       ,1,i,l)
      t(izz,15+l) = q(izz+1*jdim,1,i,l)
      end do
c
c     left boundary face 1
c
      do ipl=1,npl
      ii = i+ipl-1
      jc = (ipl-1)*(kdim*jdim)
      do jj=1,jdim1
      jc = jc + 1
      !cell-face type, t(1~4) = qk0(1)
      tc = 1.0-bck(jj,ii,1)
      tf =     bck(jj,ii,1)
      t(jc,0 +l) = tc*qk0(jj,ii,l,2)
     .             + tf*qk0(jj,ii,l,1)
      t(jc,5 +l) = tc*qk0(jj,ii,l,1)
     .             + tf*qk0(jj,ii,l,1)
      t(jc,10+l) = tc*q(jj,1,ii,l)  
     .             + tf*qk0(jj,ii,l,1)
      t(jc,15+l) = tc*q(jj,2,ii,l)
     .             + tf*qk0(jj,ii,l,1)
      end do
      end do
c
c     left boundary face 2 if kdim1.ne.2
c
      if(kdim1.ne.2) then
      do ipl=1,npl
      ii = i+ipl-1
      jc = (ipl-1)*(kdim*jdim) + jdim
      do jj=1,jdim1
      jc = jc + 1
      !cell-face type, t(1) = 2.*qk0(1)-q(1)
      tc = 1.0-bck(jj,ii,1)
      tf =     bck(jj,ii,1)
      t(jc,0 +l) = tc*qk0(jj,ii,l,1)
     .             + tf*(2.*qk0(jj,ii,l,1)-q(jj,1,ii,l))
      t(jc,5 +l) = tc*q(jj,1,ii,l)
     .             + tf*q(jj,1,ii,l)
      t(jc,10+l) = tc*q(jj,2,ii,l)
     .             + tf*q(jj,2,ii,l)
      t(jc,15+l) = tc*q(jj,3,ii,l)
     .             + tf*q(jj,3,ii,l)
      end do
      end do   
      end if
c
c     right boundary face kdim
c
      do ipl=1,npl
      ii = i+ipl-1
      jc = (ipl-1)*(kdim*jdim) + kdim1*jdim
      do jj=1,jdim1
      jc = jc + 1
      !cell-face type, t(1~4) = qk0(3)
      tc = 1.0-bck(jj,ii,2)
      tf =     bck(jj,ii,2)
      t(jc,0 +l) = tc*q(jj,kdim1-1,ii,l)
     .             + tf*qk0(jj,ii,l,3)
      t(jc,5 +l) = tc*q(jj,kdim1,ii,l)
     .             + tf*qk0(jj,ii,l,3)
      t(jc,10+l) = tc*qk0(jj,ii,l,3)
     .             + tf*qk0(jj,ii,l,3)
      t(jc,15+l) = tc*qk0(jj,ii,l,4)
     .             + tf*qk0(jj,ii,l,3)
      end do
      end do      
c
c     right boundary face kdim1=kdim-1 if kdim1.ne.2
c     
      if(kdim1.ne.2) then
      do ipl=1,npl
      ii = i+ipl-1
      jc = (ipl-1)*(kdim*jdim) + kdim1*jdim - jdim
      do jj=1,jdim1
      jc = jc + 1
      !cell-face type, t(4) = 2.*qk0(3)-q(kdim1)
      tc = 1.0-bck(jj,ii,2)
      tf =     bck(jj,ii,2)
      t(jc,0 +l) = tc*q(jj,kdim1-2,ii,l)        
     .             + tf*q(jj,kdim1-2,ii,l)
      t(jc,5 +l) = tc*q(jj,kdim1-1,ii,l)
     .             + tf*q(jj,kdim1-1,ii,l)
      t(jc,10+l) = tc*q(jj,kdim1  ,ii,l)
     .             + tf*q(jj,kdim1  ,ii,l)
      t(jc,15+l) = tc*qk0(jj,ii,l,3)
     .             + tf*(2.*qk0(jj,ii,l,3)-q(jj,kdim1,ii,l))
      end do
      end do
      end if
c
c     special treatment for face(2)/face(kdim1) if kdim1=2
c          
      if(kdim1.eq.2) then
      do ipl=1,npl
      ii = i+ipl-1
      jc = (ipl-1)*(kdim*jdim) + jdim
      do jj=1,jdim1
      jc = jc + 1
      !q(0)/q(3) replaced by ghost cells qj0(1)/qj0(3)
      !q(1)/q(2) from above
      tc = 1.0-bck(jj,ii,1)
      tf =     bck(jj,ii,1)
      t(jc,0 +l) = tc*qk0(jj,ii,l,1)
     .             + tf*(2.*qk0(jj,ii,l,1)-q(jj,1,ii,l))
      t(jc,5 +l) = tc*q(jj,1,ii,l)
     .             + tf*q(jj,1,ii,l)
      tc = 1.0-bck(jj,ii,2)
      tf =     bck(jj,ii,2)
      t(jc,10+l) = tc*q(jj,2,ii,l)
     .             + tf*q(jj,2,ii,l)
      t(jc,15+l) = tc*qk0(jj,ii,l,3)
     .             + tf*(2.*qk0(jj,ii,l,3)-q(jj,2,ii,l))
      end do
      end do
      endif
c
c     fill end points for safety (just for completeness)
      t(jdim,0 +l) = t(jdim-1,0 +l)
      t(jdim,5 +l) = t(jdim-1,5 +l)
      t(jdim*kdim*npl,10+l) = t(jdim*kdim*npl-1,10+l)
      t(jdim*kdim*npl,15+l) = t(jdim*kdim*npl-1,15+l)
c
      end do
c
      call fcd(sk(1,i,1),sk(1,i,2),sk(1,i,3),sk(1,i,4),sk(1,i,5),
     .         t(1,21),t(1,1),t(1,6),t(1,11),t(1,16),
     .         n,nvtq)
            
      return
      end

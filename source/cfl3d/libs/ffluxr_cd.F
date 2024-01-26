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
      subroutine ffluxr_cd(k,npl,jdim,kdim,idim,q,qi0,si,t,nvtq,bci)
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
      dimension si(jdim*kdim,idim,5)
      dimension q(jdim,kdim,idim,5),qi0(jdim,kdim,5,4)
      dimension bci(jdim,kdim,2)
      !dimension t(nvtq,25),tsi(nvtq,5)
      dimension t(nvtq,30)
      !t(26:30): si
c
      idim1 = idim-1
      jdim1 = jdim-1
      kdim1 = kdim-1
c         
      if (npl.eq.kdim1 .and. nvtq.ge.jdim*kdim*idim) npl = kdim
c
      jv = npl*jdim
      n  = jv*idim
      nr = n-jv
      nr2= n-2*jv
c
      if (npl.ne.kdim) then   
c
c        fill temp arrays with metrics if all k planes are not done at once
c
         do i=1,idim
         jk  = (i-1)*jv+1
         jk1 = (k-1)*jdim+1
         do l=1,5
cdir$ ivdep
         do izz=1,jv
         t(izz+jk-1,25+l) = si(izz+jk1-1,i,l)
         end do
         end do
         end do
      end if
c
c     all orders
c
c     fill temp arrays with interior values
c
      do l=1,5
      if (npl.eq.kdim) then
cdir$ ivdep
         do 50 izz=1,nr
         t(izz+jv,20+l) = q(izz,1,1,l)
   50    continue
      else                 
         do 70 i=1,idim1
         jk = (i-1)*jv+1
cdir$ ivdep
         do 60 izz=1,jv
         t(izz+jv+jk-1,20+l) = q(izz,k,i,l)
   60    continue
   70    continue
      end if
c
c      fill temp array with boundary values
c
cdir$ ivdep
      do izz=1,jv
      t(izz,20+l)    = qi0(izz,k,l,1)
      end do
      end do
      
c      higher order
c
      do l=1,5
c
c     interior face (1+2*jv)~(nr2=n-2*jv)
c
cdir$ ivdep
      do izz=1+2*jv,nr2
      t(izz,0 +l) = t(izz-jv  ,20+l)
      t(izz,5 +l) = t(izz     ,20+l)
      t(izz,10+l) = t(izz+jv  ,20+l)
      t(izz,15+l) = t(izz+2*jv,20+l)
      end do
c
c     left boundary face 1
c
      do kpl=1,npl
      kk = k+kpl-1
      jc = (kpl-1)*jdim
      do jj=1,jdim1
      jc = jc + 1
      !cell-face type, t(1~4) = qi0(1)
      tc = 1.0-bci(jj,kk,1)
      tf =     bci(jj,kk,1)
      t(jc,0 +l) = tc*qi0(jj,kk,l,2)
     .             + tf*qi0(jj,kk,l,1)
      t(jc,5 +l) = tc*qi0(jj,kk,l,1)
     .             + tf*qi0(jj,kk,l,1)
      t(jc,10+l) = tc*q(jj,kk,1,l)
     .             + tf*qi0(jj,kk,l,1) 
      t(jc,15+l) = tc*q(jj,kk,2,l)  
     .             + tf*qi0(jj,kk,l,1)
      end do
      end do
c
c     left boundary face 2 if idim1.ne.2
c
      if(idim2.ne.2) then
      do kpl=1,npl
      kk = k+kpl-1
      jc = (kpl-1)*jdim + jv
      do jj=1,jdim1
      jc = jc + 1
      !cell-face type, t(1) = 2.*qi0(1)-q(1)
      tc = 1.0-bci(jj,kk,1)
      tf =     bci(jj,kk,1)
      t(jc,0 +l) = tc*qi0(jj,kk,l,1)
     .             + tf*(2.*qi0(jj,kk,l,1)-q(jj,kk,1,l))
      t(jc,5 +l) = tc*q(jj,kk,1,l)
     .             + tf*q(jj,kk,1,l)
      t(jc,10+l) = tc*q(jj,kk,2,l)
     .             + tf*q(jj,kk,2,l)
      t(jc,15+l) = tc*q(jj,kk,3,l)
     .             + tf*q(jj,kk,3,l)
      end do
      end do  
      end if
c
c     right boundary face idim
c
      do kpl=1,npl
      kk = k+kpl-1
      jc = idim1*jv + (kpl-1)*jdim
      do jj=1,jdim1
      jc = jc + 1
      !cell-face type, t(1~4) = qi0(3)
      tc = 1.0-bci(jj,kk,2)
      tf =     bci(jj,kk,2)
      t(jc,0 +l) = tc*q(jj,kk,idim1-1,l)
     .             + tf*qi0(jj,kk,l,3)
      t(jc,5 +l) = tc*q(jj,kk,idim1  ,l)
     .             + tf*qi0(jj,kk,l,3)
      t(jc,10+l) = tc*qi0(jj,kk,l,3)
     .             + tf*qi0(jj,kk,l,3)
      t(jc,15+l) = tc*qi0(jj,kk,l,4)
     .             + tf*qi0(jj,kk,l,3)
      end do
      end do            
c
c     right boundary face idim1=idim-1
c     
      if(idim2.ne.2) then
      do kpl=1,npl
      kk = k+kpl-1
      jc = idim1*jv + (kpl-1)*jdim - jv
      do jj=1,jdim1
      jc = jc + 1
      !cell-face type, t(4) = 2.*qi0(3)-q(idim1)
      tc = 1.0-bci(jj,kk,2)
      tf =     bci(jj,kk,2)
      t(jc,0 +l) = tc*q(jj,kk,idim1-2,l)
     .             + tf*q(jj,kk,idim1-2,l)
      t(jc,5 +l) = tc*q(jj,kk,idim1-1,l)
     .             + tf*q(jj,kk,idim1-1,l)
      t(jc,10+l) = tc*q(jj,kk,idim1  ,l)
     .             + tf*q(jj,kk,idim1  ,l)
      t(jc,15+l) = tc*qi0(jj,kk,l,3)
     .             + tf*(2.*qi0(jj,kk,l,3)-q(jj,kk,idim1,l))
      end do
      end do  
      end if
c
c     special treatment for face(2)/face(idim1) if idim1=2
c          
      if(idim1.eq.2) then
      do kpl=1,npl
      kk = k+kpl-1
      jc = (kpl-1)*jdim + jv
      do jj=1,jdim1
      jc = jc + 1
      !q(0)/q(3) replaced by ghost cells qi0(1)/qi0(3)
      !q(1)/q(2) from above
      tc = 1.0-bci(jj,kk,1)
      tf =     bci(jj,kk,1)
      t(jc,0 +l) = tc*qi0(jj,kk,l,1)
     .             + tf*(2.*qi0(jj,kk,l,1)-q(jj,kk,1,l))
      t(jc,5 +l) = tc*q(jj,kk,1,l)
     .             + tf*q(jj,kk,1,l)
      tc = 1.0-bci(jj,kk,2)
      tf =     bci(jj,kk,2)
      t(jc,10+l) = tc*q(jj,kk,2,l)
     .             + tf*q(jj,kk,2,l)
      t(jc,15+l) = tc*qi0(jj,kk,l,3)
     .             + tf*(2.*qi0(jj,kk,l,3)-q(jj,kk,2,l))
      end do
      end do
      endif      
c
c     fill end points for safety (just for completeness)
      t(jdim,0 +l) = t(jdim-1,0 +l)
      t(jdim,5 +l) = t(jdim-1,5 +l)
      t(jdim*kdim*npl,10+l) = t(jdim*kdim*npl-1,10+l)
      t(jdim*kdim*npl,15+l) = t(jdim*kdim*npl-1,15+l)
      
c     fill end points for safety
      do kpl=1,npl
      jc = (kpl-1)*jdim+jdim
      jl = idim1*jv + (kpl-1)*jdim + jdim
      t(jc,0 +l) = t(jc-1,0 +l)
      t(jc,5 +l) = t(jc-1,5 +l)
      t(jl,10+l) = t(jl-1,10+l)
      t(jl,15+l) = t(jl-1,15+l)
      end do
c
      end do
     
      if (npl.eq.kdim) then
          call fcd(si(1,1,1),si(1,1,2),si(1,1,3),si(1,1,4),si(1,1,5),
     .         t(1,21),t(1,1),t(1,6),t(1,11),t(1,16),
     .         n,nvtq)
      else
          call fcd(t(1,26),t(1,27),t(1,28),t(1,29),t(1,30),
     .         t(1,21),t(1,1),t(1,6),t(1,11),t(1,16),
     .         n,nvtq)
      end if
      
      return
      end

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
      subroutine bc_vdsp(jdim,kdim,idim,dum,dumj0,dumk0,dumi0,
     .                   bcj,bck,bci,blank,ldim)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Install ghost values in vdj0,vdk0,vdi0
c           to obtain vdsp values at grid points.
c           Assume all boundaries are cell-center type. There are errors
c           on wall boundaries, but it is unimportant.
c
c     dum(j=1,jdim;k=1,kdim;i=1,idim;l=1,ldim)
c
c     dumj0(k=1,kdim;i=1,idim-1;l=1,ldim;mm=1,4)
c
c     dumk0(j=1,jdim;i=1,idim-1;l=1,ldim;mm=1,4)
c
c     dumi0(j=1,jdim;k=1,kdim;  l=1,ldim;mm=1,4)
c
c     ldim is the dimension covering the number of solution variables:
c
c     ldim=nvdsp for vdsp/vdj0/vdk0/vdi0
c
c     The bcj/bck/bci arrays indicate whether the data stored in the
c     j0/k0/i0 boundary data arrays in the mm=1 and mm=3 arrays are
c     already face-center based (array entry = 1) or at ghost cell
c     centers (array entry = 1).
c
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension bcj(kdim,idim-1,2),bck(jdim,idim-1,2),bci(jdim,kdim,2)
      dimension dum(jdim,kdim,idim,ldim), dumj0(kdim,idim-1,ldim,4),
     .          dumk0(jdim,idim-1,ldim,4),dumi0(jdim,kdim,ldim,4)
      dimension blank(jdim,kdim,idim)
c
      common /twod/ i2d
c
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
c
      do 9000 l=1,ldim
c
c     j - boundaries
c
      do 6400 m=1,2
c
      imin = 1
      imax = idim-1
      kmin = 1
      kmax = kdim-1
c
      if (m.eq.1) then
c
c        left boundary
c
         do ii=imin,imax
         do kk=kmin,kmax
         dumj0(kk,ii,l,1) = 2.*dum(1,kk,ii,l)
     .                    - 1.*dum(2,kk,ii,l)
         end do
         end do
c
      else
c
c        right boundary
c
         do ii=imin,imax
         do kk=kmin,kmax
         dumj0(kk,ii,l,3) = 2.*dum(jdim-1,kk,ii,l) 
     .                    - 1.*dum(jdim-2,kk,ii,l)
         end do
         end do
c
      end if
 6400 continue
c
c     k - boundaries
c
      do 7400 m=1,2
c
      imin = 1
      imax = idim-1
      jmin = 1
      jmax = jdim-1
c
      if (m.eq.1) then
c
c        left boundary
c
         do ii=imin,imax
         do jj=jmin,jmax
         dumk0(jj,ii,l,1) = 2.*dum(jj,1,ii,l)
     .                    - 1.*dum(jj,2,ii,l)
         end do
         end do
c
      else
c
c        right boundary
c
         do ii=imin,imax
         do jj=jmin,jmax
         dumk0(jj,ii,l,3) = 2.*dum(jj,kdim-1,ii,l)
     .                    - 1.*dum(jj,kdim-2,ii,l)
         end do
         end do
      end if
 7400 continue
c
c     i - boundaries
c
      do 8400 m=1,2
      kmin = 1
      kmax = kdim-1
      jmin = 1
      jmax = jdim-1
c
      if (m.eq.1) then
c
c        left boundary
c
         if(i2d.eq.1) then
           do kk=kmin,kmax
           do jj=jmin,jmax
           dumi0(jj,kk,l,1) = dum(jj,kk,1,l)
           end do
           end do
         else
           do kk=kmin,kmax
           do jj=jmin,jmax
           dumi0(jj,kk,l,1) = 2.*dum(jj,kk,1,l)
     .                      - 1.*dum(jj,kk,2,l)
           end do
           end do
         end if
c
      else
c
c      right boundary
c
         if(i2d.eq.1) then        
           do kk=kmin,kmax
           do jj=jmin,jmax
           dumi0(jj,kk,l,3)  = dum(jj,kk,idim-1,l)
           end do
           end do
         else
           do kk=kmin,kmax
           do jj=jmin,jmax
           dumi0(jj,kk,l,3)  = 2.*dum(jj,kk,idim-1,l)
     .                       - 1.*dum(jj,kk,idim-2,l)
           end do
           end do
         end if
      end if
 8400 continue
c
 9000 continue   
      return
      end 

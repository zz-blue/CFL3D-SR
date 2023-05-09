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
      subroutine xlsfree(jdim,kdim,idim,q,sj,sk,si,vol,ux,xlesfr)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Compute LES length scale away from walls
c
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
      dimension q(jdim,kdim,idim,5),sj(jdim,kdim,idim-1,5),
     + sk(jdim,kdim,idim-1,5),si(jdim,kdim,idim,5),vol(jdim,kdim,idim-1)
      dimension ux(jdim-1,kdim-1,idim-1,9)
      dimension xlesfr(jdim-1,kdim-1,idim-1)
      dimension a(3,3),uxt(3,3),s(3,3),delta(3),w(3),b(3),c(3)
c
      common /twod/ i2d
      common /zero/ iexp
      common /axisym/ iaxi2plane,iaxi2planeturb,istrongturbdis,iforcev0
      common /des/ cdes,ides,cddes,cdesamd,lfe,cdt1,clim,lnwsc
c
      if(ides.eq.11) then
c     maximum length of edges
        do i=1,idim-1
          do j=1,jdim-1
            do k=1,kdim-1
              deltaj = 2.*vol(j,k,i)/(sj(j,k,i,4)+sj(j+1,k,i,4))
              deltak = 2.*vol(j,k,i)/(sk(j,k,i,4)+sk(j,k+1,i,4))
              deltai = 2.*vol(j,k,i)/(si(j,k,i,4)+si(j,k,i+1,4))
              deltam = ccmax(deltaj,deltak)
              if( i2d .ne. 1 .and. iaxi2planeturb .ne. 1 ) then
                  deltam = ccmax(deltam,deltai)  !Modification to allow DES to run in 2d
              end if
              xlesfr(j,k,i) = cdes*deltam
            enddo
          enddo
        enddo
c
      else if(ides.eq.12) then
c     AMD anisotropic length scale, developed by Xiao, Maochao 
c     based on the AMD SGS model of Rozema, W. et al,
c     Minimum-Dissipation Models for Large-Eddy Simulation
c     implemetation of anisotropic models/length scales, refer to
c     Haering, S. W., Resolution-Induced Anisotropy in Large-Eddy Simulation
        do i=1,idim-1
          do j=1,jdim-1
            do k=1,kdim-1
            
              delta(1) = 2.*vol(j,k,i)/(si(j,k,i,4)+si(j,k,i+1,4))
              delta(2) = 2.*vol(j,k,i)/(sj(j,k,i,4)+sj(j+1,k,i,4))
              delta(3) = 2.*vol(j,k,i)/(sk(j,k,i,4)+sk(j,k+1,i,4))
              if( i2d .eq. 1 .or. iaxi2planeturb .eq. 1 ) then
                  delta(1) = 0.0  !Modification to allow 2d
              end if
              
c             direction cosines of face i,j,k
              a(1:3,1) = si(j,k,i,1:3)+si(j,k,i+1,1:3)
              a(1:3,2) = sj(j,k,i,1:3)+sj(j+1,k,i,1:3)
              a(1:3,3) = sk(j,k,i,1:3)+sk(j,k+1,i,1:3)
              xtmp1 = sqrt(a(1,1)*a(1,1)+a(2,1)*a(2,1)+a(3,1)*a(3,1))
              xtmp2 = sqrt(a(1,2)*a(1,2)+a(2,2)*a(2,2)+a(3,2)*a(3,2))
              xtmp3 = sqrt(a(1,3)*a(1,3)+a(2,3)*a(2,3)+a(3,3)*a(3,3))
              a(1:3,1) = a(1:3,1)/xtmp1
              a(1:3,2) = a(1:3,2)/xtmp2
              a(1:3,3) = a(1:3,3)/xtmp3
              
c             velocity gradients
              uxt(1:3,1) = ux(j,k,i,1:3)    !du/dx,du/dy,du/dz
              uxt(1:3,2) = ux(j,k,i,4:6)
              uxt(1:3,3) = ux(j,k,i,7:9)
              s(1,1) = uxt(1,1)
              s(2,2) = uxt(2,2)
              s(3,3) = uxt(3,3)
              s(2,1) = 0.5*(uxt(1,2) + uxt(2,1))
              s(3,1) = 0.5*(uxt(1,3) + uxt(3,1))
              s(3,2) = 0.5*(uxt(2,3) + uxt(3,2))
              s(1,2) = s(2,1)
              s(1,3) = s(3,1)
              s(2,3) = s(3,2)
              
c             numerator
              dds = 0.0
              do kk = 1,3
                do m = 1,3
                do n = 1,3
                do ii= 1,3
                do jj= 1,3
                  dds = dds + delta(kk)*a(ii,kk)*uxt(ii,m)
     .                       *delta(kk)*a(jj,kk)*uxt(jj,n)*s(m,n)
                end do
                end do        
                end do
                end do
              end do
              
c             denominator
              dd  = uxt(1,1)*uxt(1,1)
     .            + uxt(2,1)*uxt(2,1)
     .            + uxt(3,1)*uxt(3,1)
     .            + uxt(1,2)*uxt(1,2)
     .            + uxt(2,2)*uxt(2,2)
     .            + uxt(3,2)*uxt(3,2)
     .            + uxt(1,3)*uxt(1,3)
     .            + uxt(2,3)*uxt(2,3)
     .            + uxt(3,3)*uxt(3,3)
              s0  = s(1,1)*s(1,1) + s(2,2)*s(2,2) + s(3,3)*s(3,3) +
     .            2.*s(1,2)*s(1,2) + 2.*s(1,3)*s(1,3) + 2.*s(2,3)*s(2,3)
              dd  = dd*sqrt(2.*s0)
              !dd  = dd*sqrt(dd)  
              
c             length scale
c             xls = sqrt(abs(-dds)/ccmax(dd,1.e-10))
              xls = -dds/ccmax(dd,1.e-10)
c             spatial average desired as in Shur et al. (2015)
              xlesfr(j,k,i) = cdesamd*cdesamd*xls
c       
            enddo
          enddo
        enddo
c
      else if(ides.eq.13) then
c     simplified/modified implementation of SLA length scale by Xiao, Maochao
c     do not incorporate face/body diagnals, so 1/sqrt(3) is left out
c     farfield treatment done outside FKH, together with boundary layer treatment
c     refer to Shur, 2015 and Guseva, 2017 for the original version
        do i=1,idim-1
          do j=1,jdim-1
            do k=1,kdim-1
c
              delta(1) = 2.*vol(j,k,i)/(si(j,k,i,4)+si(j,k,i+1,4))
              delta(2) = 2.*vol(j,k,i)/(sj(j,k,i,4)+sj(j+1,k,i,4))
              delta(3) = 2.*vol(j,k,i)/(sk(j,k,i,4)+sk(j,k+1,i,4))
              if( i2d .eq. 1 .or. iaxi2planeturb .eq. 1 ) then
                  delta(1) = 0.0  !Modification to allow 2d
              end if
              
c             edge vectors
              a(1:3,1) = si(j,k,i,1:3)+si(j,k,i+1,1:3)
              a(1:3,2) = sj(j,k,i,1:3)+sj(j+1,k,i,1:3)
              a(1:3,3) = sk(j,k,i,1:3)+sk(j,k+1,i,1:3)
              xtmp1 = sqrt(a(1,1)*a(1,1)+a(2,1)*a(2,1)+a(3,1)*a(3,1))
              xtmp2 = sqrt(a(1,2)*a(1,2)+a(2,2)*a(2,2)+a(3,2)*a(3,2))
              xtmp3 = sqrt(a(1,3)*a(1,3)+a(2,3)*a(2,3)+a(3,3)*a(3,3))
              a(1:3,1) = a(1:3,1)/xtmp1*delta(1)
              a(1:3,2) = a(1:3,2)/xtmp2*delta(2)
              a(1:3,3) = a(1:3,3)/xtmp3*delta(3)
              
c             velocity gradients
              uxt(1:3,1) = ux(j,k,i,1:3)    !du/dx,du/dy,du/dz
              uxt(1:3,2) = ux(j,k,i,4:6)
              uxt(1:3,3) = ux(j,k,i,7:9)
              s(1,1) = uxt(1,1)
              s(2,2) = uxt(2,2)
              s(3,3) = uxt(3,3)
              s(2,1) = 0.5*(uxt(1,2) + uxt(2,1)) !0.5*(dvdx+dudy)
              s(3,1) = 0.5*(uxt(1,3) + uxt(3,1))
              s(3,2) = 0.5*(uxt(2,3) + uxt(3,2))
              s(1,2) = s(2,1)
              s(1,3) = s(3,1)
              s(2,3) = s(3,2)
              tr_s  = s(1,1) + s(2,2) + s(3,3)
              tr_s2 = s(1,1)*s(1,1) + s(2,2)*s(2,2) + s(3,3)*s(3,3)
     .              + 2.*s(1,2)*s(1,2) + 2.*s(2,3)*s(2,3) 
     .              + 2.*s(1,3)*s(1,3)
              w(1) = uxt(2,3) - uxt(3,2)
              w(2) = uxt(3,1) - uxt(1,3)
              w(3) = uxt(1,2) - uxt(2,1)
              w_mod = sqrt(w(1)*w(1)+w(2)*w(2)+w(3)*w(3))
              b(1:3) = w(1:3)/w_mod
c
              xls = 0.
              do kk = 1,3
                call cross_prod(a(1,kk),b(1),c(1))
                xtmp = c(1)*c(1) + c(2)*c(2) + c(3)*c(3)
                xls = ccmax(xls,xtmp)
              enddo
              xls = sqrt(xls)
c
c             VTM
              sw_1 = s(1,1)*w(1) + s(1,2)*w(2) + s(1,3)*w(3)
              sw_2 = s(2,1)*w(1) + s(2,2)*w(2) + s(2,3)*w(3)
              sw_3 = s(3,1)*w(1) + s(3,2)*w(2) + s(3,3)*w(3)
              sww_1 = sw_2*w(3) - sw_3*w(2)
              sww_2 = sw_3*w(1) - sw_1*w(3)
              sww_3 = sw_1*w(2) - sw_2*w(1)
              sww_mod = sqrt(sww_1*sww_1+sww_2*sww_2+sww_3*sww_3)
              vtm = 2.45*sww_mod/(w_mod*w_mod*sqrt(3.*tr_s2-tr_s*tr_s))
              fkh = ccmax(0.1,ccmin(1.,0.1+6.0*(vtm-0.15)))             
c
c             length scale
              xls = fkh*xls
              xlesfr(j,k,i) = cdes*xls
c       
            enddo
          enddo
        enddo
c
      end if
      
      return
      end
c
c
      subroutine cross_prod(v1,v2,v3)
      implicit none
      real, intent(in) ::v1(3), v2(3)
      real, intent(out) :: v3(3)
      
      v3(1) = v1(2)*v2(3) - v1(3)*v2(2)
      v3(2) = v1(3)*v2(1) - v1(1)*v2(3)
      v3(3) = v1(1)*v2(2) - v1(2)*v2(1)
      
      return
      end subroutine

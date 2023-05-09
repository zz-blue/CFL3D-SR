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
      subroutine fdelay(jdim,kdim,idim,q,sj,sk,si,vol,vist3d,smin,fnu,
     +              ux,iex3,vk,re,fdsav,felv)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Compute delayed function fd and elevating function fe
c     
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
      dimension q(jdim,kdim,idim,5),sj(jdim,kdim,idim-1,5),
     + sk(jdim,kdim,idim-1,5),si(jdim,kdim,idim,5),
     + vol(jdim,kdim,idim-1),vist3d(jdim,kdim,idim),
     + smin(jdim-1,kdim-1,idim-1)
      dimension fnu(0:jdim,0:kdim,0-iex3:idim+iex3)
      dimension ux(jdim-1,kdim-1,idim-1,9)
      dimension fdsav(jdim-1,kdim-1,idim-1),
     +          felv(jdim-1,kdim-1,idim-1)
c
      common /twod/ i2d
      common /axisym/ iaxi2plane,iaxi2planeturb,istrongturbdis,iforcev0
      common /des/ cdes,ides,cddes,cdesamd,lfe,cdt1,clim,lnwsc
c
      if(ides.eq.11 .or. ides.eq.12 .or. ides.eq.13) then
c
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
            dist = ccabs(smin(j,k,i))
            velterm= ux(j,k,i,1)*ux(j,k,i,1) + ux(j,k,i,2)*ux(j,k,i,2)
     .             + ux(j,k,i,3)*ux(j,k,i,3) + ux(j,k,i,4)*ux(j,k,i,4)
     .             + ux(j,k,i,5)*ux(j,k,i,5) + ux(j,k,i,6)*ux(j,k,i,6)
     .             + ux(j,k,i,7)*ux(j,k,i,7) + ux(j,k,i,8)*ux(j,k,i,8)
     .             + ux(j,k,i,9)*ux(j,k,i,9)
            rdt = vist3d(j,k,i)/(q(j,k,i,1)*      
     .           sqrt(velterm)*vk*vk*dist*dist*re)
            fdt = 1.0-cctanh((cdt1*rdt)*(cdt1*rdt)*(cdt1*rdt))
            rb  = 0.25-dist/deltam
            fb  = ccmincr(2.*exp(-9.*rb*rb),1.0)
            fd  = ccmax(1.0-fdt,fb)   !delayed function
            
            if( lfe.eq.1 ) then
                rdl = fnu(j,k,i)/(q(j,k,i,1)*
     .                sqrt(velterm)*vk*vk*dist*dist*re)
                fl  = cctanh((25.*rdl)**10)
                ft  = cctanh((3.4969*rdt)**3)
                fe2 = 1.0-ccmax(ft,fl)
                if( rb .ge. 0.0) then
                    fe1 = 2.*exp(-11.09*rb*rb)
                else
                    fe1 = 2.*exp(-9.*rb*rb)
                end if
                fe  = fe2*ccmaxcr(fe1-1.0,0.0)  !elevating function
           else
                fe  = 0.0
           end if
            felv(j,k,i) = fe
            fdsav(j,k,i)= fd
          end do
        end do
      end do
c
      else
c     please add other delayed/elevating functions
c
      end if
      
      return
      end

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
      subroutine plot3d_avg(jdim,kdim,idim,i1,i2,i3,j1,j2,j3,k1,k2,k3,
     .                     q,qi0,qj0,qk0,x,y,z,xw,blank2,blank,xg,iflag,
     .                     vist3d,iover,nblk,nmap,bcj,bck,bci,
     .                     vj0,vk0,vi0,ifunc,iplot,jdw,kdw,idw,
     .                     nplots,jdimg,kdimg,idimg,nblcg,jsg,ksg,isg,
     .                     jeg,keg,ieg,ninter,iindex,intmax,nsub1,
     .                     maxxe,nblkk,nbli,limblk,isva,nblon,mxbli,
     .                     thetay,maxbl,maxgr,myid,myhost,mycomm,
     .                     mblk2nd,inpl3d,nblock,nblkpt,xv,
     .                     sj,sk,si,vol,nummem,ifuncdim,qavg,q2avg,
     .                     qsavg,qs2avg,vdavg,vd2avg,nvdsp,
     .                     vsj0,vsjdim,vsk0,vskdim,vsi0,vsidim)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Write the output file at the grid points in function (PLOT3D)
c     format and print solution data
c
c     outputs grid/solution in single precision for use with FAST/PLOT3D
c***********************************************************************
c
#if defined ADP_OFF
#   ifdef CMPLX
#     ifdef DBLE_PRECSN
      implicit complex*8(a-h,o-z)
#     else
      implicit complex(a-h,o-z)
#     endif
#   else
#     ifdef DBLE_PRECSN
      implicit real*8 (a-h,o-z)
#     endif
#   endif
#else
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
#endif
c
#if defined DIST_MPI
#     include "mpif.h"
      dimension istat(MPI_STATUS_SIZE)
#   ifdef P3D_SINGLE
#     define MY_MPI_REAL MPI_REAL
#   else
#     define MY_MPI_REAL MPI_DOUBLE_PRECISION
#   endif
#endif
#if defined P3D_SINGLE
      real*4    xw(jdw,kdw,idw,2*ifuncdim+3),xg(jdw,kdw,idw,4)
      real*4    xmachw,alphww,reuew,timew
#else
      real xw(jdw,kdw,idw,2*ifuncdim+3),xg(jdw,kdw,idw,4)
#endif
c
      dimension xv(jdw,kdw,idw,2*ifuncdim+3)
      dimension x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim)
      dimension q(jdim,kdim,idim,5), qi0(jdim,kdim,5,4),
     .          qj0(kdim,idim-1,5,4),qk0(jdim,idim-1,5,4)
      dimension bcj(kdim,idim-1,2),bck(jdim,idim-1,2),bci(jdim,kdim,2)
      dimension vist3d(jdim,kdim,idim)
      dimension vj0(kdim,idim-1,1,4),vk0(jdim,idim-1,1,4),
     .          vi0(jdim,kdim,1,4)
      dimension blank(jdim,kdim,idim),blank2(jdim,kdim,idim,2)
      dimension nmap(maxbl),jdimg(maxbl),kdimg(maxbl),idimg(maxbl),
     .          nblcg(maxbl),jsg(maxbl),ksg(maxbl),isg(maxbl),
     .          jeg(maxbl),keg(maxbl),ieg(maxbl),
     .          iindex(intmax,6*nsub1+9),nblkpt(maxxe),nblkk(2,mxbli),
     .          limblk(2,6,mxbli),isva(2,2,mxbli),nblon(mxbli), 
     .          inpl3d(nplots,11),thetay(maxbl),mblk2nd(maxbl)
      dimension sk(jdim,kdim,idim-1,5),sj(jdim,kdim,idim-1,5),
     .          si(jdim,kdim,idim,5),vol(jdim,kdim,idim-1)
      dimension turre(jdim,kdim,idim,nummem)
      dimension tj0(kdim,idim-1,nummem,4),tk0(jdim,idim-1,nummem,4),
     .          ti0(jdim,kdim,nummem,4)
      dimension vdsp(jdim,kdim,idim,nvdsp),vdj0(kdim,idim-1,nvdsp,4),
     .          vdk0(jdim,idim-1,nvdsp,4),vdi0(jdim,kdim,nvdsp,4)
      dimension qavg(jdim,kdim,idim,5),q2avg(jdim,kdim,idim,5),
     .          qsavg(jdim,kdim,idim,nummem+7),
     .          qs2avg(jdim,kdim,idim,nummem+7),
     .          vdavg(jdim,kdim,idim,nvdsp),
     .          vd2avg(jdim,kdim,idim,nvdsp)
      dimension vsj0(kdim,idim,ifuncdim,2),vsjdim(kdim,idim,ifuncdim,2),
     .          vsk0(jdim,idim,ifuncdim,2),vskdim(jdim,idim,ifuncdim,2),
     .          vsi0(jdim,kdim,ifuncdim,2),vsidim(jdim,kdim,ifuncdim,2)
      allocatable :: xvist3d(:,:,:),xturre(:,:,:,:),xvdsp(:,:,:,:)
c
      common /bin/ ibin,iblnk,iblnkfr,ip3dgrad
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /fluid2/ pr,prt,cbar
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /reyue/ reue,tinf,ivisc(3)
      common /twod/ i2d
      common /unst/ time,cfltau,ntstep,ita,iunst,cfltau0,cfltauMax
      common /wallfun/ iwf(3)
      common /igrdtyp/ ip3dgrd,ialph
      common /moov/movie,nframes,icall1,lhdr,icoarsemovie,i2dmovie,
     .        icallsf
      common /conversion/ radtodeg
      common /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),
     .        tur10(7)
      common /avgdata/ xnumavg,iteravg,xnumavg2,ipertavg,npertavg,
     .      iclcd,isubit_r,icallavg
c
#if defined DIST_MPI
c     set baseline tag values
c
      ioffset = maxbl
      itag_grid = 1
      itag_q    = itag_grid + ioffset
#endif
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
c
c     initialize xw, xv and xg arrays
c
      jw = (j2-j1)/j3 + 1
      kw = (k2-k1)/k3 + 1
      iw = (i2-i1)/i3 + 1
      do j=1,jw
         do k=1,kw
            do i=1,iw
               do l= 1,5
               xw(j,k,i,l) = 0.
               xv(j,k,i,l) = 0.
               end do
               do l= 1,4
               xg(j,k,i,l) = 0.
               end do
            end do
         end do
      end do
c
c     assign single precision scalars
c
      alphaw = radtodeg*(alpha+thetay(nblk))
      xmachw = xmach
      alphww = alphaw
      reuew  = reue
      timew  = time
c
c***********************************************************************
c     load grid into first 3 locations of the single precision xg array
c***********************************************************************
c
#if defined DIST_MPI
      if (mblk2nd(nblk).eq.myid) then
c
#endif
      iw = 0
      do 2000 i=i1,i2,i3
      iw = iw + 1
      jw = 0
      do 2000 j=j1,j2,j3
      jw = jw + 1
      kw = 0
      do 2000 k=k1,k2,k3
      kw = kw + 1
      xg(jw,kw,iw,1) = x(j,k,i)
      xg(jw,kw,iw,2) = y(j,k,i)
      xg(jw,kw,iw,3) = z(j,k,i)
 2000 continue
c
#if defined DIST_MPI
      end if
c
#endif
c***********************************************************************
c     set iblank (blank2) array
c***********************************************************************
c
c     currently don't need blanking data for printout option (iflag=2),
c     so only compute blank2 array for plot3d output option
c     also, don't compute blanking info if iblnk = 0
c
      if (iflag.eq.1 .and. iblnk.gt.0) then
c
#if defined DIST_MPI
      if (mblk2nd(nblk).eq.myid) then
c
#endif
c     assign default iblank (blank2) array 
c
      do 2049 i=1,idim
      do 2049 j=1,jdim
      do 2049 k=1,kdim
      blank2(j,k,i,1) = 1.
      blank2(j,k,i,2) = 1.
 2049 continue
c
c     zero out edges and corners if desired (izero > 0)
c
      izero = 0
c     
      if (izero.gt.0) then
         do 2050 i=1,idim
         do 2050 j=1,jdim,jdim1
         do 2050 k=1,kdim,kdim1
         blank2(j,k,i,1) = 0.
 2050    continue
c
         if (idim.gt.2) then
            do 2051 i=1,idim,idim1
            do 2151 j=1,jdim
            blank2(j,1,i,1)    = 0.
            blank2(j,kdim,i,1) = 0.
 2151       continue
            do 2251 k=1,kdim
            blank2(1,k,i,1)    = 0.
            blank2(jdim,k,i,1) = 0.
 2251       continue
 2051       continue
         else
            do 2052 j=1,jdim
            do 2052 k=1,kdim
            blank2(j,k,idim,1) = 0.
 2052       continue
         end if
c
         if (jdim.eq.2) then
            do 2054 i=1,idim
            do 2054 k=1,kdim
            blank2(jdim,k,i,1) = 0.
 2054       continue
         end if
c
         if (kdim.eq.2) then
            do 2056 i=1,idim
            do 2056 k=1,kdim
            blank2(j,kdim,i,1) = 0.
 2056       continue
         end if
      end if
c
c     solid surface iblank (blank2) values (iblank=2 for solid surface)
c
      j11 = 1
      j22 = jdim1
      if(jdim.eq.2) then
        j11 = 1
        j22 = 1
      end if
      k11 = 1
      k22 = kdim1
      if(kdim.eq.2) then
        k11 = 1
        k22 = 1
      end if
      i11 = 1
      i22 = idim1
      if(idim.eq.2) then
        i11 = 1
        i22 = 1
      end if
      i = 1
      do 1190 nnn=1,2
      do 1191 j=j11,j22
      do 1191 k=k11,k22
      kk = min(k+1,k22)
      jj = min(j+1,j22)
      blank2(j,k,i,1) = 1. + ccmax4( bci(j,k,nnn),  bci(jj,k,nnn),
     .                               bci(j,kk,nnn), bci(jj,kk,nnn) )
 1191 continue
      i = idim
 1190 continue
      j = 1
      do 2190 nnn=1,2
      do 2191 i=i11,i22
      do 2191 k=k11,k22
      kk = min(k+1,k22)
      ii = min(i+1,i22)
      blank2(j,k,i,1) = 1. + ccmax4( bcj(k,i,nnn),  bcj(kk,i,nnn),
     .                               bcj(k,ii,nnn), bcj(kk,ii,nnn) )
 2191 continue
      j = jdim
 2190 continue
      k = 1
      do 3190 nnn=1,2
      do 3191 j=j11,j22
      do 3191 i=i11,i22
      jj = min(j+1,j22)
      ii = min(i+1,i22)
      blank2(j,k,i,1) = 1. + ccmax4( bck(j,i,nnn),  bck(jj,i,nnn),
     .                               bck(j,ii,nnn), bck(jj,ii,nnn) )
 3191 continue
      k = kdim
 3190 continue
c
c     iblank (blank2) array for patch interface boundaries
c  
c     Set blank2 array on block edges from nblkpt array.  Nblkpt array
c     corresponds to cell center grid.  For ambiguous points which arise
c     when going from cell-center data to grid point data (e.g. along block
c     edges, corners, or for interior points on a patch interface for which 
c     surrounding cell-centers lie in different blocks), blank2 value 
c     remains equal to previously set values.
c
      if (abs(ninter).gt.0) then
         do 1600 inter=1,abs(ninter)
         lmax1  = iindex(inter,1)
         nbl    = iindex(inter,lmax1+2)
         if (nbl.ne.nblk) go to 1600
         lst    = iindex(inter,2*lmax1+5)
         lcoord = iindex(inter,2*lmax1+3)/10
         lend   = iindex(inter,2*lmax1+3)-lcoord*10
         j21    = iindex(inter,2*lmax1+6)
         j22    = iindex(inter,2*lmax1+7)
         k21    = iindex(inter,2*lmax1+8)
         k22    = iindex(inter,2*lmax1+9)
c
         if (lcoord.eq.1) then
            if (lend.eq.1) i = 1
            if (lend.eq.2) i = idim
            if (jdim.gt.2 .and. kdim.gt.2) then
               do 1610 j=j21+1,j22-1
               do 1610 k=k21+1,k22-1
c              connecting block numbers of the four cell-centered points 
c              surrounding the current grid point
               ll = lst + (j22-j21)*(k-1-k21) + (j-1-j21)
               mblk1 = iindex(inter,nblkpt(ll)+1)
               ll = lst + (j22-j21)*(k-k21-1) + (j-j21)
               mblk2 = iindex(inter,nblkpt(ll)+1)
               ll = lst + (j22-j21)*(k-k21) + (j-j21)
               mblk3 = iindex(inter,nblkpt(ll)+1)
               ll = lst + (j22-j21)*(k-k21) + (j-j21-1)
               mblk4 = iindex(inter,nblkpt(ll)+1)
               if (mblk1.eq.mblk2 .and. mblk1.eq.mblk3 .and.
     .             mblk1.eq.mblk4) blank2(j,k,i,1) = -float(nmap(mblk1))
 1610          continue
            else if (jdim.gt.2) then
c              connecting block numbers of the two cell-centered points
c              surrounding the current grid point
               k = k21
               do 1611 j=j21+1,j22-1
               ll = lst + (j22-j21)*(k-k21) + (j-j21)
               mblk1 = iindex(inter,nblkpt(ll)+1)
               ll = lst + (j22-j21)*(k-k21) + (j-j21-1)
               mblk2 = iindex(inter,nblkpt(ll)+1)
               if (mblk1.eq.mblk2) blank2(j,k,i,1) = -float(nmap(mblk1))
 1611          continue
            else if (kdim.gt.2) then
               j = j21
               do 1612 k=k21+1,k22-1
               ll = lst + (j22-j21)*(k-k21-1) + (j-j21)
               mblk1 = iindex(inter,nblkpt(ll)+1)
               ll = lst + (j22-j21)*(k-k21) + (j-j21)
               mblk2 = iindex(inter,nblkpt(ll)+1)
               if (mblk1.eq.mblk2) blank2(j,k,i,1) = -float(nmap(mblk1))
 1612          continue
            end if
         end if
c
         if (lcoord.eq.2) then
            if (lend.eq.1) j = 1
            if (lend.eq.2) j = jdim
            if (idim.gt.2 .and .kdim.gt.2) then
               do 1620 i=k21+1,k22-1
               do 1620 k=j21+1,j22-1
c              connecting block numbers of the four cell-centered points 
c              surrounding the current grid point
               ll = lst + (j22-j21)*(i-1-k21) + (k-1-j21)
               mblk1 = iindex(inter,nblkpt(ll)+1)
               ll = lst + (j22-j21)*(i-k21-1) + (k-j21)
               mblk2 = iindex(inter,nblkpt(ll)+1)
               ll = lst + (j22-j21)*(i-k21) + (k-j21)
               mblk3 = iindex(inter,nblkpt(ll)+1)
               ll = lst + (j22-j21)*(i-k21) + (k-j21-1)
               mblk4 = iindex(inter,nblkpt(ll)+1)
               if (mblk1.eq.mblk2 .and. mblk1.eq.mblk3 .and.
     .             mblk1.eq.mblk4) blank2(j,k,i,1) = -float(nmap(mblk1))
 1620          continue
            else if (idim.gt.2) then
               k = j21
               do 1621 i=k21+1,k22-1
               ll = lst + (j22-j21)*(i-k21-1) + (k-j21)
               mblk1 = iindex(inter,nblkpt(ll)+1)
               ll = lst + (j22-j21)*(i-k21) + (k-j21)
               mblk2 = iindex(inter,nblkpt(ll)+1)
               if (mblk1.eq.mblk2) blank2(j,k,i,1) = -float(nmap(mblk1))
 1621          continue
            else if (kdim.gt.2) then
               i = k21
               do 1622 k=j21+1,j22-1
               ll = lst + (j22-j21)*(i-k21) + (k-j21)
               mblk1 = iindex(inter,nblkpt(ll)+1)
               ll = lst + (j22-j21)*(i-k21) + (k-j21-1)
               mblk2 = iindex(inter,nblkpt(ll)+1)
               if (mblk1.eq.mblk2) blank2(j,k,i,1) = -float(nmap(mblk1))
 1622          continue
            end if
         end if
c
         if (lcoord.eq.3) then
            if (lend.eq.1) k = 1
            if (lend.eq.2) k = kdim
            if (idim.gt.2 .and. jdim.gt.2) then
               do 1630 i=k21+1,k22-1
               do 1630 j=j21+1,j22-1
c              connecting block numbers of the four cell-centered points 
c              surrounding the current grid point
               ll = lst + (j22-j21)*(i-1-k21) + (j-1-j21)
               mblk1 = iindex(inter,nblkpt(ll)+1)
               ll = lst + (j22-j21)*(i-k21-1) + (j-j21)
               mblk2 = iindex(inter,nblkpt(ll)+1)
               ll = lst + (j22-j21)*(i-k21) + (j-j21)
               mblk3 = iindex(inter,nblkpt(ll)+1)
               ll = lst + (j22-j21)*(i-k21) + (j-j21-1)
               mblk4 = iindex(inter,nblkpt(ll)+1)
               if (mblk1.eq.mblk2 .and. mblk1.eq.mblk3 .and.
     .             mblk1.eq.mblk4) blank2(j,k,i,1) = -float(nmap(mblk1))
 1630          continue
            else if (idim.gt.2) then
               j = j21
               do 1631 i=k21+1,k22-1
               ll = lst + (j22-j21)*(i-k21-1) + (j-j21)
               mblk1 = iindex(inter,nblkpt(ll)+1)
               ll = lst + (j22-j21)*(i-k21) + (j-j21)
               mblk2 = iindex(inter,nblkpt(ll)+1)
               if (mblk1.eq.mblk2) blank2(j,k,i,1) = -float(nmap(mblk1))
 1631          continue
            else if (jdim.gt.2) then
               i = k21
               do 1632 j=j21+1,j22-1
               ll = lst + (j22-j21)*(i-k21) + (j-j21)
               mblk1 = iindex(inter,nblkpt(ll)+1)
               ll = lst + (j22-j21)*(i-k21) + (j-j21-1)
               mblk2 = iindex(inter,nblkpt(ll)+1)
               if (mblk1.eq.mblk2) blank2(j,k,i,1) = -float(nmap(mblk1))
 1632          continue
            end if
         end if
c
 1600    continue
      end if
c
c     iblank (blank2) array for 1:1 interface boundaries
c
      if(nbli.gt.0) then
        do 100 n=1,abs(nbli)
        if(nblon(n).ge.0) then
          if(nblk.eq.nblkk(1,n) .or. nblk.eq.nblkk(2,n)) then
            it = 1
            ir = 2
            if(nblk.eq.nblkk(2,n)) then
              it = 2
              ir = 1
            end if
c
c           allow for 1-1 blocking in same grid
c
            itime = 1
            if (nblkk(1,n).eq.nblkk(2,n)) itime = 2
            do 101 iti = 1,itime
            if (iti.gt.1) then
               it = 1
               ir = 2
            end if
c
            is = limblk(it,1,n)
            ie = limblk(it,4,n)
            js = limblk(it,2,n)
            je = limblk(it,5,n)
            ks = limblk(it,3,n)
            ke = limblk(it,6,n)
c
c           cell center indices ---> grid point indices
c
            if(isva(it,1,n)+isva(it,2,n).eq.5) then
c             i = constant interface
              if(js.gt.je) js = js+1
              if(js.lt.je) je = je+1
              if(ks.gt.ke) ks = ks+1
              if(ks.lt.ke) ke = ke+1
c             2d cases
              if(jdim.eq.2) then
                js = 1
                je = 2
              end if
              if(kdim.eq.2) then
                ks = 1
                ke = 2
              end if
            end if
            if(isva(it,1,n)+isva(it,2,n).eq.4) then
c             j = constant interface
              if(is.gt.ie) is = is+1
              if(is.lt.ie) ie = ie+1
              if(ks.gt.ke) ks = ks+1
              if(ks.lt.ke) ke = ke+1
c             2d cases
              if(idim.eq.2) then
                is = 1
                ie = 2
              end if
              if(kdim.eq.2) then
                ks = 1
                ke = 2
              end if
            end if
            if(isva(it,1,n)+isva(it,2,n).eq.3) then
c             k = constant interface
              if(js.gt.je) js = js+1
              if(js.lt.je) je = je+1
              if(is.gt.ie) is = is+1
              if(is.lt.ie) ie = ie+1
c             2d cases
              if(jdim.eq.2) then
                js = 1
                je = 2
              end if
              if(idim.eq.2) then
                is = 1
                ie = 2
              end if
            end if
c
            is1 = min(is,ie)
            ie1 = max(is,ie)
            js1 = min(js,je)
            je1 = max(js,je)
            ks1 = min(ks,ke)
            ke1 = max(ks,ke)
c
            do 110 i = is1,ie1
            do 110 j = js1,je1
            do 110 k = ks1,ke1
            blank2(j,k,i,1) = -float(nmap(nblkk(ir,n)))
  110       continue
  101       continue
          end if
        end if
  100   continue
      end if
c
c     iblank (blank2) array for embedded grids - the underlying
c     coarse grid areas are blanked out if the parameter ibembed > 0
c
      ibembed = 1
c
      if (ibembed.gt.0) then
         do 7500 nblc=1,nblock
         if (nblk.eq.nblc) go to 7500
         nblcc    = nblcg(nblc)
         if (nblcc.eq.nblk) then
            js = jsg(nblc)
            if (js.lt.jdimg(nblcc) .and. js.gt.1) js = js + 1
            ks = ksg(nblc)
            if (ks.lt.kdimg(nblcc) .and. ks.gt.1) ks = ks + 1
            is = isg(nblc)
            if (is.lt.idimg(nblcc) .and. is.gt.1) is = is + 1
            je = jeg(nblc)
            if (je.gt.2 .and. je.lt.jdimg(nblcc)) je = je - 1
            ke = keg(nblc)
            if (ke.gt.2 .and. ke.lt.kdimg(nblcc)) ke = ke - 1
            ie = ieg(nblc)
            if (ie.gt.2 .and. ie.lt.idimg(nblcc)) ie = ie - 1
            do 7501 i=is,ie
            do 7501 j=js,je
            do 7501 k=ks,ke
            blank2(j,k,i,1) = 0.
 7501       continue
         end if
 7500    continue
      end if
c
c     iblank (blank2) array for overlapped grids
c
      if (iover.eq.1) then
c
c        interior of faces (i=1 and i=idim)
c
         do 413 j=2,jdim1
         do 413 k=2,kdim1
         blank2(j,k,1,2)    = ccmin4(blank(j,k,1),
     .                               blank(j-1,k,1),
     .                               blank(j,k-1,1),
     .                               blank(j-1,k-1,1))
         blank2(j,k,idim,2) = ccmin4(blank(j,k,idim1),  
     .                               blank(j-1,k,idim1),
     .                               blank(j,k-1,idim1),
     .                               blank(j-1,k-1,idim1))
  413    continue
c
c        edges and corners (i=1 and i=idim)
c
         do 414 m=1,2
         k  = 1
         kk = 1
         if (m.eq.2) then
            k  = kdim
            kk = kdim1
         end if
         do 414 j=2,jdim1
         blank2(j,k,1,2)    = ccmin(blank(j,kk,1),
     .                              blank(j-1,kk,1))
         blank2(j,k,idim,2) = ccmin(blank(j,kk,idim1),
     .                              blank(j-1,kk,idim1))
  414    continue
         do 415 m=1,2
         j  = 1
         jj = 1
         if (m.eq.2) then
            j  = jdim
            jj = jdim1
         end if
         do 415 k=2,kdim1
         blank2(j,k,1,2)    = ccmin(blank(jj,k,1),
     .                              blank(jj,k-1,1))
         blank2(j,k,idim,2) = ccmin(blank(jj,k,idim1),
     .                              blank(jj,k-1,idim1))
  415    continue
         blank2(jdim,kdim,1,2)    = blank2(jdim-1,kdim-1,1,2)
         blank2(1,1,1,2)          = blank2(2,2,1,2)
         blank2(jdim,kdim,idim,2) = blank2(jdim-1,kdim-1,idim,2)
         blank2(1,1,idim,2)       = blank2(2,2,idim,2)
c
c        interior of faces (j=1 and j=jdim)
c
         do 513 i=2,idim1
         do 513 k=2,kdim1
         blank2(1,k,i,2)    = ccmin4(blank(1,k,i),
     .                               blank(1,k,i-1),
     .                               blank(1,k-1,i),
     .                               blank(1,k-1,i-1))
         blank2(jdim,k,i,2) = ccmin4(blank(jdim1,k,i),
     .                               blank(jdim1,k,i-1),
     .                               blank(jdim1,k-1,i),
     .                               blank(jdim1,k-1,i-1))
  513    continue
c
c        edges and corners (j=1 and j=jdim)
c  
         do 514 m=1,2
         k  = 1
         kk = 1
         if (m.eq.2) then
            k  = kdim
            kk = kdim1
         end if
         do 514 i=2,idim1
         blank2(1,k,i,2)    = ccmin(blank(1,kk,i),
     .                              blank(1,kk,i-1))
         blank2(jdim,k,i,2) = ccmin(blank(jdim1,kk,i),
     .                              blank(jdim1,kk,i-1))
  514    continue
         do 515 m=1,2
         i  = 1
         ii = 1
         if (m.eq.2) then
            i  = idim
            ii = idim1
         end if
         do 515 k=2,kdim1
         blank2(1,k,i,2)    = ccmin(blank(1,k,ii),
     .                              blank(1,k-1,ii))
         blank2(jdim,k,i,2) = ccmin(blank(jdim1,k,ii),
     .                              blank(jdim1,k-1,ii))
  515    continue
         blank2(1,kdim,idim,2)    = blank2(1,kdim-1,idim-1,2)
         blank2(1,1,1,2)          = blank2(1,2,2,2)
         blank2(jdim,kdim,idim,2) = blank2(jdim,kdim-1,idim-1,2)
         blank2(jdim,1,1,2)       = blank2(jdim,2,2,2)
c
c        interior of faces (k=1 and k=kdim)
c
         do 613 i=2,idim1
         do 613 j=2,jdim1
         blank2(j,1,i,2)    = ccmin4(blank(j,1,i),
     .                               blank(j,1,i-1),
     .                               blank(j-1,1,i),
     .                               blank(j-1,1,i-1))
         blank2(j,kdim,i,2) = ccmin4(blank(j,kdim1,i),
     .                               blank(j,kdim1,i-1),
     .                               blank(j-1,kdim1,i),
     .                               blank(j-1,kdim1,i-1))
  613    continue
c
c        edges and corners (k=1 and k=kdim)
c
         do 614 m=1,2
         j  = 1
         jj = 1
         if (m.eq.2) then
            j  = jdim
            jj = jdim1
         end if
         do 614 i=2,idim1
         blank2(j,1,i,2)    = ccmin(blank(jj,1,i),
     .                              blank(jj,1,i-1))
         blank2(j,kdim,i,2) = ccmin(blank(jj,kdim1,i),
     .                              blank(jj,kdim1,i-1))
  614    continue
         do 615 m=1,2
         i  = 1
         ii = 1
         if (m.eq.2) then
            i  = idim
            ii = idim1
         end if
         do 615 j=2,jdim1
         blank2(j,1,i,2)    = ccmin(blank(j,1,ii),
     .                              blank(j-1,1,ii))
         blank2(j,kdim,i,2) = ccmin(blank(j,kdim1,ii),
     .                              blank(j-1,kdim1,ii))
  615   continue
         blank2(jdim,1,idim,2)    = blank2(jdim-1,1,idim-1,2)
         blank2(1,1,1,2)          = blank2(2,1,2,2)
         blank2(jdim,kdim,idim,2) = blank2(jdim-1,kdim,idim-1,2)
         blank2(1,kdim,1,2)       = blank2(2,kdim,2,2)
c
c       interior cells
c
         do 713 i=2,idim1
         do 713 j=2,jdim1
         do 713 k=2,kdim1
         blank2(j,k,i,2) = ccmin8(blank(j,k,i),
     .                            blank(j-1,k,i),
     .                            blank(j,k-1,i),
     .                            blank(j-1,k-1,i),
     .                            blank(j,k,i-1),
     .                            blank(j-1,k,i-1),
     .                            blank(j,k-1,i-1),
     .                            blank(j-1,k-1,i-1))
  713    continue
c
c        combining topology and overset grid blankings
c
         do 813 i=1,idim
         do 813 j=1,jdim
         do 813 k=1,kdim
         blank2(j,k,i,1) = ccmin(blank2(j,k,i,1),blank2(j,k,i,2))
  813    continue
c
      end if
c
#if defined DIST_MPI
      end if
c
#endif
      end if
c
      if (iflag.eq.1) then
c
c***********************************************************************
c      plot3d data
c***********************************************************************
c
      if (lhdr.gt.0 .and. myid.eq.myhost) then
         if (i2d .eq. 1) then
            write(11,'(''writing plot3d file for JDIM X KDIM ='',i5,
     .      '' x '',i5,'' grid'')') jdim,kdim
         else
            write(11,93)idim,jdim,kdim
   93       format(45h writing plot3d file for IDIM X JDIM X KDIM =,
     .      i5,3h x ,i5,3h x ,i5,5h grid)
         end if
         if (inpl3d(iplot,2).eq.0) then
            write(11,'(''   plot3dg file is an xyz file'',
     .      '' at grid points'')')
            write(11,'(''   plot3dq file is a    q file'',
     .      '' at grid points'')')
#           ifdef CMPLX
            if (ip3dgrad .ne. 0) then 
               write(11,'(''   NOTE: contains dq/d(), rather than q'')')
               call gradinfo(delh,11)
            end if
#           endif
         end if
         if (inpl3d(iplot,2).eq.-1) then
            write(11,'(''   plot3dg file is an xyz file'',
     .      '' at grid points'')')
            write(11,'(''   plot3dq file is a    q file'',
     .      '' at grid points with off-surface data used when output'',
     .      '' is on body only'')')
#           ifdef CMPLX
            if (ip3dgrad .ne. 0) then
               write(11,'(''   NOTE: contains dq/d(), rather than q'')')
               call gradinfo(delh,11)
            end if
#           endif
         end if
         if (ifunc.eq.-2) then
            write(11,'(''   plot3dg file is an xyz file'',
     .      '' at grid points'')')
            write(11,'(''   plot3dq file is a function'',
     .      '' file at grid points'')')
#           ifdef CMPLX
            if (ip3dgrad .ne. 0) then
               write(11,'(''   NOTE: contains d()/d(), rather'',
     .                    '' than amut'')')
               call gradinfo(delh,11)
            end if
#           endif
         end if
         if (ifunc.eq.-5) then
            write(11,'(''   plot3dg file is an xyz file'',
     .      '' at grid points'')')
            write(11,'(''   plot3dq file is a function'',
     .      '' file of Cp at grid points'')')
#           ifdef CMPLX
            if (ip3dgrad .ne. 0) then 
               write(11,'(''   NOTE: contains dCp/d(), rather than'',
     .                    '' Cp'')')
               call gradinfo(delh,11)
            end if
#           endif
         end if
         if (ifunc.eq.-6) then
            write(11,'(''   plot3dg file is an xyz file'',
     .      '' at grid points'')')
            write(11,'(''   plot3dq file is a function'',
     .      '' file of P/Pinf at grid points'')')
#           ifdef CMPLX
            if (ip3dgrad .ne. 0) then 
               write(11,'(''   NOTE: contains d(P/Pinf)/d(), rather'',
     .                    '' than P/Pinf'')')
               call gradinfo(delh,11)
            end if
#           endif
         end if
         if (ifunc.eq.-4) then
            write(11,'(''   plot3dg file is an xyz file'',
     .      '' at grid points'')')
            write(11,'(''   plot3dq file is a function'',
     .      '' file of amut at grid points'')')
#           ifdef CMPLX
            if (ip3dgrad .ne. 0) then
               write(11,'(''   NOTE: contains d(amut)/d(), rather'',
     .                    '' than amut'')')
               call gradinfo(delh,11)
            end if
#           endif
         end if
         if (i2d .eq. 1) then
            if (iblnk .eq. 0) then
               write(11,2041)
 2041          format(3x,35hplot3d files to be read with /mgrid,
     .         14h/2d qualifiers)
            else
               write(11,2141)
 2141          format(3x,41hplot3d files to be read with /mgrid/blank,
     .         14h/2d qualifiers)
            end if
         else
            if (iblnk .eq. 0) then
               write(11,2042)
 2042          format(3x,35hplot3d files to be read with /mgrid,
     .         11h qualifiers)
            else
               write(11,2142)
 2142          format(3x,41hplot3d files to be read with /mgrid/blank,
     .         11h qualifiers)
            end if
         end if
      end if
c
#if defined DIST_MPI
      if (mblk2nd(nblk).eq.myid) then
c
#endif
c     load blank2 into 4th location of the single precision xg array
c
      iw = 0
      do 9000 i=i1,i2,i3
      iw = iw + 1
      jw = 0
      do 9000 j=j1,j2,j3
      jw = jw + 1
      kw = 0
      do 9000 k=k1,k2,k3
      kw = kw + 1
      xg(jw,kw,iw,4) = blank2(j,k,i,1)
 9000 continue
#if defined DIST_MPI
c
      end if
#endif
#if defined DIST_MPI
c
c     send/receive plot3d grid data
c
      jw = (j2-j1)/j3 + 1
      kw = (k2-k1)/k3 + 1
      iw = (i2-i1)/i3 + 1
c
      nng = jw*kw*iw*4
      nd_srce = mblk2nd(nblk)
      mytag = itag_grid + nblk
c
      if (mblk2nd(nblk).eq.myid) then
         call MPI_Send (xg, nng, MY_MPI_REAL, myhost, mytag,
     &                  mycomm, ierr)
      else if (myid.eq.myhost) then
         call MPI_Recv (xg, nng, MY_MPI_REAL, nd_srce,
     &                  mytag, mycomm, istat, ierr)
      end if
c
#endif
c
      if (myid.eq.myhost .and. icallavg.eq.0) then
c
c     ialph > 0 for a grid that was read in plot3d format with alpha measured
c               in the xy plane (TLNS3D convention)
c
      if (ialph.ne.0 .and. i2d.ne.1) then
         do i=1,iw
            do j=1,jw
               do k=1,kw
                  temp        = xg(j,k,i,2)
                  xg(j,k,i,2) = xg(j,k,i,3)
                  xg(j,k,i,3) = -temp
               end do
            end do
         end do
      end if
c
c     output grid
c
      if (ibin.eq.0) then
         if (i2d.eq.0) then
            if (iblnk.eq.0) then
               write(96,'(5e14.6)')
     .                   (((xg(j,k,i,1),i=1,iw),j=1,jw),k=1,kw),
     .                   (((xg(j,k,i,2),i=1,iw),j=1,jw),k=1,kw),
     .                   (((xg(j,k,i,3),i=1,iw),j=1,jw),k=1,kw)
            else
                 write(96,'(5e14.6)')  !used
     .                     (((xg(j,k,i,1),i=1,iw),j=1,jw),k=1,kw),
     .                     (((xg(j,k,i,2),i=1,iw),j=1,jw),k=1,kw),
     .                     (((xg(j,k,i,3),i=1,iw),j=1,jw),k=1,kw)
                 write(96,'(5i14)')
     .                 (((int(xg(j,k,i,4)),i=1,iw),j=1,jw),k=1,kw)
            end if
         else
            if (iblnk.eq.0) then
               write(96,'(5e14.6)')
     .                   (((xg(j,k,i,1),i=1,iw),j=1,jw),k=1,kw),
     .                   (((xg(j,k,i,3),i=1,iw),j=1,jw),k=1,kw)
            else
               write(96,'(5e14.6)')    !used
     .                   (((xg(j,k,i,1),i=1,iw),j=1,jw),k=1,kw),
     .                   (((xg(j,k,i,3),i=1,iw),j=1,jw),k=1,kw)
               write(96,'(5i14)')
     .                 (((int(xg(j,k,i,4)),i=1,iw),j=1,jw),k=1,kw)
            end if
         end if
      else
         if (i2d.eq.0) then
            if (iblnk.eq.0) then
               write(96)
     .                 (((xg(j,k,i,1),i=1,iw),j=1,jw),k=1,kw),
     .                 (((xg(j,k,i,2),i=1,iw),j=1,jw),k=1,kw),
     .                 (((xg(j,k,i,3),i=1,iw),j=1,jw),k=1,kw)
            else
               write(96)
     .                 (((xg(j,k,i,1),i=1,iw),j=1,jw),k=1,kw),
     .                 (((xg(j,k,i,2),i=1,iw),j=1,jw),k=1,kw),
     .                 (((xg(j,k,i,3),i=1,iw),j=1,jw),k=1,kw),
     .             (((int(xg(j,k,i,4)),i=1,iw),j=1,jw),k=1,kw)
            end if
         else
            if (iblnk.eq.0) then
               write(96)
     .                 (((xg(j,k,i,1),i=1,iw),j=1,jw),k=1,kw),
     .                 (((xg(j,k,i,3),i=1,iw),j=1,jw),k=1,kw) 
            else
               write(96)
     .                 (((xg(j,k,i,1),i=1,iw),j=1,jw),k=1,kw),
     .                 (((xg(j,k,i,3),i=1,iw),j=1,jw),k=1,kw),
     .             (((int(xg(j,k,i,4)),i=1,iw),j=1,jw),k=1,kw)
            end if
         end if
      end if
c
      end if
c
#if defined DIST_MPI
      if (mblk2nd(nblk).eq.myid) then
c
#endif
c     determine q values at grid points and load into xv array
c     jdw = jdim,kdw = kdim,idw = idim
c
      if(0.eq.1) then
c     no need to do this if loop for plot3d_avg
      jdw = (j2-j1)/j3 + 1
      kdw = (k2-k1)/k3 + 1
      idw = (i2-i1)/i3 + 1
c
c     ifunc...flag to determine whether a plot3d q file or plot3d function
c             file is written. function files contain only one variable to
c             be compatable with the current version of FAST
c           = 0 plot3d q file is written at grid points
c             (if inpl3d(iplot,2)=-1 then same, except off-body data used 
c             when output is on body only)
c           < 0 function file is written at grid points
c            -4 vist3d at grid points
c            -5 pressure coefficient (cp) at grid points
c            -6 p/pinf at grid points
c
      if (ifunc .eq. -4) then
        allocate(xvist3d(jdw,kdw,idw))
        call cctogp(jdim,kdim,idim,i1,i2,i3,j1,j2,j3,k1,k2,k3,vist3d,
     .              vi0,vj0,vk0,jdw,kdw,idw,xvist3d,1)
      else
      
        if (inpl3d(iplot,2) .eq. -1) then
        
        if (i1 .eq. i2 .and. i1 .eq. 1) then
          call cctogp(jdim,kdim,idim,i1+1,i2+1,i3,j1,j2,j3,k1,k2,k3,q,
     .              qi0,qj0,qk0,jdw,kdw,idw,xv,ldw)
        end if
        if (i1 .eq. i2 .and. i1 .eq. idim) then
          call cctogp(jdim,kdim,idim,i1-1,i2-1,i3,j1,j2,j3,k1,k2,k3,q,
     .              qi0,qj0,qk0,jdw,kdw,idw,xv,ldw)
        end if
        if (j1 .eq. j2 .and. j1 .eq. 1) then
          call cctogp(jdim,kdim,idim,i1,i2,i3,j1+1,j2+1,j3,k1,k2,k3,q,
     .              qi0,qj0,qk0,jdw,kdw,idw,xv,ldw)
        end if
        if (j1 .eq. j2 .and. j1 .eq. jdim) then
          call cctogp(jdim,kdim,idim,i1,i2,i3,j1-1,j2-1,j3,k1,k2,k3,q,
     .              qi0,qj0,qk0,jdw,kdw,idw,xv,ldw)
        end if
        if (k1 .eq. k2 .and. k1 .eq. 1) then
          call cctogp(jdim,kdim,idim,i1,i2,i3,j1,j2,j3,k1+1,k2+1,k3,q,
     .              qi0,qj0,qk0,jdw,kdw,idw,xv,ldw)
        end if
        if (k1 .eq. k2 .and. k1 .eq. kdim) then
          call cctogp(jdim,kdim,idim,i1,i2,i3,j1,j2,j3,k1-1,k2-1,k3,q,
     .              qi0,qj0,qk0,jdw,kdw,idw,xv,ldw)
        end if
        else
          call cctogp(jdim,kdim,idim,i1,i2,i3,j1,j2,j3,k1,k2,k3,q,
     .                qi0,qj0,qk0,jdw,kdw,idw,xv,ldw)
        end if
      end if
      
      if (ifunc .eq. -2) then
c
c       The size of xvist3d is (jdim,kdim,idim) in CFL3D 6.7.
c       It should be (jdw,kdw,idw). Otherwise, problematic in i2d=1 cases
c       because idim(2).ne.idw(1) 
c
        allocate(xvist3d(jdw,kdw,idw))
        call cctogp(jdim,kdim,idim,i1,i2,i3,j1,j2,j3,k1,k2,k3,vist3d,
     .              vi0,vj0,vk0,jdw,kdw,idw,xvist3d,1)
        allocate(xturre(jdw,kdw,idw,nummem))
        call cctogp(jdim,kdim,idim,i1,i2,i3,j1,j2,j3,k1,k2,k3,turre,
     .              ti0,tj0,tk0,jdw,kdw,idw,xturre,nummem)
        allocate(xvdsp(jdw,kdw,idw,nvdsp))
        call cctogp(jdim,kdim,idim,i1,i2,i3,j1,j2,j3,k1,k2,k3,vdsp,
     .              vdi0,vdj0,vdk0,jdw,kdw,idw,xvdsp,nvdsp)
      end if
c
      end if
c
      jw = (j2-j1)/j3 + 1
      kw = (k2-k1)/k3 + 1
      iw = (i2-i1)/i3 + 1
c
#     ifdef CMPLX
c
c     get derivative setp size on the local processor
c
      call gradinfo(delh,-1)
#     endif
c
      if (-2.eq.0) then
c       special: if inpl3d(iplot,2) is -1, find approx wall-parallel components
        if (inpl3d(iplot,2) .eq. -1) then
          if (i1 .eq. i2) then
            ix=0
            do i=i1,i2,i3
              ix=ix+1
              im1=max(1,i-1)
              kx=0
              do k=k1,k2,k3
                kx=kx+1
                km1=max(1,k-1)
                jx=0
                do j=j1,j2,j3
                  jx=jx+1
                  jm1=max(1,j-1)
                  sim1=(si(j  ,k  ,i  ,1)+si(jm1,k  ,i  ,1)+
     +                  si(j  ,km1,i  ,1)+si(jm1,km1,i  ,1))*.25
                  sim2=(si(j  ,k  ,i  ,2)+si(jm1,k  ,i  ,2)+
     +                  si(j  ,km1,i  ,2)+si(jm1,km1,i  ,2))*.25
                  sim3=(si(j  ,k  ,i  ,3)+si(jm1,k  ,i  ,3)+
     +                  si(j  ,km1,i  ,3)+si(jm1,km1,i  ,3))*.25
                  vnorm = xv(jx,kx,ix,2)*sim1 +
     +                    xv(jx,kx,ix,3)*sim2 +
     +                    xv(jx,kx,ix,4)*sim3
                  xv(jx,kx,ix,2) = xv(jx,kx,ix,2)-vnorm*sim1
                  xv(jx,kx,ix,3) = xv(jx,kx,ix,3)-vnorm*sim2
                  xv(jx,kx,ix,4) = xv(jx,kx,ix,4)-vnorm*sim3
                enddo
              enddo
            enddo
          end if
          if (j1 .eq. j2) then
            ix=0
            do i=i1,i2,i3
              ix=ix+1
              im1=max(1,i-1)
              kx=0
              do k=k1,k2,k3
                kx=kx+1
                km1=max(1,k-1)
                jx=0
                do j=j1,j2,j3
                  jx=jx+1
                  jm1=max(1,j-1)
                  sjm1=(sj(j  ,k  ,i  ,1)+sj(j  ,k  ,im1,1)+
     +                  sj(j  ,km1,i  ,1)+sj(j  ,km1,im1,1))*.25
                  sjm2=(sj(j  ,k  ,i  ,2)+sj(j  ,k  ,im1,2)+
     +                  sj(j  ,km1,i  ,2)+sj(j  ,km1,im1,2))*.25
                  sjm3=(sj(j  ,k  ,i  ,3)+sj(j  ,k  ,im1,3)+
     +                  sj(j  ,km1,i  ,3)+sj(j  ,km1,im1,3))*.25
                  vnorm = xv(jx,kx,ix,2)*sjm1 +
     +                    xv(jx,kx,ix,3)*sjm2 +
     +                    xv(jx,kx,ix,4)*sjm3
                  xv(jx,kx,ix,2) = xv(jx,kx,ix,2)-vnorm*sjm1
                  xv(jx,kx,ix,3) = xv(jx,kx,ix,3)-vnorm*sjm2
                  xv(jx,kx,ix,4) = xv(jx,kx,ix,4)-vnorm*sjm3
                enddo
              enddo
            enddo
          end if
          if (k1 .eq. k2) then
            ix=0
            do i=i1,i2,i3
              ix=ix+1
              im1=max(1,i-1)
              kx=0
              do k=k1,k2,k3
                kx=kx+1
                km1=max(1,k-1)
                jx=0
                do j=j1,j2,j3
                  jx=jx+1
                  jm1=max(1,j-1)
                  skm1=(sk(j  ,k  ,i  ,1)+sk(j  ,k  ,im1,1)+
     +                  sk(jm1,k  ,i  ,1)+sk(jm1,k  ,im1,1))*.25
                  skm2=(sk(j  ,k  ,i  ,2)+sk(j  ,k  ,im1,2)+
     +                  sk(jm1,k  ,i  ,2)+sk(jm1,k  ,im1,2))*.25
                  skm3=(sk(j  ,k  ,i  ,3)+sk(j  ,k  ,im1,3)+
     +                  sk(jm1,k  ,i  ,3)+sk(jm1,k  ,im1,3))*.25
                  vnorm = xv(jx,kx,ix,2)*skm1 +
     +                    xv(jx,kx,ix,3)*skm2 +
     +                    xv(jx,kx,ix,4)*skm3
                  xv(jx,kx,ix,2) = xv(jx,kx,ix,2)-vnorm*skm1
                  xv(jx,kx,ix,3) = xv(jx,kx,ix,3)-vnorm*skm2
                  xv(jx,kx,ix,4) = xv(jx,kx,ix,4)-vnorm*skm3
                enddo
              enddo
            enddo
          end if
        end if
c
c        load solution (q) in xv array (conserved variables)
c
         do 4000 i=1,iw
         do 4000 k=1,kw
         do 4000 j=1,jw
         xv(j,k,i,5) = xv(j,k,i,5)/gm1+0.5*(xv(j,k,i,2)**2
     .               + xv(j,k,i,3)**2+xv(j,k,i,4)**2)*xv(j,k,i,1)
         xv(j,k,i,2) = xv(j,k,i,1)*xv(j,k,i,2)
         xv(j,k,i,3) = xv(j,k,i,1)*xv(j,k,i,3)
         xv(j,k,i,4) = xv(j,k,i,1)*xv(j,k,i,4)
 4000    continue
c
      else
c
c        load desired function in xv array
c
         if (-2.eq.-5) then
            cpc = 2.e0/(gamma*xmach*xmach)
            do 4030 i=1,iw
            do 4030 k=1,kw
            do 4030 j=1,jw
            xv(j,k,i,1) = (xv(j,k,i,5)/p0-1)*cpc
 4030       continue
         end if
c
         if (-2.eq.-6) then
            do 4032 i=1,iw
            do 4032 k=1,kw
            do 4032 j=1,jw
            xv(j,k,i,1) = xv(j,k,i,5)/p0
 4032       continue
         end if
c
         if (-2.eq.-4) then
            do i=1,iw
              do k=1,kw
                do j=1,jw
                  xv(j,k,i,1) = xvist3d(j,k,i)
                end do
              end do
            end do
         end if
c
         if (-2.eq.-2) then
         iw = 0
         do i=i1,i2,i3
           iw = iw + 1
           jw = 0
           do j=j1,j2,j3
             jw = jw + 1
             kw = 0
             do k=k1,k2,k3
               kw = kw + 1
c
c               rho_avg,rho_rms
c
                xv(jw,kw,iw,1) = qavg(j,k,i,1)
                xv(jw,kw,iw,2) = sqrt(q2avg(j,k,i,1)
     .                         - qavg(j,k,i,1)*qavg(j,k,i,1))
c
c               u_avg_x,u_avg_y,u_avg_z
c
                xv(jw,kw,iw,3) = qavg(j,k,i,2)
                xv(jw,kw,iw,4) = qavg(j,k,i,3)
                xv(jw,kw,iw,5) = qavg(j,k,i,4)
c
c               u_rms_x,u_rms_y,u_rms_z
c
                xv(jw,kw,iw,6) = sqrt(q2avg(j,k,i,2)
     .                         - qavg(j,k,i,2)*qavg(j,k,i,2))
                xv(jw,kw,iw,7) = sqrt(q2avg(j,k,i,3)
     .                         - qavg(j,k,i,3)*qavg(j,k,i,3))
                xv(jw,kw,iw,8) = sqrt(q2avg(j,k,i,4)
     .                         - qavg(j,k,i,4)*qavg(j,k,i,4))
c
c               u_rey_x,u_rey_y,u_rey_z
c                    
                l=5+nummem
                xv(jw,kw,iw,9 ) = qs2avg(j,k,i,l)
     .                          - qavg(j,k,i,3)*qavg(j,k,i,4)
                l=6+nummem
                xv(jw,kw,iw,10) = qs2avg(j,k,i,l)
     .                          - qavg(j,k,i,2)*qavg(j,k,i,4)
                l=7+nummem
                xv(jw,kw,iw,11) = qs2avg(j,k,i,l)
     .                          - qavg(j,k,i,2)*qavg(j,k,i,3)
c
c               p,p_rms
c
                xv(jw,kw,iw,12) = qavg(j,k,i,5)
                xv(jw,kw,iw,13) = sqrt(q2avg(j,k,i,5)
     .                          - qavg(j,k,i,5)*qavg(j,k,i,5))
c
c               T,mul,mut,omega,k,Ma
c
                ixv = 13
                do l = 1,4+nummem
                ixv = ixv + 1
                xv(jw,kw,iw,ixv) = qsavg(j,k,i,l)
                ixv = ixv + 1
                xv(jw,kw,iw,ixv) = sqrt(qs2avg(j,k,i,l)
     .                           - qsavg(j,k,i,l)*qsavg(j,k,i,l))
                end do
c
c               vort_avg_x,vort_avg_y,vort_avg_z
c
                ixv = ixv + 1
                xv(jw,kw,iw,ixv) = vdavg(j,k,i,1)
                ixv = ixv + 1
                xv(jw,kw,iw,ixv) = vdavg(j,k,i,2)
                ixv = ixv + 1
                xv(jw,kw,iw,ixv) = vdavg(j,k,i,3)
c
c               vort_rms_x,vort_rms_y,vort_rms_z
c
                ixv = ixv + 1
                xv(jw,kw,iw,ixv) = sqrt(vd2avg(j,k,i,1)
     .                           - vdavg(j,k,i,1)*vdavg(j,k,i,1))
                ixv = ixv + 1
                xv(jw,kw,iw,ixv) = sqrt(vd2avg(j,k,i,2)
     .                           - vdavg(j,k,i,2)*vdavg(j,k,i,2))
                ixv = ixv + 1
                xv(jw,kw,iw,ixv) = sqrt(vd2avg(j,k,i,3)
     .                           - vdavg(j,k,i,3)*vdavg(j,k,i,3))
c
c               q0_avg,gradr,gradp...
c
                do l = 4,nvdsp
                ixv = ixv + 1
                xv(jw,kw,iw,ixv) = vdavg(j,k,i,l)
                ixv = ixv + 1
                xv(jw,kw,iw,ixv) = sqrt(vd2avg(j,k,i,l)
     .                           - vdavg(j,k,i,l)*vdavg(j,k,i,l))
                end do
c
             end do
           end do
         end do
         end if
c
      end if
c
c     load xv into single precision xw array
c     (for complex code, load either real or complex part
c     of xv into single precision xw array, depending
c     on whether solution or derivative is to be output)
c
      leq = 5
      if (ifunc .ne. 0 ) leq = 1
      if (-2 .eq. -2) leq = 2*ifuncdim + 3
c
      jw = (j2-j1)/j3 + 1
      kw = (k2-k1)/k3 + 1
      iw = (i2-i1)/i3 + 1
#     ifdef CMPLX
c
      if (ip3dgrad .eq. 1) then
c        output derivative of solution
         do l=1,leq
            do i=1,iw
               do k=1,kw
                  do j=1,jw
                     xw(j,k,i,l) = imag(xv(j,k,i,l))/real(delh)
                  end do
               end do
            end do
         end do
      else
c        output solution
         do l=1,leq
            do i=1,iw
               do k=1,kw
                  do j=1,jw
                     xw(j,k,i,l) = real(xv(j,k,i,l))
                  end do
               end do
            end do
         end do
      end if
#     else
      do l=1,leq
         do i=1,iw
            do k=1,kw
               do j=1,jw
                  xw(j,k,i,l) = xv(j,k,i,l)
               end do
            end do
         end do
      end do
#     endif
c
#if defined DIST_MPI
      end if
c
c     send/receive plot3d q data
c
      jw = (j2-j1)/j3 + 1
      kw = (k2-k1)/k3 + 1
      iw = (i2-i1)/i3 + 1
c
      if (-2.eq.0) then
         nnq = jw*kw*iw*5
      else if (-2.eq.-2) then
         nnq = jw*kw*iw*(2*ifuncdim+3)
      else
         nnq = jw*kw*iw*1
      end if
      nd_srce = mblk2nd(nblk)
      mytag = itag_q + nblk
c
      if (mblk2nd(nblk).eq.myid) then
         call MPI_Send (xw, nnq, MY_MPI_REAL, myhost, mytag,
     &                  mycomm, ierr)
      else if (myid .eq. myhost) then
         call MPI_Recv (xw, nnq, MY_MPI_REAL, nd_srce,
     &                  mytag, mycomm, istat, ierr)
      end if
#endif
c
      if (myid .eq. myhost) then
c
c     output solution
c
      if (-2 .eq. 0) then
         if (ibin.eq.0) then
            write(97,'(5e14.6)') xmachw,alphww,reuew,timew
            if (i2d.eq.0) then
               if (ialph.eq.0) then
                  write(97,'(5e14.6)')
     .            ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=1,5)
               else
                  do i=1,iw
                     do j=1,jw
                        do k=1,kw
                           xw(j,k,i,3) = -xw(j,k,i,3)
                        end do
                     end do
                  end do
                  write(97,'(5e14.6)')
     .            ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=1,2),
     .             (((xw(j,k,i,4),i=1,iw),j=1,jw),k=1,kw),
     .             (((xw(j,k,i,3),i=1,iw),j=1,jw),k=1,kw),
     .             (((xw(j,k,i,5),i=1,iw),j=1,jw),k=1,kw) 
               end if
            else
               write(97,'(5e14.6)')
     .         ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=1,2),
     .         ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=4,5)
            end if
         else
            write(97)  xmachw,alphww,reuew,timew
            if (i2d.eq.0) then
               if (ialph.eq.0) then
                  write(97)
     .            ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=1,5)
               else
                  do i=1,iw
                     do j=1,jw
                        do k=1,kw
                           xw(j,k,i,3) = -xw(j,k,i,3)
                        end do
                     end do
                  end do
                  write(97)
     .            ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=1,2),
     .             (((xw(j,k,i,4),i=1,iw),j=1,jw),k=1,kw),
     .             (((xw(j,k,i,3),i=1,iw),j=1,jw),k=1,kw),
     .             (((xw(j,k,i,5),i=1,iw),j=1,jw),k=1,kw)
               end if
            else
               write(97)
     .         ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=1,2),
     .         ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=4,5)
            end if
         end if
      else if (-2 .eq. -2) then
c     function file at grid points
c     2D/3D share the same flow variables
         lt2  = 2*nummem   
         leq2 = 2*ifuncdim+3
         if (ibin.eq.0) then
            if (ialph.eq.0) then
              write(97,'(5e14.6)')
     .         ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=1,5)
            else
              do i=1,iw
                do j=1,jw
                  do k=1,kw
                      xw(j,k,i,4)      = -xw(j,k,i,4)      !u_avg_y in solver
                      xw(j,k,i,9)      = -xw(j,k,i,9)      !u_rey_x in solver
                      xw(j,k,i,11)     = -xw(j,k,i,11)     !u_rey_z in solver
                      xw(j,k,i,23+lt2) = -xw(j,k,i,23+lt2) !vort_avg_y in solver
                  end do
                end do
              end do
              write(97,'(5e14.6)')
c             rho
     .        ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=1,2),
c             u
     .        ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=3,3),
     .        ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=5,5),
     .        ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=4,4),
     .        ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=6,6),
     .        ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=8,8),
     .        ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=7,7),
c             u_rey
     .        ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=9,9),
     .        ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=11,11),
     .        ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=10,10),
c             p,T,mul,mut,omega,k,Ma
     .        ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=12,21+lt2),
c             vort
     .        ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=22+lt2,22+lt2),
     .        ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=24+lt2,24+lt2),
     .        ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=23+lt2,23+lt2),
     .        ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=25+lt2,25+lt2),
     .        ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=27+lt2,27+lt2),
     .        ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=26+lt2,26+lt2),
c             q0,gradr,gradp...
     .        ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=28+lt2,leq2)
            end if
         else
            if (ialph.eq.0) then
              write(97)
     .         ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=1,5)
            else
              do i=1,iw
                do j=1,jw
                  do k=1,kw
                      xw(j,k,i,4)      = -xw(j,k,i,4)      !u_avg_y
                      xw(j,k,i,9)      = -xw(j,k,i,9)      !u_rey_x
                      xw(j,k,i,11)     = -xw(j,k,i,11)     !u_rey_z
                      xw(j,k,i,23+lt2) = -xw(j,k,i,23+lt2) !vort_avg_y
                  end do
                end do
              end do
              write(97)
c             rho
     .        ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=1,2),
c             u
     .        ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=3,3),
     .        ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=5,5),
     .        ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=4,4),
     .        ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=6,6),
     .        ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=8,8),
     .        ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=7,7),
c             u_rey
     .        ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=9,9),
     .        ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=11,11),
     .        ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=10,10),
c             p,T,mul,mut,omega,k,Ma
     .        ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=12,21+lt2),  
c             vort
     .        ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=22+lt2,22+lt2),
     .        ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=24+lt2,24+lt2),
     .        ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=23+lt2,23+lt2),
     .        ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=25+lt2,25+lt2),
     .        ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=27+lt2,27+lt2),
     .        ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=26+lt2,26+lt2),
c             q0,gradr,gradp...
     .        ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=28+lt2,leq2)
            end if
         end if
      else
         if (ibin.eq.0) then
            write(97,'(5e14.6)')
     .      (((xw(j,k,i,1),i=1,iw),j=1,jw),k=1,kw)
         else
            write(97)
     .      (((xw(j,k,i,1),i=1,iw),j=1,jw),k=1,kw)
         end if
      end if
c
      end if
c
      end if
c
      if (iflag.eq.2) then
c
c***********************************************************************
c     print out data
c***********************************************************************
c
c     NOTE: for the first section of printed output (up to the end of
c           the 5000 loop), turb. viscosity data are stored in the 4th
c           position of the xg array as follows
c
c           xg(j,k,i,4) = vist3d
c
c           the solution q is stored in the xw array as follows
c
c           xw(j,k,i,1) = q(1) (rho)
c           xw(j,k,i,2) = q(2) (u)
c           xw(j,k,i,3) = q(3) (v)
c           xw(j,k,i,4) = q(4) (w)
c           xw(j,k,i,5) = q(5) (p)
c
c           for the subsequent sections of printed output (loops 6002,
c           7002 and 8002), the following storage is used:
c
c           xg(j,k,i,4) = dn (normal spacing away from a surface)
c
c           xw(j,k,i,1) = pres1 (pressure)
c           xw(j,k,i,2) = temp1 (temperature)
c           xw(j,k,i,3) = Cf    (skin friction coeff.)
c           xw(j,k,i,4) = Ch    (heat transfer coeff.)
c           xw(j,k,i,5) = yplus (turbulent length scale y+)
c
c           in all sections xg(1-3) contain x,y,z, as set at the top
c           of the subroutine
c
      if (lhdr.gt.0 .and. myid.eq.myhost) then
         write(11,94)idim,jdim,kdim
   94    format(47h writing printout file for IDIM X JDIM X KDIM =,
     .   i5,3h x ,i5,3h x ,i5,5h grid)
      end if
c
#if defined DIST_MPI
      if(mblk2nd(nblk).eq.myid) then
#endif
      if(j1.eq.1 .and. j2.eq.1) then
c
        do l = 1,ifuncdim
          iw = 0
          do i=i1,i2,i3
          iw = iw + 1
          kw = 0
          do k=k1,k2,k3
          kw = kw + 1
          jw = 0
          do j=j1,j2,j3
          jw = jw + 1
            xv(jw,kw,iw,l*2-1) = vsj0(k,i,l,1)
            xv(jw,kw,iw,l*2)   = vsj0(k,i,l,2) 
     .                         - vsj0(k,i,l,1)*vsj0(k,i,l,1)
          end do
          end do
          end do
        end do
c
      else if(j1.eq.jdim .and. j2.eq.jdim) then
c
        do l = 1,ifuncdim
          iw = 0
          do i=i1,i2,i3
          iw = iw + 1
          kw = 0
          do k=k1,k2,k3
          kw = kw + 1
          jw = 0
          do j=j1,j2,j3
          jw = jw + 1
            xv(jw,kw,iw,l*2-1) = vsjdim(k,i,l,1)
            xv(jw,kw,iw,l*2)   = vsjdim(k,i,l,2) 
     .                         - vsjdim(k,i,l,1)*vsjdim(k,i,l,1)
          end do
          end do
          end do
        end do
c
      else if(k1.eq.1 .and. k2.eq.1) then
c
        do l = 1,ifuncdim
          iw = 0
          do i=i1,i2,i3
          iw = iw + 1
          kw = 0
          do k=k1,k2,k3
          kw = kw + 1
          jw = 0
          do j=j1,j2,j3
          jw = jw + 1
            xv(jw,kw,iw,l*2-1) = vsk0(j,i,l,1)
            xv(jw,kw,iw,l*2)   = vsk0(j,i,l,2)
     .                         - vsk0(j,i,l,1)*vsk0(j,i,l,1)
          end do
          end do
          end do
        end do
c
      else if(k1.eq.kdim .and. k2.eq.kdim) then
c
        do l = 1,ifuncdim
          iw = 0
          do i=i1,i2,i3
          iw = iw + 1
          kw = 0
          do k=k1,k2,k3
          kw = kw + 1
          jw = 0
          do j=j1,j2,j3
          jw = jw + 1
            xv(jw,kw,iw,l*2-1) = vskdim(j,i,l,1)
            xv(jw,kw,iw,l*2)   = vskdim(j,i,l,2)
     .                         - vskdim(j,i,l,1)*vskdim(j,i,l,1)
          end do
          end do
          end do
        end do
c    
      else if(i1.eq.1 .and. i2.eq.1) then
c
        do l = 1,ifuncdim
          iw = 0
          do i=i1,i2,i3
          iw = iw + 1
          kw = 0
          do k=k1,k2,k3
          kw = kw + 1
          jw = 0
          do j=j1,j2,j3
          jw = jw + 1
            xv(jw,kw,iw,l*2-1) = vsi0(j,k,l,1)
            xv(jw,kw,iw,l*2)   = vsi0(j,k,l,2)
     .                         - vsi0(j,k,l,1)*vsi0(j,k,l,1)
          end do
          end do
          end do
        end do
c
      else if(i1.eq.idim .and. i2.eq.idim) then
c
        do l = 1,ifuncdim
          iw = 0
          do i=i1,i2,i3
          iw = iw + 1
          kw = 0
          do k=k1,k2,k3
          kw = kw + 1
          jw = 0
          do j=j1,j2,j3
          jw = jw + 1
            xv(jw,kw,iw,l*2-1) = vsidim(j,k,l,1)
            xv(jw,kw,iw,l*2)   = vsidim(j,k,l,2)
     .                         - vsidim(j,k,l,1)*vsidim(j,k,l,1)
          end do
          end do
          end do
        end do
c
      end if
c
      ldw = 2*ifuncdim+3
c     set as 2*ifuncdim+3 for consistency with average flowfield
c     it is enough to set as 2*ifuncdim
      do l=1,ldw
         do i=1,idw
            do k=1,kdw
               do j=1,jdw
                  xw(j,k,i,l) = xv(j,k,i,l)
               end do
            end do
         end do
      end do
c
#if defined DIST_MPI
      end if
#endif
c
#if defined DIST_MPI
c
c     send/receive x,y,z,q and vist3d data
c
      jw = (j2-j1)/j3 + 1
      kw = (k2-k1)/k3 + 1
      iw = (i2-i1)/i3 + 1
c
      nnq = jw*kw*iw*(2*ifuncdim+3)
      nd_srce = mblk2nd(nblk)
      mytag = itag_q + nblk
c
      if (mblk2nd(nblk).eq.myid) then
         call MPI_Send (xw, nnq, MY_MPI_REAL, myhost, mytag,
     &                  mycomm, ierr)
      else if (myid.eq.myhost) then
         call MPI_Recv (xw, nnq, MY_MPI_REAL, nd_srce,
     &                  mytag, mycomm, istat, ierr)
     
      end if
c
      nng = jw*kw*iw*4
      nd_srce = mblk2nd(nblk)
      mytag = itag_grid + nblk
c
      if (mblk2nd(nblk).eq.myid) then
         call MPI_Send (xg, nng, MY_MPI_REAL, myhost, mytag,
     &                  mycomm, ierr)
      else if (myid.eq.myhost) then
         call MPI_Recv (xg, nng, MY_MPI_REAL, nd_srce,
     &                  mytag, mycomm, istat, ierr)
      end if
c
#endif
c
      if (myid .eq. myhost) then
         if (ibin.eq.0) then
            if (ialph.eq.0) then
c             to be considered
            else
              if(icallavg.eq.0) then
                if(i2d.eq.0) then
                  do i=1,iw
                    do k=1,kw
                      do j=1,jw
                        xg(j,k,i,2) = -xg(j,k,i,2)
                      end do
                    end do
                  end do
                  write(176,'(5e17.9)')
c                 x,y,z
     .            ((((xg(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=1,1),
     .            ((((xg(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=3,3),
     .            ((((xg(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=2,2)
                else
                   write(176,'(5e17.9)')
c                  x,y
     .             ((((xg(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=1,1),
     .             ((((xg(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=3,3)
                end if
              end if
c
              do i=1,iw
                do k=1,kw
                  do j=1,jw
                      xw(j,k,i,3) = -xw(j,k,i,3)
                      xw(j,k,i,25)= -xw(j,k,i,25)
                  end do
                end do
              end do
              write(177,'(5e17.9)')            
c              u (off-body by default)
     .         ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=1,1),
     .         ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=5,5),
     .         ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=3,3),
     .         ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=2,2),
     .         ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=6,6),
     .         ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=4,4),
c              p,T,M,Cp,mut,dn,Ch,yplus
     .         ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=7,22),
c              Cf
     .         ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=23,23),
     .         ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=27,27),  
     .         ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=25,25),
     .         ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=24,24),
     .         ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=28,28),
     .         ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=26,26)
            end if
         else
            if (ialph.eq.0) then
c             to be considered
            else
              if(icallavg.eq.0) then
                if(i2d.eq.0) then
                  do i=1,iw
                    do k=1,kw
                      do j=1,jw
                        xg(j,k,i,2) = -xg(j,k,i,2)
                      end do
                    end do
                  end do
                  write(176)
c                 x,y,z
     .            ((((xg(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=1,1),
     .            ((((xg(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=3,3),
     .            ((((xg(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=2,2)
                else
                  write(176)
c                 x,y
     .            ((((xg(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=1,1),
     .            ((((xg(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=3,3)
                end if
              end if
c
              do i=1,iw
                do k=1,kw
                  do j=1,jw
                      xw(j,k,i,3) = -xw(j,k,i,3)
                      xw(j,k,i,25)= -xw(j,k,i,25)
                  end do
                end do
              end do
              write(177)
c              u (off-body by default)
     .         ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=1,1),
     .         ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=5,5),
     .         ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=3,3),
     .         ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=2,2),
     .         ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=6,6),
     .         ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=4,4),
c              p,T,M,Cp,mut,dn,Ch,yplus
     .         ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=7,22),
c              Cf
     .         ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=23,23),
     .         ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=27,27),  
     .         ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=25,25),
     .         ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=24,24),
     .         ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=28,28),
     .         ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=26,26)
            end if
         end if
      end if
c
      end if
      return
      end

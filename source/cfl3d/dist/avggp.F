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
      subroutine avggp(jdim,kdim,idim,                            
     .                 i1,i2,i3,j1,j2,j3,k1,k2,k3,q,                 
     .                 qi0,qj0,qk0,x,y,z,xw,blank2,blank,xg,iflag,
     .                 vist3d,iover,nblk,nmap,bcj,bck,bci,
     .                 vj0,vk0,vi0,ifunc,iplot,jdw,kdw,idw,
     .                 nplots,jdimg,kdimg,idimg,nblcg,jsg,ksg,isg,
     .                 jeg,keg,ieg,ninter,iindex,intmax,nsub1,
     .                 maxxe,nblkk,nbli,limblk,isva,nblon,mxbli,
     .                 thetay,maxbl,maxgr,myid,myhost,mycomm,
     .                 mblk2nd,inpl3d,nblock,nblkpt,xv,
     .                 sj,sk,si,vol,nummem,ifuncdim,
     .                 turre,tj0,tk0,ti0,vdsp,nvdsp,vdj0,vdk0,vdi0,
     .                 qavg,q2avg,qsavg,qs2avg,vdavg,vd2avg,
     .                 vsj0,vsjdim,vsk0,vskdim,vsi0,vsidim)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Write the output file at the grid points in PLOT3D
c     format and print solution data.
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
#   ifdef P3D_SINGLE
#     define MY_MPI_REAL MPI_REAL
#   else
#     define MY_MPI_REAL MPI_DOUBLE_PRECISION
#   endif
#endif
#if defined P3D_SINGLE
      real*4    xw(jdw,kdw,idw,ifuncdim),xg(jdw,kdw,idw,4)
      real*4    xmachw,alphww,reuew,timew
#else
      real xw(jdw,kdw,idw,ifuncdim),xg(jdw,kdw,idw,4)
#endif

      dimension xv(jdw,kdw,idw,ifuncdim)
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
      dimension turre(jdim,kdim,idim,nummem),tj0(kdim,idim-1,nummem,4),
     .          tk0(jdim,idim-1,nummem,4),ti0(jdim,kdim,nummem,4)
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
      common /unclass/ ioffbody,isample
c
c     initialize xw, xv and xg arrays
c
      idim1 = idim - 1
      jdim1 = jdim - 1
      kdim1 = kdim - 1
      jdw = (j2-j1)/j3 + 1
      kdw = (k2-k1)/k3 + 1
      idw = (i2-i1)/i3 + 1
      do iw=1,idw
         do kw=1,kdw
            do jw=1,jdw
               do l= 1,5
               xw(jw,kw,iw,l) = 0.
               xv(jw,kw,iw,l) = 0.
               end do
               do l= 1,4
               xg(jw,kw,iw,l) = 0.
               end do
            end do
         end do
      end do
c
c     average field solutions at grid points
c       
      if (iflag.eq.1) then
c       determine primitive values at grid points and load into xv array
        ldw = 5
        call cctogp(jdim,kdim,idim,i1,i2,i3,j1,j2,j3,k1,k2,k3,q,
     .              qi0,qj0,qk0,jdw,kdw,idw,xv,ldw)
c       eddy viscosity
        allocate(xvist3d(jdw,kdw,idw))
        call cctogp(jdim,kdim,idim,i1,i2,i3,j1,j2,j3,k1,k2,k3,vist3d,
     .                vi0,vj0,vk0,jdw,kdw,idw,xvist3d,1)
c       turbulent variables
        allocate(xturre(jdw,kdw,idw,nummem))
        call cctogp(jdim,kdim,idim,i1,i2,i3,j1,j2,j3,k1,k2,k3,turre,
     .                ti0,tj0,tk0,jdw,kdw,idw,xturre,nummem)
c       vdsp variables
        allocate(xvdsp(jdw,kdw,idw,nvdsp))
        call cctogp(jdim,kdim,idim,i1,i2,i3,j1,j2,j3,k1,k2,k3,vdsp,
     .                vdi0,vdj0,vdk0,jdw,kdw,idw,xvdsp,nvdsp)
c
        do iw=1,idw
        do kw=1,kdw
        do jw=1,jdw
          q1 = xv(jw,kw,iw,1)
          q2 = xv(jw,kw,iw,2)
          q3 = xv(jw,kw,iw,3)
          q4 = xv(jw,kw,iw,4)
          q5 = xv(jw,kw,iw,5)
c         temperature (ixv=6)
          ixv = 6
          tt = gamma*q5/q1
          xv(jw,kw,iw,ixv) = tt
c         laminar viscosity (7)
          ixv = ixv + 1
          c2b = cbar/tinf
          c2bp= c2b+1.0
          fnu = c2bp*tt*sqrt(tt)/(c2b+tt)
          xv(jw,kw,iw,ixv) = fnu
c         turbulent viscosity (8)
          ixv = ixv + 1
          xv(jw,kw,iw,ixv) = xvist3d(jw,kw,iw)
c         turbulent quantities (9,10 if two-equation model)
          do l = 1,nummem
          ixv = ixv + 1
          xv(jw,kw,iw,ixv) = xturre(jw,kw,iw,l)
          end do
c         mach number (11)
          ixv = ixv + 1
          cc  = gamma*q5/q1
          xma = sqrt((q2*q2+q3*q3+q4*q4)/cc)
          xv(jw,kw,iw,ixv) = xma
        end do
        end do
        end do
c       vortx,vorty,vortz,q0,gradr,gradp,blank... (12-17,18,19,20,21)
c       xvdsp(jw,kw,iw,1:nvdsp)
c       
c       qavg/q2avg
c       
        do l = 1,5
          do iw=1,idw
          do kw=1,kdw
          do jw=1,jdw
            qavg(jw,kw,iw,l)=(qavg(jw,kw,iw,l)*(xnumavg-1.)+
     +        xv(jw,kw,iw,l))/xnumavg
            q2avg(jw,kw,iw,l)=(q2avg(jw,kw,iw,l)*(xnumavg-1.)+          
     +        xv(jw,kw,iw,l)*xv(jw,kw,iw,l))/xnumavg
          end do
          end do
          end do
        end do
c       
c       qavg/q2avg augmentation
c       
        do l = 1,4+nummem
          ixv = l + 5
          do iw=1,idw
          do kw=1,kdw
          do jw=1,jdw  
            qsavg (jw,kw,iw,l)=(qsavg (jw,kw,iw,l)*(xnumavg-1.)+
     +        xv(jw,kw,iw,ixv))/xnumavg
            qs2avg(jw,kw,iw,l)=(qs2avg(jw,kw,iw,l)*(xnumavg-1.)+
     +        xv(jw,kw,iw,ixv)*xv(jw,kw,iw,ixv))/xnumavg
          end do
          end do
          end do
        end do
c
c       u_rey_x,u_rey_y,u_rey_z
c       qs2avg has the same size as qsavg for consistency
c
        l = 5+nummem
        do iw=1,idw
        do kw=1,kdw
        do jw=1,jdw  
            !qsavg (jw,kw,iw,l)=0.0
            qs2avg(jw,kw,iw,l)=(qs2avg(jw,kw,iw,l)*(xnumavg-1.)
     +       +xv(jw,kw,iw,3)*xv(jw,kw,iw,4))/xnumavg
        end do
        end do
        end do
        l = 6+nummem
        do iw=1,idw
        do kw=1,kdw
        do jw=1,jdw  
            !qsavg (jw,kw,iw,l)=0.0
            qs2avg(jw,kw,iw,l)=(qs2avg(jw,kw,iw,l)*(xnumavg-1.)
     +       +xv(jw,kw,iw,2)*xv(jw,kw,iw,4))/xnumavg
        end do
        end do
        end do 
        l = 7+nummem
        do iw=1,idw
        do kw=1,kdw
        do jw=1,jdw  
            !qsavg (jw,kw,iw,l)=0.0
            qs2avg(jw,kw,iw,l)=(qs2avg(jw,kw,iw,l)*(xnumavg-1.)
     +       +xv(jw,kw,iw,2)*xv(jw,kw,iw,3))/xnumavg
        end do
        end do
        end do
c
c       vdsp array
c       
        do l = 1,nvdsp
          do iw=1,idw
          do kw=1,kdw
          do jw=1,jdw  
              vdavg (jw,kw,iw,l)=(vdavg (jw,kw,iw,l)*(xnumavg-1.)+      
     +          xvdsp(jw,kw,iw,l))/xnumavg
              vd2avg(jw,kw,iw,l)=(vd2avg(jw,kw,iw,l)*(xnumavg-1.)+
     +          xvdsp(jw,kw,iw,l)*vdsp(jw,kw,iw,l))/xnumavg
          end do
          end do
          end do
        end do
c
      end if
c
c     average surface solutions at grid points
c
      if (iflag.eq.2) then
c
c       determine values at grid points and load into xv
c
        ldw = 5
        call cctogp(jdim,kdim,idim,i1,i2,i3,j1,j2,j3,k1,k2,k3,q,
     .              qi0,qj0,qk0,jdw,kdw,idw,xv,ldw)
c
c       eddy viscosity at grid points and load
c       into xvist3d
c       
        if (ivisc(3).gt.1 .or. ivisc(2).gt.1 .or. ivisc(1).gt.1) then
           allocate(xvist3d(jdw,kdw,idw))
           call cctogp(jdim,kdim,idim,i1,i2,i3,j1,j2,j3,k1,k2,k3,vist3d,
     .                 vi0,vj0,vk0,jdw,kdw,idw,xvist3d,1)
        end if
c      
c       inviscid-flow data
c
        term3 = 1./( (1.+0.5*gm1*xmach*xmach)**(gamma/gm1) )
        iw = 0
        do i=i1,i2,i3
        iw = iw + 1
        kw = 0
        do k=k1,k2,k3
        kw = kw + 1
        jw = 0
        do j=j1,j2,j3
        jw = jw + 1
c
          q1    = xv(jw,kw,iw,1)
          q2    = xv(jw,kw,iw,2)/xmach
          q3    = xv(jw,kw,iw,3)/xmach
          q4    = xv(jw,kw,iw,4)/xmach
          q5    = gamma*xv(jw,kw,iw,5)
          t1    = q5/q1
          xm1   = sqrt(xmach**2*(q2**2+q3**2+q4**2)/t1)
c         cp or pitot pressure
          icp   = 1
          if(icp .eq. 1) then
             pitot = 2.*(q5-1.)/(gamma*xmach*xmach)
          else
             if (real(xm1).gt.1.0) then
                term1 = (0.5*gp1*xm1*xm1)**(gamma/gm1)
                term2 = (2.*gamma*xm1*xm1/gp1 - gm1/gp1)**(1./gm1)
                pitot = q5*term1*term3/term2
             else
                term1 = (1.0+0.5*gm1*xm1*xm1)**(gamma/gm1)
                pitot = q5*term1*term3
             end if
          end if
c         turbulent viscosity
          edvis = xvist3d(jw,kw,iw)
c         
          xv(jw,kw,iw,1) = real(q2)
          xv(jw,kw,iw,2) = real(q3)
          xv(jw,kw,iw,3) = real(q4)
          xv(jw,kw,iw,4) = real(q5)
          xv(jw,kw,iw,5) = real(t1)
          xv(jw,kw,iw,6) = real(xm1)
          xv(jw,kw,iw,7) = real(pitot)
          xv(jw,kw,iw,8) = real(edvis)
c
        end do
        end do
        end do
c
c       output viscous-flow data (skin friction, heat transfer, etc) on
c       j=1 and/or j=jdim surfaces
c
        if (ivisc(2).ne.0 .and. j1.eq.j2) then
c       
        iw = 0
        do 7002 i=i1,i2,i3
        iw = iw + 1
        id = min(i,idim1)
        id1= max(i-1,1)
        jw = 0
        do 7001 j=j1,j2,j3
        jw = jw + 1
        kw = 0
        do 7000 k=k1,k2,k3
        kw = kw + 1
        kd = min(k,kdim1)
        kd1= max(k-1,1)
        if (j.eq.1) then
           jd  = 1
           jdx = 1
           m   = 2
        else
           jd  = jdim1
           jdx = jdim
           m   = 4
        end if
c       
c       wall value - average in K and I (after QFACE, m=2 & 4 are interface
c       (wall) values
c       
        q1j1 =       .25*(qj0(kd,id,1,m)  +qj0(kd1,id,1,m)
     .                  + qj0(kd,id1,1,m) +qj0(kd1,id1,1,m))
        q2j1 =       .25*(qj0(kd,id,2,m)  +qj0(kd1,id,2,m)
     .                  + qj0(kd,id1,2,m) +qj0(kd1,id1,2,m))
        q3j1 =       .25*(qj0(kd,id,3,m)  +qj0(kd1,id,3,m)
     .                  + qj0(kd,id1,3,m) +qj0(kd1,id1,3,m))
        q4j1 =       .25*(qj0(kd,id,4,m)  +qj0(kd1,id,4,m)
     .                  + qj0(kd,id1,4,m) +qj0(kd1,id1,4,m))
        q5j1 = gamma*.25*(qj0(kd,id,5,m)  +qj0(kd1,id,5,m)
     .                  + qj0(kd,id1,5,m) +qj0(kd1,id1,5,m))
        t1j1 = q5j1/q1j1
c       
c       first cell center location - average in K and I
c       
        q1j2 =       .25*(q(jd,kd,id,1)   +q(jd,kd1,id,1)
     .                  + q(jd,kd,id1,1)  +q(jd,kd1,id1,1))
        q2j2 =       .25*(q(jd,kd,id,2)   +q(jd,kd1,id,2)
     .                  + q(jd,kd,id1,2)  +q(jd,kd1,id1,2))
        q3j2 =       .25*(q(jd,kd,id,3)   +q(jd,kd1,id,3)
     .                  + q(jd,kd,id1,3)  +q(jd,kd1,id1,3))
        q4j2 =       .25*(q(jd,kd,id,4)   +q(jd,kd1,id,4)
     .                  + q(jd,kd,id1,4)  +q(jd,kd1,id1,4))
        q5j2 = gamma*.25*(q(jd,kd,id,5)   +q(jd,kd1,id,5)
     .               + q(jd,kd,id1,5)  +q(jd,kd1,id1,5))
        t1j2 = q5j2/q1j2
c       
        if (sj(jdx,kd ,id ,4) .eq. 0. .or. 
     +      sj(jdx,kd1,id ,4) .eq. 0. .or.
     +      sj(jdx,kd ,id1,4) .eq. 0. .or.
     +      sj(jdx,kd1,id1,4) .eq. 0.) then
        dx = x(jd+1,kd,i)-x(jd,kd,i)
        dy = y(jd+1,kd,i)-y(jd,kd,i)
        dz = z(jd+1,kd,i)-z(jd,kd,i)
        dn = sqrt(dx**2+dy**2+dz**2)
        else
        dn = (vol(jd,kd ,id )/sj(jdx,kd ,id ,4)+
     +        vol(jd,kd1,id )/sj(jdx,kd1,id ,4)+
     +        vol(jd,kd ,id1)/sj(jdx,kd ,id1,4)+
     +        vol(jd,kd1,id1)/sj(jdx,kd1,id1,4))/4.
        end if
c       
c       Get turb viscosity at wall (0 unless wall fn used)
        if (iwf(2) .eq. 1) then
        avgmut =       .25*(vj0(kd,id,1,m)  +vj0(kd1,id,1,m)
     .                 + vj0(kd,id1,1,m) +vj0(kd1,id1,1,m))
        else
        avgmut = 0.
        end if
        emuka = (t1j1**1.5)*((1.0+cbar/tinf)/(t1j1+cbar/tinf))
c       
c       Use component of velocity parallel to wall
        urel = q2j2-q2j1
        vrel = q3j2-q3j1
        wrel = q4j2-q4j1
        sj1  = (sj(jdx,kd ,id ,1)+sj(jdx,kd1,id ,1)+
     +          sj(jdx,kd ,id1,1)+sj(jdx,kd1,id1,1))/4.
        sj2  = (sj(jdx,kd ,id ,2)+sj(jdx,kd1,id ,2)+
     +          sj(jdx,kd ,id1,2)+sj(jdx,kd1,id1,2))/4.
        sj3  = (sj(jdx,kd ,id ,3)+sj(jdx,kd1,id ,3)+
     +          sj(jdx,kd ,id1,3)+sj(jdx,kd1,id1,3))/4.
        vnorm = urel*sj1 + vrel*sj2 + wrel*sj3
        upar = urel-vnorm*sj1
        vpar = vrel-vnorm*sj2
        wpar = wrel-vnorm*sj3
        Cf   = 2.*sqrt(upar**2 + vpar**2 + wpar**2)
        Cf   = 2.0*(emuka+avgmut)/(reue*xmach)*Cf/dn
        if (real(q2j2).lt.0.) cf = -cf
c       Cf vector
        cfx  = 4.*(emuka+avgmut)/(reue*xmach)*upar/dn
        cfy  = 4.*(emuka+avgmut)/(reue*xmach)*vpar/dn
        cfz  = 4.*(emuka+avgmut)/(reue*xmach)*wpar/dn
c       
c       Note:  Ch definition not standard
c       
        eps  = 1.0e-6
        tty1 = 1.0 +gm1*0.5*xmach*xmach
        t1m1 = t1j1-tty1
        if (abs(real(t1m1)).le.real(eps)) then
           Ch = 999.99999
        else
c          2-point formula
              Ch = 2.*(t1j2-t1j1)
              Ch = (emuka+avgmut)/(reue*pr*(t1j1-tty1))*Ch/dn
        end if
c       
        pres1 = q5j1
        temp1 = t1j1
        yplus = 0.
c       yplus is value at first GRIDPOINT away from
        if (abs(real(Cf)).gt.0.) then
           yplus = dn*reue*q1j1*sqrt(ccabs(Cf*0.5/q1j1))/emuka
        end if
c       
c       store data in xv array
c       
        xv(jw,kw,iw,9)  = dn
        xv(jw,kw,iw,10) = Ch
        xv(jw,kw,iw,11) = yplus
        xv(jw,kw,iw,12) = cfx
        xv(jw,kw,iw,13) = cfy
        xv(jw,kw,iw,14) = cfz
c       
 7000   continue
 7001   continue
 7002   continue
c       
        end if
c       
c       viscous-flow data (skin friction, heat transfer, etc) on
c       k=1 and/or k=kdim surfaces
c       
        if (ivisc(3).ne.0 .and. k1.eq.k2) then
c       
        iw = 0
        do 6002 i=i1,i2,i3
        iw = iw + 1
        id = min(i,idim1)
        id1= max(i-1,1)
        kw = 0
        do 6001 k=k1,k2,k3
        kw = kw + 1
        jw = 0
        do 6000 j=j1,j2,j3
        jw = jw + 1
        jd = min(j,jdim1)
        jd1= max(j-1,1)
        if (k.eq.1) then
           kd  = 1
           kdx = 1
           m   = 2
        else
           kd  = kdim1
           kdx = kdim
           m   = 4
        end if
c       
c       wall value - average in J and I (after QFACE, m=2 & 4 are interface
c       (wall) values
        q1k1 =       .25*(qk0(jd,id,1,m)  +qk0(jd1,id,1,m)
     .                  + qk0(jd,id1,1,m) +qk0(jd1,id1,1,m))
        q2k1 =       .25*(qk0(jd,id,2,m)  +qk0(jd1,id,2,m)
     .                  + qk0(jd,id1,2,m) +qk0(jd1,id1,2,m))
        q3k1 =       .25*(qk0(jd,id,3,m)  +qk0(jd1,id,3,m)
     .                  + qk0(jd,id1,3,m) +qk0(jd1,id1,3,m))
        q4k1 =       .25*(qk0(jd,id,4,m)  +qk0(jd1,id,4,m)
     .                  + qk0(jd,id1,4,m) +qk0(jd1,id1,4,m))
        q5k1 = gamma*.25*(qk0(jd,id,5,m)  +qk0(jd1,id,5,m)
     .                  + qk0(jd,id1,5,m) +qk0(jd1,id1,5,m))
        t1k1 = q5k1/q1k1   
c       
c       first cell center location - average in J and I
        q1k2 =       .25*(q(jd,kd,id,1)   +q(jd1,kd,id,1)
     .                  + q(jd,kd,id1,1)  +q(jd1,kd,id1,1))
        q2k2 =       .25*(q(jd,kd,id,2)   +q(jd1,kd,id,2)
     .                  + q(jd,kd,id1,2)  +q(jd1,kd,id1,2))
        q3k2 =       .25*(q(jd,kd,id,3)   +q(jd1,kd,id,3)
     .                  + q(jd,kd,id1,3)  +q(jd1,kd,id1,3))
        q4k2 =       .25*(q(jd,kd,id,4)   +q(jd1,kd,id,4)
     .                  + q(jd,kd,id1,4)  +q(jd1,kd,id1,4))
        q5k2 = gamma*.25*(q(jd,kd,id,5)   +q(jd1,kd,id,5)
     .                  + q(jd,kd,id1,5)  +q(jd1,kd,id1,5))
        t1k2 = q5k2/q1k2    
c       
c       height of first-layer cell
        if (sk(jd ,kdx,id ,4) .eq. 0. .or. 
     +      sk(jd1,kdx,id ,4) .eq. 0. .or.
     +      sk(jd ,kdx,id1,4) .eq. 0. .or.
     +      sk(jd1,kdx,id1,4) .eq. 0.) then
        dx = x(jd,kd+1,i)-x(jd,kd,i)
        dy = y(jd,kd+1,i)-y(jd,kd,i)
        dz = z(jd,kd+1,i)-z(jd,kd,i)
        dn = sqrt(dx**2+dy**2+dz**2)       
        else
        dn = (vol(jd ,kd,id )/sk(jd ,kdx,id ,4)+
     +        vol(jd1,kd,id )/sk(jd1,kdx,id ,4)+
     +        vol(jd ,kd,id1)/sk(jd ,kdx,id1,4)+
     +        vol(jd1,kd,id1)/sk(jd1,kdx,id1,4))/4. 
        end if
c       
c       Get turb viscosity at wall (0 unless wall fn used)
        if (iwf(3) .eq. 1) then
        avgmut =       .25*(vk0(jd,id,1,m)  +vk0(jd1,id,1,m)
     .                    + vk0(jd,id1,1,m) +vk0(jd1,id1,1,m)) 
        else
        avgmut = 0.
        end if
        emuka = (t1k1**1.5)*((1.0+cbar/tinf)/(t1k1+cbar/tinf))
c       
c       Use component of velocity parallel to wall (with 2-point formula)
        urel = q2k2-q2k1   
        vrel = q3k2-q3k1   
        wrel = q4k2-q4k1
        sk1  = (sk(jd ,kdx,id ,1)+sk(jd1,kdx,id ,1)+
     +          sk(jd ,kdx,id1,1)+sk(jd1,kdx,id1,1))/4.
        sk2  = (sk(jd ,kdx,id ,2)+sk(jd1,kdx,id ,2)+
     +          sk(jd ,kdx,id1,2)+sk(jd1,kdx,id1,2))/4.
        sk3  = (sk(jd ,kdx,id ,3)+sk(jd1,kdx,id ,3)+
     +          sk(jd ,kdx,id1,3)+sk(jd1,kdx,id1,3))/4.
        vnorm = urel*sk1 + vrel*sk2 + wrel*sk3
        upar = urel-vnorm*sk1
        vpar = vrel-vnorm*sk2
        wpar = wrel-vnorm*sk3
        Cf   = 2.*sqrt(upar**2 + vpar**2 + wpar**2)
        Cf = 2.0*(emuka+avgmut)/(reue*xmach)*Cf/dn
        if (real(q2k2).lt.0.) cf = -cf
c       Cf vector
        cfx  = 4.*(emuka+avgmut)/(reue*xmach)*upar/dn
        cfy  = 4.*(emuka+avgmut)/(reue*xmach)*vpar/dn
        cfz  = 4.*(emuka+avgmut)/(reue*xmach)*wpar/dn
c       
c       Note:  Ch definition not standard
c       
        eps  = 1.0e-6
        tty1 = 1.0 +gm1*0.5*xmach*xmach
        t1m1 = t1k1-tty1
        if (abs(real(t1m1)).le.real(eps)) then
           Ch = 999.99999
        else
c          2-point formula
        Ch = 2.*(t1k2-t1k1)   
        Ch = (emuka+avgmut)/(reue*pr*(t1k1-tty1))*Ch/dn
        end if
c       
        pres1 = q5k1
        temp1 = t1k1
        yplus = 0.
c       yplus is value at first GRIDPOINT away from wall
        if (abs(real(Cf)).gt.0.) then
           yplus = dn*reue*q1k1*sqrt(ccabs(Cf*0.5/q1k1))/emuka
        end if
c       
c       store data in xv array
c       
        xv(jw,kw,iw,9)  = dn
        xv(jw,kw,iw,10) = Ch
        xv(jw,kw,iw,11) = yplus
        xv(jw,kw,iw,12) = cfx
        xv(jw,kw,iw,13) = cfy
        xv(jw,kw,iw,14) = cfz
c       
 6000   continue
 6001   continue
 6002   continue
c       
        end if
c         
c       output viscous-flow data (skin friction, heat transfer, etc) on
c       i=1 and/or i=idim surfaces
c       
        if (idim.ne.2) then
        if (ivisc(1).ne.0 .and. i1.eq.i2) then
c       
        jw = 0
        do 8002 j=j1,j2,j3
        jw = jw + 1
        jd = min(j,jdim1)
        jd1= max(j-1,1)
        iw = 0
        do 8001 i=i1,i2,i3
        iw = iw + 1
        kw = 0
        do 8000 k=k1,k2,k3
        kw = kw + 1
        kd = min(k,kdim1)
        kd1= max(k-1,1)
        if (i.eq.1) then
           id  = 1
           idx = 1
           m   = 2
        else
           id  = idim1
           idx = idim
           m   = 4
        end if
c       
c       wall value - average in J and k (after QFACE, m=2 & 4 are interface
c       (wall) values
c       
        q1i1 =       .25*(qi0(jd,kd,1,m)   +qi0(jd1,kd,1,m)
     .                  + qi0(jd,kd1,1,m) +qi0(jd1,kd1,1,m))
        q2i1 =       .25*(qi0(jd,kd,2,m)   +qi0(jd1,kd,2,m)
     .                  + qi0(jd,kd1,2,m) +qi0(jd1,kd1,2,m))
        q3i1 =       .25*(qi0(jd,kd,3,m)   +qi0(jd1,kd,3,m)
     .                  + qi0(jd,kd1,3,m) +qi0(jd1,kd1,3,m))
        q4i1 =       .25*(qi0(jd,kd,4,m)   +qi0(jd1,kd,4,m)
     .                  + qi0(jd,kd1,4,m) +qi0(jd1,kd1,4,m))
        q5i1 = gamma*.25*(qi0(jd,kd,5,m)   +qi0(jd1,kd,5,m)
     .                  + qi0(jd,kd1,5,m) +qi0(jd1,kd1,5,m))
        t1i1 = q5i1/q1i1
c       
c       first cell center location - average in J and K
c       
        q1i2 =       .25*(q(jd,kd,id,1)    +q(jd,kd1,id,1)
     .                  + q(jd1,kd,id,1)  +q(jd1,kd1,id,1))
        q2i2 =       .25*(q(jd,kd,id,2)    +q(jd,kd1,id,2)
     .                  + q(jd1,kd,id,2)  +q(jd1,kd1,id,2))
        q3i2 =       .25*(q(jd,kd,id,3)    +q(jd,kd1,id,3)
     .                  + q(jd1,kd,id,3)  +q(jd1,kd1,id,3))
        q4i2 =       .25*(q(jd,kd,id,4)    +q(jd,kd1,id,4)
     .                  + q(jd1,kd,id,4)  +q(jd1,kd1,id,4))
        q5i2 = gamma*.25*(q(jd,kd,id,5)    +q(jd,kd1,id,5)
     .                  + q(jd1,kd,id,5)  +q(jd1,kd1,id,5))
        t1i2 = q5i2/q1i2
c       
        if (si(jd ,kd ,idx,4) .eq. 0. .or.                 
     +      si(jd ,kd1,idx,4) .eq. 0. .or.
     +      si(jd1,kd ,idx,4) .eq. 0. .or.
     +      si(jd1,kd1,idx,4) .eq. 0.) then
        dx = x(jd,kd,id+1)-x(jd,kd,id)
        dy = y(jd,kd,id+1)-y(jd,kd,id)
        dz = z(jd,kd,id+1)-z(jd,kd,id)
        dn = sqrt(dx**2+dy**2+dz**2)
        else
        dn = (vol(jd ,kd ,id)/si(jd ,kd ,idx,4)+
     +        vol(jd ,kd1,id)/si(jd ,kd1,idx,4)+
     +        vol(jd1,kd ,id)/si(jd1,kd ,idx,4)+
     +        vol(jd1,kd1,id)/si(jd1,kd1,idx,4))/4.
        end if
c       
c       Get turb viscosity at wall (0 unless wall fn used)
        if (iwf(1) .eq. 1) then
        avgmut =       .25*(vi0(jd,kd,1,m)  +vi0(jd1,kd,1,m)
     .                    + vi0(jd,kd1,1,m) +vi0(jd1,kd1,1,m))
        else
        avgmut = 0.
        end if
        emuka = (t1i1**1.5)*((1.0+cbar/tinf)/(t1i1+cbar/tinf))
c       
c       Use component of velocity parallel to wall
        urel = q2i2-q2i1
        vrel = q3i2-q3i1
        wrel = q4i2-q4i1
        si1  = (si(jd ,kd ,idx,1)+si(jd ,kd1,idx,1)+
     +          si(jd1,kd ,idx,1)+si(jd1,kd1,idx,1))/4.
        si2  = (si(jd ,kd ,idx,2)+si(jd ,kd1,idx,2)+
     +          si(jd1,kd ,idx,2)+si(jd1,kd1,idx,2))/4.
        si3  = (si(jd ,kd ,idx,3)+si(jd ,kd1,idx,3)+
     +          si(jd1,kd ,idx,3)+si(jd1,kd1,idx,3))/4.
        vnorm = urel*si1 + vrel*si2 + wrel*si3
        upar = urel-vnorm*si1
        vpar = vrel-vnorm*si2
        wpar = wrel-vnorm*si3
        Cf   = 2.*sqrt(upar**2 + vpar**2 + wpar**2)
        Cf = 2.0*(emuka+avgmut)/(reue*xmach)*Cf/dn
        if (real(q2i2).lt.0.) cf = -cf
c       Cf vector
        cfx  = 4.*(emuka+avgmut)/(reue*xmach)*upar/dn
        cfy  = 4.*(emuka+avgmut)/(reue*xmach)*vpar/dn
        cfz  = 4.*(emuka+avgmut)/(reue*xmach)*wpar/dn
c       
c       Note:  Ch definition not standard
c       
        eps  = 1.0e-6
        tty1 = 1.0 +gm1*0.5*xmach*xmach
        t1m1 = t1i1-tty1
        if (abs(real(t1m1)).le.real(eps)) then
           Ch = 999.99999
        else
c          2-point formula
           Ch = 2.*(t1i2-t1i1)
           Ch = (emuka+avgmut)/(reue*pr*(t1i1-tty1))*Ch/dn
        end if
c       
        pres1 = q5i1
        temp1 = t1i1
        yplus = 0.
c       yplus is valus at first GRIDPOINT away from wall
        if (abs(real(Cf)).gt.0.) then
           yplus = dn*reue*q1i1*sqrt(ccabs(Cf*0.5/q1i1))/emuka
        end if
c       
c       store data in xv array
c       
        xv(jw,kw,iw,9)  = dn
        xv(jw,kw,iw,10) = Ch
        xv(jw,kw,iw,11) = yplus
        xv(jw,kw,iw,12) = cfx
        xv(jw,kw,iw,13) = cfy
        xv(jw,kw,iw,14) = cfz
c       
 8000   continue
 8001   continue
 8002   continue
c       
        end if
        end if
c
c       off-body wall paralell velocity at grid points (first-layer cell height)
c       easier to plot wall surface-restricted streamlines
c       
        if(ioffbody.eq.1) then
          ldw = 5
          if (i1 .eq. i2 .and. i1 .eq. 1) then
            call cctogp(jdim,kdim,idim,i1+1,i2+1,i3,j1,j2,j3,k1,k2,k3,
     .                q,qi0,qj0,qk0,jdw,kdw,idw,xv(1,1,1,16),ldw)
          end if
          if (i1 .eq. i2 .and. i1 .eq. idim) then
            call cctogp(jdim,kdim,idim,i1-1,i2-1,i3,j1,j2,j3,k1,k2,k3,
     .                q,qi0,qj0,qk0,jdw,kdw,idw,xv(1,1,1,16),ldw)
          end if
          if (j1 .eq. j2 .and. j1 .eq. 1) then
            call cctogp(jdim,kdim,idim,i1,i2,i3,j1+1,j2+1,j3,k1,k2,k3,
     .                q,qi0,qj0,qk0,jdw,kdw,idw,xv(1,1,1,16),ldw)
          end if
          if (j1 .eq. j2 .and. j1 .eq. jdim) then
            call cctogp(jdim,kdim,idim,i1,i2,i3,j1-1,j2-1,j3,k1,k2,k3,
     .                q,qi0,qj0,qk0,jdw,kdw,idw,xv(1,1,1,16),ldw)
          end if
          if (k1 .eq. k2 .and. k1 .eq. 1) then
            call cctogp(jdim,kdim,idim,i1,i2,i3,j1,j2,j3,k1+1,k2+1,k3,
     .                q,qi0,qj0,qk0,jdw,kdw,idw,xv(1,1,1,16),ldw)
          end if
          if (k1 .eq. k2 .and. k1 .eq. kdim) then
            call cctogp(jdim,kdim,idim,i1,i2,i3,j1,j2,j3,k1-1,k2-1,k3,
     .                q,qi0,qj0,qk0,jdw,kdw,idw,xv(1,1,1,16),ldw)
          end if
c         
          if (j1 .eq. j2) then
            ix=0
            do i=i1,i2,i3
              ix=ix+1
              id1=max(i-1,1)
              id =min(i,idim1)
              kx=0
              do k=k1,k2,k3
                kx=kx+1
                kd1=max(k-1,1)
                kd =min(k,kdim1)
                jx=0
                do j=j1,j2,j3
                  jx=jx+1
                  sjm1=(sj(j  ,kd ,id ,1)+sj(j  ,kd ,id1,1)+
     +                  sj(j  ,kd1,id ,1)+sj(j  ,kd1,id1,1))*.25
                  sjm2=(sj(j  ,kd ,id ,2)+sj(j  ,kd ,id1,2)+
     +                  sj(j  ,kd1,id ,2)+sj(j  ,kd1,id1,2))*.25
                  sjm3=(sj(j  ,kd ,id ,3)+sj(j  ,kd ,id1,3)+
     +                  sj(j  ,kd1,id ,3)+sj(j  ,kd1,id1,3))*.25
                  vnorm = xv(jx,kx,ix,17)*sjm1 +
     +                    xv(jx,kx,ix,18)*sjm2 +
     +                    xv(jx,kx,ix,19)*sjm3
                  xv(jx,kx,ix,1) = (xv(jx,kx,ix,17)-vnorm*sjm1)/xmach
                  xv(jx,kx,ix,2) = (xv(jx,kx,ix,18)-vnorm*sjm2)/xmach
                  xv(jx,kx,ix,3) = (xv(jx,kx,ix,19)-vnorm*sjm3)/xmach
                enddo
              enddo
            end do
          else if (k1 .eq. k2) then
            ix=0
            do i=i1,i2,i3
              ix=ix+1
              id1=max(i-1,1)
              id =min(i,idim1)
              kx=0
              do k=k1,k2,k3
                kx=kx+1
                jx=0
                do j=j1,j2,j3
                  jx=jx+1
                  jd1=max(j-1,1)
                  jd =min(j,jdim1)
                  skm1=(sk(jd ,k  ,id ,1)+sk(jd ,k  ,id1,1)+
     +                  sk(jd1,k  ,id ,1)+sk(jd1,k  ,id1,1))*.25
                  skm2=(sk(jd ,k  ,id ,2)+sk(jd ,k  ,id1,2)+
     +                  sk(jd1,k  ,id ,2)+sk(jd1,k  ,id1,2))*.25
                  skm3=(sk(jd ,k  ,id ,3)+sk(jd ,k  ,id1,3)+
     +                  sk(jd1,k  ,id ,3)+sk(jd1,k  ,id1,3))*.25
                  vnorm = xv(jx,kx,ix,17)*skm1 +
     +                    xv(jx,kx,ix,18)*skm2 +
     +                    xv(jx,kx,ix,19)*skm3
                  xv(jx,kx,ix,1) = (xv(jx,kx,ix,17)-vnorm*skm1)/xmach
                  xv(jx,kx,ix,2) = (xv(jx,kx,ix,18)-vnorm*skm2)/xmach
                  xv(jx,kx,ix,3) = (xv(jx,kx,ix,19)-vnorm*skm3)/xmach
                enddo
              enddo
            enddo
          else if (i1 .eq. i2) then
            ix=0
            do i=i1,i2,i3
              ix=ix+1
              kx=0
              do k=k1,k2,k3
                kx=kx+1
                kd1=max(k-1,1)
                kd =min(k,kdim1)
                jx=0
                do j=j1,j2,j3
                  jx=jx+1
                  jd1=max(j-1,1)
                  jd =min(j,jdim1)
                  sim1=(si(jd ,kd ,i  ,1)+si(jd1,kd ,i  ,1)+
     +                  si(jd ,kd1,i  ,1)+si(jd1,kd1,i  ,1))*.25
                  sim2=(si(jd ,kd ,i  ,2)+si(jd1,kd ,i  ,2)+
     +                  si(jd ,kd1,i  ,2)+si(jd1,kd1,i  ,2))*.25
                  sim3=(si(jd ,kd ,i  ,3)+si(jd1,kd ,i  ,3)+
     +                  si(jd ,kd1,i  ,3)+si(jd1,kd1,i  ,3))*.25
                  vnorm = xv(jx,kx,ix,17)*sim1 +
     +                    xv(jx,kx,ix,18)*sim2 +
     +                    xv(jx,kx,ix,19)*sim3
                  xv(jx,kx,ix,1) = (xv(jx,kx,ix,17)-vnorm*sim1)/xmach
                  xv(jx,kx,ix,2) = (xv(jx,kx,ix,18)-vnorm*sim2)/xmach
                  xv(jx,kx,ix,3) = (xv(jx,kx,ix,19)-vnorm*sim3)/xmach
                enddo
              enddo
            enddo 
          end if
        end if
c
c       install average surface solutions into vsj0,vsjdim,vsk0...
c       the storage array size is 20, but only the first 15 is used in fact
c       
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
              vsj0(k,i,l,1)=(vsj0(k,i,l,1)*(xnumavg-1.)+
     +          xv(jw,kw,iw,l))/xnumavg
              vsj0(k,i,l,2)=(vsj0(k,i,l,2)*(xnumavg-1.)+
     +          xv(jw,kw,iw,l)*xv(jw,kw,iw,l))/xnumavg
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
              vsjdim(k,i,l,1)=(vsjdim(k,i,l,1)*(xnumavg-1.)+
     +          xv(jw,kw,iw,l))/xnumavg
              vsjdim(k,i,l,2)=(vsjdim(k,i,l,2)*(xnumavg-1.)+
     +          xv(jw,kw,iw,l)*xv(jw,kw,iw,l))/xnumavg
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
              vsk0(j,i,l,1)=(vsk0(j,i,l,1)*(xnumavg-1.)+
     +          xv(jw,kw,iw,l))/xnumavg
              vsk0(j,i,l,2)=(vsk0(j,i,l,2)*(xnumavg-1.)+
     +          xv(jw,kw,iw,l)*xv(jw,kw,iw,l))/xnumavg
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
              vskdim(j,i,l,1)=(vskdim(j,i,l,1)*(xnumavg-1.)+
     +          xv(jw,kw,iw,l))/xnumavg
              vskdim(j,i,l,2)=(vskdim(j,i,l,2)*(xnumavg-1.)+
     +          xv(jw,kw,iw,l)*xv(jw,kw,iw,l))/xnumavg
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
              vsi0(j,k,l,1)=(vsi0(j,k,l,1)*(xnumavg-1.)+
     +          xv(jw,kw,iw,l))/xnumavg
              vsi0(j,k,l,2)=(vsi0(j,k,l,2)*(xnumavg-1.)+
     +          xv(jw,kw,iw,l)*xv(jw,kw,iw,l))/xnumavg
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
              vsidim(j,k,l,1)=(vsidim(j,k,l,1)*(xnumavg-1.)+
     +          xv(jw,kw,iw,l))/xnumavg
              vsidim(j,k,l,2)=(vsidim(j,k,l,2)*(xnumavg-1.)+
     +          xv(jw,kw,iw,l)*xv(jw,kw,iw,l))/xnumavg
            end do
            end do
            end do
          end do
c
        end if
c
      end if
      return
      end
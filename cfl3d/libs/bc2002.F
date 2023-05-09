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
      subroutine bc2002(jdim,kdim,idim,q,qj0,qk0,qi0,sj,sk,si,bcj,bck,
     .                  bci,xtbj,xtbk,xtbi,atbj,atbk,atbi,ista,iend,
     .                  jsta,jend,ksta,kend,nface,tursav,tj0,tk0,
     .                  ti0,vist3d,vj0,vk0,vi0,mdim,ndim,bcdata,
     .                  filname,iuns,nou,bou,nbuf,ibufdim,myid,nummem,
     .                  x,y,z)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Set pressure ratio; extrapolate other quantities
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
      character*80 filname,char
c
      dimension nou(nbuf)
      dimension q(jdim,kdim,idim,5), qi0(jdim,kdim,5,4),
     .          qj0(kdim,idim-1,5,4),qk0(jdim,idim-1,5,4),
     .          x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim)  
      dimension sk(jdim,kdim,idim-1,5),si(jdim,kdim,idim,5),
     .          sj(jdim,kdim,idim-1,5)
      dimension bcj(kdim,idim-1,2),bck(jdim,idim-1,2),bci(jdim,kdim,2)
      dimension xtbj(kdim,idim-1,3,2),xtbk(jdim,idim-1,3,2),
     .          xtbi(jdim,kdim,3,2),atbj(kdim,idim-1,3,2),
     .          atbk(jdim,idim-1,3,2),atbi(jdim,kdim,3,2)
      dimension bcdata(mdim,ndim,2,12)
      dimension tursav(jdim,kdim,idim,nummem),tj0(kdim,idim-1,nummem,4),
     .          tk0(jdim,idim-1,nummem,4),ti0(jdim,kdim,nummem,4),
     .          vj0(kdim,idim-1,1,4),vk0(jdim,idim-1,1,4),
     .          vi0(jdim,kdim,1,4),vist3d(jdim,kdim,idim)
c
      common /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),
     .        tur10(7)
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
      common /reyue/ reue,tinf,ivisc(3)
      common /sklton/ isklton
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!张宇飞用来定义进出口边界 zyf
      dimension dataoutlet(10000,10),temp(2)
      common /outletfile/ outletcount,ndataoutlet
      common /outletdata/ dataoutlet
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!定义流量控制边界
      common /massflow/ xmassflow_now(10) , dpratio(10),ptotratio(10)   !zyf
      common /mass_group/      xmassflow_target(10)    !zyf      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!进出口剖面边界
      twopi = 8.0*atan(1.0)
      if(abs(bcdata(1,1,1,2)-10000.0) < 1e-2 )then        !如果bcdata(1,1,1,2)=10000，那么表示给定出口剖面速度
          outletcount = outletcount + 1
          if(outletcount==1)then
              open(51,file='outlet.dat')
              read(51,*)char,char,ndataoutlet
              do idata=1,ndataoutlet
                  read(51,*)dataoutlet(idata,1:2)  
                ! 1到2的输入数据分别为半径、静压
              enddo
              close(51)
          endif
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
c
      jend1 = jend-1
      kend1 = kend-1
      iend1 = iend-1
c
c     this bc makes use of only one plane of data 
c
      ip    = 1
c
c            * * * * * * * * * * * * * * * * * * * * * *
c            * standard boundary condition bctype=2002 *
c            * * * * * * * * * * * * * * * * * * * * * *
c
c******************************************************************************
c      j=1 boundary        set pressure, extrapolate others         bctype 2002
c******************************************************************************
c
      if (nface.eq.3) then
c
      do 100 l=1,4
      do 100 i=ista,iend1
      do 100 k=ksta,kend1
      qj0(k,i,l,1) = q(1,k,i,l)
      qj0(k,i,l,2) = qj0(k,i,l,1)
      bcj(k,i,1) = 0.0
  100 continue
      do 102 i=ista,iend1
      ii = i-ista+1
      do 102 k=ksta,kend1
      kk = k-ksta+1
      pratio = bcdata(kk,ii,ip,1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!进出口剖面边界 zyf
      if(abs(bcdata(kk,ii,ip,2)-10000.0) < 1e-2 )then
          xc =0.25*(x(1,k,i)+x(1,k,i+1)+x(1,k+1,i)+x(1,k+1,i+1))
          yc =0.25*(y(1,k,i)+y(1,k,i+1)+y(1,k+1,i)+y(1,k+1,i+1))
          zc =0.25*(z(1,k,i)+z(1,k,i+1)+z(1,k+1,i)+z(1,k+1,i+1))
          radius = sqrt(yc*yc + zc*zc)
          do idata=1,ndataoutlet-1
              r1=abs(dataoutlet(idata,1))
              r2=abs(dataoutlet(idata+1,1))
              if ((radius>=r1 .and. radius<=r2) .or.
     .            (idata  ==1 .and. radius< r1) .or. 
     .   (idata==ndataoutlet-1 .and. radius> r2))then
                  temp(1)=radius
                  do itemp=2,2
       temp(itemp)=(dataoutlet(idata+1,itemp)-dataoutlet(idata,itemp))
     .                /(r2-r1)*(radius-r1)+dataoutlet(idata,itemp)
                  end do
                  cycle
              endif
          end do
          pratio = temp(2)
          !write(*,*)xc,yc,zc,radius,pratio
          !pause
      elseif(bcdata(kk,ii,ip,2).ne.0)then        !控制组的流量控制设置
          igroup=int(bcdata(kk,ii,ip,2))
          if(igroup.ge.1)then
          if(abs(xmassflow_target(igroup))>1e-8)then
          pratio=pratio+ptotratio(igroup)
          endif
          endif
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      pset   = pratio*p0
      qj0(k,i,5,1) = pset
      qj0(k,i,5,2) = pset
  102 continue
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 191 i=ista,iend1
        do 191 k=ksta,kend1
          vj0(k,i,1,1) = vist3d(1,k,i)
          vj0(k,i,1,2) = vist3d(1,k,i)
  191   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do l=1,nummem
        do 101 i=ista,iend1
        do 101 k=ksta,kend1
          tj0(k,i,l,1) = tursav(1,k,i,l)
          tj0(k,i,l,2) = tursav(1,k,i,l)
  101   continue
        enddo
      end if
      end if
      end if
c
c******************************************************************************
c      j=jdim boundary        set pressure, extrapolate others      bctype 2002
c******************************************************************************
c 
      if (nface.eq.4) then
c
      do 200 l=1,4
      do 200 i=ista,iend1
      do 200 k=ksta,kend1
      qj0(k,i,l,3) = q(jdim1,k,i,l)
      qj0(k,i,l,4) = qj0(k,i,l,3)
      bcj(k,i,2) = 0.0
  200 continue
      do 202 i=ista,iend1
      ii = i-ista+1
      do 202 k=ksta,kend1
      kk = k-ksta+1
      pratio = bcdata(kk,ii,ip,1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!进出口剖面边界zyf
      if(abs(bcdata(kk,ii,ip,2)-10000.0) < 1e-2 )then     
      xc =0.25*(x(jdim,k,i)+x(jdim,k,i+1)+x(jdim,k+1,i)+x(jdim,k+1,i+1))
      yc =0.25*(y(jdim,k,i)+y(jdim,k,i+1)+y(jdim,k+1,i)+y(jdim,k+1,i+1))
      zc =0.25*(z(jdim,k,i)+z(jdim,k,i+1)+z(jdim,k+1,i)+z(jdim,k+1,i+1))
          radius = sqrt(yc*yc + zc*zc)
          do idata=1,ndataoutlet-1
              r1=abs(dataoutlet(idata,1))
              r2=abs(dataoutlet(idata+1,1))
              if ((radius>=r1 .and. radius<=r2) .or.
     .            (idata  ==1 .and. radius< r1) .or. 
     .   (idata==ndataoutlet-1 .and. radius> r2))then
                  temp(1)=radius
                  do itemp=2,2
       temp(itemp)=(dataoutlet(idata+1,itemp)-dataoutlet(idata,itemp))
     .                /(r2-r1)*(radius-r1)+dataoutlet(idata,itemp)
                  end do
                  cycle
              endif
          end do
          pratio = temp(2)
      elseif(bcdata(kk,ii,ip,2).ne.0)then        !控制组的流量控制设置
          igroup=int(bcdata(kk,ii,ip,2))
          if(igroup.ge.1)then
          if(abs(xmassflow_target(igroup))>1e-8)then
          pratio=pratio+ptotratio(igroup)
          endif
          endif
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      pset   = pratio*p0
      qj0(k,i,5,3) = pset
      qj0(k,i,5,4) = pset
  202 continue
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 291 i=ista,iend1
        do 291 k=ksta,kend1
          vj0(k,i,1,3) = vist3d(jdim1,k,i)
          vj0(k,i,1,4) = vist3d(jdim1,k,i)
  291   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do l=1,nummem
        do 201 i=ista,iend1
        do 201 k=ksta,kend1
          tj0(k,i,l,3) = tursav(jdim1,k,i,l)
          tj0(k,i,l,4) = tursav(jdim1,k,i,l)
  201   continue
        enddo
      end if
      end if
      end if
c
c******************************************************************************
c      k=1 boundary        set pressure, extrapolate others         bctype 2002
c******************************************************************************
c
      if (nface.eq.5) then
c
      do 300 l=1,4
      do 300 i=ista,iend1
      do 300 j=jsta,jend1
      qk0(j,i,l,1) = q(j,1,i,l)
      qk0(j,i,l,2) = qk0(j,i,l,1)
      bck(j,i,1) = 0.0
  300 continue
      do 302 i=ista,iend1
      ii = i-ista+1
      do 302 j=jsta,jend1
      jj = j-jsta+1
      pratio = bcdata(jj,ii,ip,1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!进出口剖面边界zyf
      if(abs(bcdata(jj,ii,ip,2)-10000.0) < 1e-2 )then     
      xc =0.25*(x(j,1,i)+x(j,1,i+1)+x(j+1,1,i)+x(j+1,1,i+1))
      yc =0.25*(y(j,1,i)+y(j,1,i+1)+y(j+1,1,i)+y(j+1,1,i+1))
      zc =0.25*(z(j,1,i)+z(j,1,i+1)+z(j+1,1,i)+z(j+1,1,i+1))
          radius = sqrt(yc*yc + zc*zc)
          do idata=1,ndataoutlet-1
              r1=abs(dataoutlet(idata,1))
              r2=abs(dataoutlet(idata+1,1))
              if ((radius>=r1 .and. radius<=r2) .or.
     .            (idata  ==1 .and. radius< r1) .or. 
     .   (idata==ndataoutlet-1 .and. radius> r2))then
                  temp(1)=radius
                  do itemp=2,2
       temp(itemp)=(dataoutlet(idata+1,itemp)-dataoutlet(idata,itemp))
     .                /(r2-r1)*(radius-r1)+dataoutlet(idata,itemp)
                  end do
                  cycle
              endif
          end do
          pratio = temp(2)
      elseif(bcdata(jj,ii,ip,2).ne.0)then        !控制组的流量控制设置
          igroup=int(bcdata(jj,ii,ip,2))
          if(igroup.ge.1)then
          if(abs(xmassflow_target(igroup))>1e-8)then
          pratio=pratio+ptotratio(igroup)
          endif
          endif
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      pset   = pratio*p0
      qk0(j,i,5,1) = pset
      qk0(j,i,5,2) = pset
  302 continue
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 391 i=ista,iend1
        do 391 j=jsta,jend1
          vk0(j,i,1,1) = vist3d(j,1,i)
          vk0(j,i,1,2) = vist3d(j,1,i)
  391   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do l=1,nummem
        do 301 i=ista,iend1
        do 301 j=jsta,jend1
          tk0(j,i,l,1) = tursav(j,1,i,l)
          tk0(j,i,l,2) = tursav(j,1,i,l)
  301   continue
        enddo
      end if
      end if
      end if
c
c******************************************************************************
c      k=kdim boundary        set pressure, extrapolate others      bctype 2002
c******************************************************************************
c
      if (nface.eq.6) then
c
      do 400 l=1,4
      do 400 i=ista,iend1
      do 400 j=jsta,jend1
      qk0(j,i,l,3) = q(j,kdim1,i,l)
      qk0(j,i,l,4) = qk0(j,i,l,3)
      bck(j,i,2) = 0.0
  400 continue
      do 402 i=ista,iend1
      ii = i-ista+1
      do 402 j=jsta,jend1
      jj = j-jsta+1
      pratio = bcdata(jj,ii,ip,1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!进出口剖面边界zyf
      if(abs(bcdata(jj,ii,ip,2)-10000.0) < 1e-2 )then     
      xc =0.25*(x(j,kdim,i)+x(j,kdim,i+1)+x(j+1,kdim,i)+x(j+1,kdim,i+1))
      yc =0.25*(y(j,kdim,i)+y(j,kdim,i+1)+y(j+1,kdim,i)+y(j+1,kdim,i+1))
      zc =0.25*(z(j,kdim,i)+z(j,kdim,i+1)+z(j+1,kdim,i)+z(j+1,kdim,i+1))
          radius = sqrt(yc*yc + zc*zc)
          do idata=1,ndataoutlet-1
              r1=abs(dataoutlet(idata,1))
              r2=abs(dataoutlet(idata+1,1))
              if ((radius>=r1 .and. radius<=r2) .or.
     .            (idata  ==1 .and. radius< r1) .or. 
     .   (idata==ndataoutlet-1 .and. radius> r2))then
                  temp(1)=radius
                  do itemp=2,2
       temp(itemp)=(dataoutlet(idata+1,itemp)-dataoutlet(idata,itemp))
     .                /(r2-r1)*(radius-r1)+dataoutlet(idata,itemp)
                  end do
                  cycle
              endif
          end do
          pratio = temp(2)
      elseif(bcdata(jj,ii,ip,2).ne.0)then        !控制组的流量控制设置
          igroup=int(bcdata(jj,ii,ip,2))
          !write(*,*)igroup,xmassflow_target(igroup),ptotratio(igroup)
          if(igroup.ge.1)then
          if(abs(xmassflow_target(igroup))>1e-8)then
          pratio=pratio+ptotratio(igroup)
          endif
          endif
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      pset   = pratio*p0
      qk0(j,i,5,3) = pset
      qk0(j,i,5,4) = pset
  402 continue
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 491 i=ista,iend1
        do 491 j=jsta,jend1
          vk0(j,i,1,3) = vist3d(j,kdim1,i)
          vk0(j,i,1,4) = vist3d(j,kdim1,i)
  491   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do l=1,nummem
        do 401 i=ista,iend1
        do 401 j=jsta,jend1
          tk0(j,i,l,3) = tursav(j,kdim1,i,l)
          tk0(j,i,l,4) = tursav(j,kdim1,i,l)
  401   continue
        enddo
      end if
      end if
      end if
c
c******************************************************************************
c      i=1 boundary        set pressure, extrapolate others         bctype 2002
c******************************************************************************
c
      if (nface.eq.1) then
c
      do 500 l=1,4
      do 500 k=ksta,kend1
      do 500 j=jsta,jend1
      qi0(j,k,l,1) = q(j,k,1,l)
      qi0(j,k,l,2) = qi0(j,k,l,1)
      bci(j,k,1) = 0.0
  500 continue
      do 502 k=ksta,kend1
      kk = k-ksta+1
      do 502 j=jsta,jend1
      jj = j-jsta+1
      pratio = bcdata(jj,kk,ip,1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!进出口剖面边界zyf
      if(abs(bcdata(jj,kk,ip,2)-10000.0) < 1e-2 )then     
          xc =0.25*(x(j,k,1)+x(j+1,k,1)+x(j,k+1,1)+x(j+1,k+1,1))
          yc =0.25*(y(j,k,1)+y(j+1,k,1)+y(j,k+1,1)+y(j+1,k+1,1))
          zc =0.25*(z(j,k,1)+z(j+1,k,1)+z(j,k+1,1)+z(j+1,k+1,1))
          radius = sqrt(yc*yc + zc*zc)
          do idata=1,ndataoutlet-1
              r1=abs(dataoutlet(idata,1))
              r2=abs(dataoutlet(idata+1,1))
              if ((radius>=r1 .and. radius<=r2) .or.
     .            (idata  ==1 .and. radius< r1) .or. 
     .   (idata==ndataoutlet-1 .and. radius> r2))then
                  temp(1)=radius
                  do itemp=2,2
       temp(itemp)=(dataoutlet(idata+1,itemp)-dataoutlet(idata,itemp))
     .                /(r2-r1)*(radius-r1)+dataoutlet(idata,itemp)
                  end do
                  cycle
              endif
          end do
          pratio = temp(2)
      elseif(bcdata(jj,kk,ip,2).ne.0)then        !控制组的流量控制设置
          igroup=int(bcdata(jj,kk,ip,2))
          if(igroup.ge.1)then
          if(abs(xmassflow_target(igroup))>1e-8)then
          pratio=pratio+ptotratio(igroup)
          endif
          endif
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      pset   = pratio*p0
      qi0(j,k,5,1) = pset
      qi0(j,k,5,2) = pset
  502 continue
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 591 k=ksta,kend1
        do 591 j=jsta,jend1
          vi0(j,k,1,1) = vist3d(j,k,1)
          vi0(j,k,1,2) = vist3d(j,k,1)
  591   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do l=1,nummem
        do 501 k=ksta,kend1
        do 501 j=jsta,jend1
          ti0(j,k,l,1) = tursav(j,k,1,l)
          ti0(j,k,l,2) = tursav(j,k,1,l)
  501   continue
        enddo
      end if
      end if
      end if
c
c******************************************************************************
c      i=idim boundary        set pressure, extrapolate others      bctype 2002
c******************************************************************************
c
      if (nface.eq.2) then
c
      do 600 l=1,4
      do 600 k=ksta,kend1
      do 600 j=jsta,jend1
      qi0(j,k,l,3) = q(j,k,idim1,l)
      qi0(j,k,l,4) = qi0(j,k,l,3)
      bci(j,k,2) = 0.0
  600 continue
      do 602 k=ksta,kend1
      kk = k-ksta+1
      do 602 j=jsta,jend1
      jj = j-jsta+1
      pratio = bcdata(jj,kk,ip,1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!进出口剖面边界zyf
      if(abs(bcdata(jj,kk,ip,2)-10000.0) < 1e-2 )then     
      xc =0.25*(x(j,k,idim)+x(j+1,k,idim)+x(j,k+1,idim)+x(j+1,k+1,idim))
      yc =0.25*(y(j,k,idim)+y(j+1,k,idim)+y(j,k+1,idim)+y(j+1,k+1,idim))
      zc =0.25*(z(j,k,idim)+z(j+1,k,idim)+z(j,k+1,idim)+z(j+1,k+1,idim))
          radius = sqrt(yc*yc + zc*zc)
          do idata=1,ndataoutlet-1
              r1=abs(dataoutlet(idata,1))
              r2=abs(dataoutlet(idata+1,1))
              if ((radius>=r1 .and. radius<=r2) .or.
     .            (idata  ==1 .and. radius< r1) .or. 
     .   (idata==ndataoutlet-1 .and. radius> r2))then
                  temp(1)=radius
                  do itemp=2,2
       temp(itemp)=(dataoutlet(idata+1,itemp)-dataoutlet(idata,itemp))
     .                /(r2-r1)*(radius-r1)+dataoutlet(idata,itemp)
                  end do
                  cycle
              endif
          end do
          pratio = temp(2)
      elseif(bcdata(jj,kk,ip,2).ne.0)then        !控制组的流量控制设置
          igroup=int(bcdata(jj,kk,ip,2))
          if(igroup.ge.1)then
          if(abs(xmassflow_target(igroup))>1e-8)then
          pratio=pratio+ptotratio(igroup)
          endif
          endif
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      pset   = pratio*p0
      qi0(j,k,5,3) = pset
      qi0(j,k,5,4) = pset
  602 continue
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 691 k=ksta,kend1
        do 691 j=jsta,jend1
          vi0(j,k,1,3) = vist3d(j,k,idim1)
          vi0(j,k,1,4) = vist3d(j,k,idim1)
  691   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do l=1,nummem
        do 601 k=ksta,kend1
        do 601 j=jsta,jend1
          ti0(j,k,l,3) = tursav(j,k,idim1,l)
          ti0(j,k,l,4) = tursav(j,k,idim1,l)
  601   continue
        enddo
      end if
      end if
      end if
c
      return 
      end

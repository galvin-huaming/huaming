! This file is part of GTC version 4.5 
! GTC version 4.5 is released under the 3-Clause BSD license:

! Copyright (c) 2002,2010,2016, GTC Team (team leader: Zhihong Lin, zhihongl@uci.edu)
! All rights reserved.

! Redistribution and use in source and binary forms, with or without 
! modification, are permitted provided that the following conditions are met:

! 1. Redistributions of source code must retain the above copyright notice, 
!    this list of conditions and the following disclaimer.

! 2. Redistributions in binary form must reproduce the above copyright notice, 
!    this list of conditions and the following disclaimer in the documentation 
!    and/or other materials provided with the distribution.

! 3. Neither the name of the GTC Team nor the names of its contributors may be 
!    used to endorse or promote products derived from this software without 
!    specific prior written permission.
! ==============================================================================

subroutine snapshot
  use system_env, only: MPI_SUM, MPI_COMM_WORLD, ierror
  use utility, only: check_allocation, mstat
  use precision, only: mpi_Rsize, lk
  use global_parameters
  use particle_array
  use field_array
  use equilibrium
  use magnetic_island
  use conservative_scheme,only:sapara_na
#ifdef ADIOS2
  use adios2
  use adios2_common
#endif
  !use physical_meaning, only: nfield => snapshot_nfield
#ifdef _SDP
  use sdp, only: sdp_perturbation_output
#endif
  use spline_function, only: spx, spz, spb, spcpsi, spgpsi, spd, spq, dxdp, dxdt, dzdp, dzdt, dots
  implicit none
  
  integer,parameter :: nfield=3,nvgrid=65,iosnap=222,iophi=223,iob=224
  integer mtgrid,nf,ns,m,ip,jenergy,kpitch,lambdabin,nsnap,i,j,ij,jt,icount,isp,jst,ii
  integer mtgrid_tmp, mpsi_tmp, indx
  real(lk) upara,delf,rmu,fullf,ai_inv,vi_inv,af_inv,vf_inv,ae_inv,ve_inv,&
       afe_inv,vfe_inv,emax_inv,delr,rg,b,&
       emax_inv_tmp,lambmax_inv,sqdet_inv,gradapara_phi,&
       energy,pitch,wt,tdum,pdum,dpx,dp2,dtx,dt2,psitmp,thetatmp,ztmp,&
       profile(0:mpsi,6,nspecies),pdf(nvgrid,4,nspecies),q,gpp,gpt,gtt,gzz,tmp1,tmp,&
       dxdptmp,dxdttmp,dzdptmp,dzdttmp,g,ri,dzx,ksz,dx(spdim),majorr,& 
       dfield(mgrid),proftmp(0:mpsi,6,nspecies),pdftmp(nvgrid,4,nspecies),niflux(0:mpsi),eqmeshni(0:mpsi),sbratio,sphiratio,dkne
  real(lk) pdf2d(nvgrid,nvgrid,2,nspecies),pdf2dtmp(nvgrid,nvgrid,2,nspecies)
  real(lk),dimension(:),allocatable :: eachflux,eachrz,bpararms,bperprms,bratio,bperpdiag,bparadiag,phieff,phieffrms,phirms,phiratio
  real(lk),dimension(:),allocatable :: density_adios
  real(lk),dimension(:,:),allocatable :: allflux,allrz
  real(lk),dimension(:,:,:),allocatable :: poloidata,fluxdata
  character(len=64) cdum0,cdum1,cdum2,cdum3
  integer :: ierr

  ! adios2 vars
#ifdef ADIOS2
  type(adios2_variable), save :: var_step_id, var_nspecies, var_nfield, var_nvgrid, var_mpsi, var_mtgrid, var_mtoroidal, var_emax_inv
  type(adios2_variable), save :: var_profile, var_pdf, var_poloidata, var_fluxdata, var_density
#endif
  integer*8 :: shape3(3)

! particle species: 1=ion, 2=electron, 3=EP
! radial profiles: density, flow, energy
  profile=0.0_lk

! distribution function: energy, pitch angle
  pdf=0.0_lk
  pdf2d=0.0_lk

! species ns=1: thermal ion
  ns=1
  ai_inv=1.0_lk/aion
  vi_inv=sqrt(aion)/rho0 !parallel velocity normalized by c_s
  delr=1.0_lk/deltar
  emax_inv=0.2_lk
  lambmax_inv=0.8_lk !so the maximum lambda=1/0.8=1.25


!!$omp parallel do private(i,pdum,isp,dpx,dp2)
  do i=0,mpsi
     pdum=psimesh(i)
     isp=max(1,min(lsp-1,ceiling(pdum*spdpsi_inv)))
     dpx=pdum-spdpsi*real(isp-1)
     dp2=dpx*dpx

     eqmeshni(i)=nipp(1,isp)+nipp(2,isp)*dpx+nipp(3,isp)*dp2
 enddo
  if(fullki==0)then  
  do m=1,mi
     psitmp=zion(1,m)
     isp=max(1,min(lsp-1,ceiling(psitmp*spdpsi_inv)))
     dpx=psitmp-spdpsi*real(isp-1)
! radaial spline of rg & b avoids sigularity near axis: y=y1+y2*sqrt(x)+y3*x for psi<spdpsi
     if(isp==1)dpx=sqrt(dpx)
     dp2=dpx*dpx

! poloidal spline is regular
     thetatmp=zion(2,m)
     jst=max(1,min(lst-1,ceiling(thetatmp*spdtheta_inv)))
     dtx=thetatmp-spdtheta*real(jst-1)
     dt2=dtx*dtx

     rg= rgpsi(1,isp)  +rgpsi(2,isp)*dpx  +rgpsi(3,isp)*dp2
#ifndef _TOROIDAL3D
     b=    bsp(1,isp,jst)+bsp(2,isp,jst)*dpx+bsp(3,isp,jst)*dp2+ &
          (bsp(4,isp,jst)+bsp(5,isp,jst)*dpx+bsp(6,isp,jst)*dp2)*dtx+ &
          (bsp(7,isp,jst)+bsp(8,isp,jst)*dpx+bsp(9,isp,jst)*dp2)*dt2
#else
     ztmp=zion(3,m)
     b=  spb(psitmp,thetatmp,ztmp)    
#endif

     fullf=zion(7,m)
     delf=fullf*zion(5,m)
     rmu=zion(6,m)
     upara=zion(4,m)*b*qion*ai_inv*vi_inv 
     energy=0.5_lk*upara*upara+rmu*rmu*b*ai_inv*vi_inv*vi_inv ! energy is normalized by etemp0
     pitch=upara/sqrt(2.0_lk*energy)

! particles sorted into bins in radius, energy, and pitch
     ip=max(1,min(mpsi,ceiling((rg-rg0)*delr)))! radial bin
     jenergy=max(1,min(nvgrid,ceiling(real(nvgrid)*energy*emax_inv*rho0*rho0/tipp(1,1))))! energy bin
     kpitch=max(1,min(nvgrid,ceiling(real(nvgrid)*0.5_lk*(pitch+1.0_lk))))! pitch bin
     !lambda=mu*B0/E
     lambdabin=max(1,min(nvgrid,ceiling(real(nvgrid)*rmu*rmu*lambmax_inv/(energy*rho0*rho0))))

     pdum=(psimesh(ip)-psitmp)/deltap(ip)
     profile(ip-1,1,ns)=profile(ip-1,1,ns)+pdum*fullf
     profile(ip,1,ns)=profile(ip,1,ns)+(1.0_lk-pdum)*fullf
     profile(ip-1,2,ns)=profile(ip-1,2,ns)+pdum*delf
     profile(ip,2,ns)=profile(ip,2,ns)+(1.0_lk-pdum)*delf
     profile(ip-1,3,ns)=profile(ip-1,3,ns)+pdum*fullf*upara
     profile(ip,3,ns)=profile(ip,3,ns)+(1.0_lk-pdum)*fullf*upara
     profile(ip-1,4,ns)=profile(ip-1,4,ns)+pdum*delf*upara
     profile(ip,4,ns)=profile(ip,4,ns)+(1.0_lk-pdum)*delf*upara
     profile(ip-1,5,ns)=profile(ip-1,5,ns)+pdum*fullf*energy
     profile(ip,5,ns)=profile(ip,5,ns)+(1.0_lk-pdum)*fullf*energy
     profile(ip-1,6,ns)=profile(ip-1,6,ns)+pdum*delf*energy
     profile(ip,6,ns)=profile(ip,6,ns)+(1.0_lk-pdum)*delf*energy

     pdf(jenergy,1,ns)=pdf(jenergy,1,ns)+fullf
     pdf(jenergy,2,ns)=pdf(jenergy,2,ns)+delf
     pdf(kpitch,3,ns)=pdf(kpitch,3,ns)+fullf
     pdf(kpitch,4,ns)=pdf(kpitch,4,ns)+delf

     pdf2d(lambdabin,jenergy,1,ns)=pdf2d(lambdabin,jenergy,1,ns)+fullf
     pdf2d(lambdabin,jenergy,2,ns)=pdf2d(lambdabin,jenergy,2,ns)+delf*delf
  enddo
  else
    do m=1,mi
      psitmp=zion(1,m)
      isp=max(1,min(lsp-1,ceiling(psitmp*spdpsi_inv)))
      dpx=psitmp-spdpsi*real(isp-1)
    !> radaial spline of rg & b avoids sigularity near axis: y=y1+y2*sqrt(x)+y3*x for psi<spdpsi
      if(isp==1)dpx=sqrt(dpx)
      dp2=dpx*dpx

    !> poloidal spline is regular
      thetatmp=zion(2,m)
      jst=max(1,min(lst-1,ceiling(thetatmp*spdtheta_inv)))
      dtx=thetatmp-spdtheta*real(jst-1)
      dt2=dtx*dtx

      rg= rgpsi(1,isp)  +rgpsi(2,isp)*dpx  +rgpsi(3,isp)*dp2
      g=    gpsi(1,isp)   +gpsi(2,isp)*dpx    +gpsi(3,isp)*dp2
      ri=   cpsi(1,isp)   +cpsi(2,isp)*dpx    +cpsi(3,isp)*dp2
#ifndef _TOROIDAL3D
    !> 2D equilibrium coordinate splines
      dx(1)=1.0_lk
      dx(2)=dpx
      dx(3)=dp2
      dx(4:6)=dx(1:3)*dtx
      dx(7:9)=dx(1:3)*dt2
      b=0.0_lk
      majorr=0.0_lk
      do ii = 1, 9
        b=b +bsp(ii,isp,jst)*dx(ii)
        majorr=majorr +xsp(ii,isp,jst)*dx(ii)
      enddo
#else
    !> 3D equilibrium coordinates splines
    !> We first calculate the spline indexes. The psi and theta indexesare the same
    !> as calculated from other equilibrium splines.
    !> dzeta needs to be calculated again because of sub-section used in 3D equilibrium
      zetatmp=zion(3,m)
      ksz=max(1,min(nzsp_sec,ceiling((zetatmp-zeta0)*spdzeta_inv)))
      dzx = zetatmp - (zeta0 + spdzeta*(ksz-1))
  
    !> Now we evaluate b and majorr using spline coefficients
      dx(1)=1.0_lk
      dx(2)=dpx
      dx(3)=dp2
      dx(4:6)=dx(1:3)*dtx
      dx(7:9)=dx(1:3)*dt2
      dx(10:18)=dx(1:9)*dzx
      dx(19:27)=dx(10:18)*dzx
      b = 0.0_lk
      majorr = 0.0_lk
      do ii=1,27
        b = b + bsp(ii,isp,jst,ksz)*dx(ii)
        majorr = majorr + xsp(ii,isp,jst,ksz)*dx(ii)
      enddo
#endif
    !yang, 06/11/2020
    !> particle velocity
      q = spq(psitmp)
      dxdptmp=dxdp(psitmp,thetatmp)
      dxdttmp=dxdt(psitmp,thetatmp)
      dzdptmp=dzdp(psitmp,thetatmp)
      dzdttmp=dzdt(psitmp,thetatmp)
      gpp=dxdptmp**2+dzdptmp**2
      gtt=dxdttmp**2+dzdttmp**2
      gzz=majorr**2
      gpt=dxdptmp*dxdttmp+dzdttmp*dzdptmp
      tmp=sqrt(gtt+q**2*gzz) ! B_0 Jacobian
      tmp1=(dxdttmp**2+dzdttmp**2)/(dxdptmp*dzdttmp-dxdttmp*dzdptmp)**2   ! gupp
      b=tmp/majorr/abs(dxdptmp*dzdttmp-dxdttmp*dzdptmp)
 
      upara=(gpt*zion(4,m)+gtt*zion(6,m)+q*gzz*zion(8,m))/tmp
    !  energy=0.5*aion*(zion(4,m)**2/tmp1+(zion(4,m)*q*gpt/gtt+q*zion(6,m)-zion(8,m))**2*tmp1/b**2+upara**2)/rho0/rho0
      energy=0.5*aion*(zion(4,m)**2*gpp+zion(6,m)**2*gtt+zion(8,m)**2*gzz+2.0_lk*zion(4,m)*zion(6,m)*gpt)/rho0/rho0
      upara=upara*vi_inv

      fullf=zion(7,m)
      delf=fullf*zion(5,m)
      pitch=upara/sqrt(2.0_lk*energy)

    !> particles sorted into bins in radius, energy, and pitch
      ip=max(1,min(mpsi,ceiling((rg-rg0)*delr)))! radial bin
      jenergy=max(1,min(nvgrid,ceiling(real(nvgrid)*energy*emax_inv*rho0*rho0/tipp(1,1))))! energy bin
      kpitch=max(1,min(nvgrid,ceiling(real(nvgrid)*0.5*(pitch+1.0))))! pitch bin

      pdum=(psimesh(ip)-psitmp)/deltap(ip)
      profile(ip-1,1,ns)=profile(ip-1,1,ns)+pdum*fullf
      profile(ip,1,ns)=profile(ip,1,ns)+(1.0_lk-pdum)*fullf
      profile(ip-1,2,ns)=profile(ip-1,2,ns)+pdum*delf
      profile(ip,2,ns)=profile(ip,2,ns)+(1.0_lk-pdum)*delf
      profile(ip-1,3,ns)=profile(ip-1,3,ns)+pdum*fullf*upara
      profile(ip,3,ns)=profile(ip,3,ns)+(1.0_lk-pdum)*fullf*upara
      profile(ip-1,4,ns)=profile(ip-1,4,ns)+pdum*delf*upara
      profile(ip,4,ns)=profile(ip,4,ns)+(1.0_lk-pdum)*delf*upara
      profile(ip-1,5,ns)=profile(ip-1,5,ns)+pdum*fullf*energy
      profile(ip,5,ns)=profile(ip,5,ns)+(1.0_lk-pdum)*fullf*energy
      profile(ip-1,6,ns)=profile(ip-1,6,ns)+pdum*delf*energy
      profile(ip,6,ns)=profile(ip,6,ns)+(1.0_lk-pdum)*delf*energy

      pdf(jenergy,1,ns)=pdf(jenergy,1,ns)+fullf
      pdf(jenergy,2,ns)=pdf(jenergy,2,ns)+delf
      pdf(kpitch,3,ns)=pdf(kpitch,3,ns)+fullf
      pdf(kpitch,4,ns)=pdf(kpitch,4,ns)+delf
    enddo
  endif
! species ns=2: electron
  if(nhybrid>0)then
  ns=2
  ae_inv=1.0/aelectron
  ve_inv=sqrt(aelectron)/rho0 ! electron velocity normalized by v_the
  do m=1,me

     psitmp=zelectron(1,m)
     isp=max(1,min(lsp-1,ceiling(psitmp*spdpsi_inv)))
     dpx=psitmp-spdpsi*real(isp-1)
! radaial spline of rg & b avoids sigularity near axis: y=y1+y2*sqrt(x)+y3*x for psi<spdpsi
     if(isp==1)dpx=sqrt(dpx)
     dp2=dpx*dpx

! poloidal spline is regular
     thetatmp=zelectron(2,m)
     jst=max(1,min(lst-1,ceiling(thetatmp*spdtheta_inv)))
     dtx=thetatmp-spdtheta*real(jst-1)
     dt2=dtx*dtx

     rg= rgpsi(1,isp)  +rgpsi(2,isp)*dpx  +rgpsi(3,isp)*dp2
#ifndef _TOROIDAL3D
     b=    bsp(1,isp,jst)+bsp(2,isp,jst)*dpx+bsp(3,isp,jst)*dp2+ &
          (bsp(4,isp,jst)+bsp(5,isp,jst)*dpx+bsp(6,isp,jst)*dp2)*dtx+ &
          (bsp(7,isp,jst)+bsp(8,isp,jst)*dpx+bsp(9,isp,jst)*dp2)*dt2
#else
     ztmp=zelectron(3,m)
     b=  spb(psitmp,thetatmp,ztmp)
#endif

     fullf=zelectron(7,m)
     delf=fullf*zelectron(5,m)
     rmu=zelectron(6,m)
     upara=zelectron(4,m)*b*qelectron*ae_inv*ve_inv 
     energy=0.5_lk*upara*upara+rmu*rmu*b*ae_inv*ve_inv*ve_inv ! energy is normalized by etemp0
     pitch=upara/sqrt(2.0_lk*energy)

     ip=max(1,min(mpsi,ceiling((rg-rg0)*delr)))! radial bin
     jenergy=max(1,min(nvgrid,ceiling(real(nvgrid)*energy*emax_inv))) !energy bin
     kpitch=max(1,min(nvgrid,ceiling(real(nvgrid)*0.5_lk*(pitch+1.0_lk))))  !pitch bin
     !lambda=mu*B0/E
     lambdabin=max(1,min(nvgrid,ceiling(real(nvgrid)*rmu*rmu*lambmax_inv/(energy*rho0*rho0))))

     pdum=(psimesh(ip)-psitmp)/deltap(ip)
     profile(ip-1,1,ns)=profile(ip-1,1,ns)+pdum*fullf
     profile(ip,1,ns)=profile(ip,1,ns)+(1.0_lk-pdum)*fullf
     profile(ip-1,2,ns)=profile(ip-1,2,ns)+pdum*delf
     profile(ip,2,ns)=profile(ip,2,ns)+(1.0_lk-pdum)*delf
     profile(ip-1,3,ns)=profile(ip-1,3,ns)+pdum*fullf*upara
     profile(ip,3,ns)=profile(ip,3,ns)+(1.0_lk-pdum)*fullf*upara
     profile(ip-1,4,ns)=profile(ip-1,4,ns)+pdum*delf*upara
     profile(ip,4,ns)=profile(ip,4,ns)+(1.0_lk-pdum)*delf*upara
     profile(ip-1,5,ns)=profile(ip-1,5,ns)+pdum*fullf*energy
     profile(ip,5,ns)=profile(ip,5,ns)+(1.0_lk-pdum)*fullf*energy
     profile(ip-1,6,ns)=profile(ip-1,6,ns)+pdum*delf*energy
     profile(ip,6,ns)=profile(ip,6,ns)+(1.0_lk-pdum)*delf*energy

     pdf(jenergy,1,ns)=pdf(jenergy,1,ns)+fullf
     pdf(jenergy,2,ns)=pdf(jenergy,2,ns)+delf
     pdf(kpitch,3,ns)=pdf(kpitch,3,ns)+fullf
     pdf(kpitch,4,ns)=pdf(kpitch,4,ns)+delf

     pdf2d(lambdabin,jenergy,1,ns)=pdf2d(lambdabin,jenergy,1,ns)+fullf
     pdf2d(lambdabin,jenergy,2,ns)=pdf2d(lambdabin,jenergy,2,ns)+delf*delf
  enddo
  endif

! species ns=3: fast ion
  if(fullkf==1)then
    ns=ns+1
    af_inv=1.0/afast
    vf_inv=sqrt(afast)/rho0
    do m=1,mf
      psitmp=zfast(1,m)
      isp=max(1,min(lsp-1,ceiling(psitmp*spdpsi_inv)))
      dpx=psitmp-spdpsi*real(isp-1)
    !> radaial spline of rg & b avoids sigularity near axis: y=y1+y2*sqrt(x)+y3*x for psi<spdpsi
      if(isp==1)dpx=sqrt(dpx)
      dp2=dpx*dpx

    !> poloidal spline is regular
      thetatmp=zfast(2,m)
      jst=max(1,min(lst-1,ceiling(thetatmp*spdtheta_inv)))
      dtx=thetatmp-spdtheta*real(jst-1)
      dt2=dtx*dtx

      rg= rgpsi(1,isp)  +rgpsi(2,isp)*dpx  +rgpsi(3,isp)*dp2
      g=    gpsi(1,isp)   +gpsi(2,isp)*dpx    +gpsi(3,isp)*dp2
      ri=   cpsi(1,isp)   +cpsi(2,isp)*dpx    +cpsi(3,isp)*dp2
#ifndef _TOROIDAL3D
    !> 2D equilibrium coordinate splines
      dx(1)=1.0_lk
      dx(2)=dpx
      dx(3)=dp2
      dx(4:6)=dx(1:3)*dtx
      dx(7:9)=dx(1:3)*dt2
      b=0.0_lk
      majorr=0.0_lk
      do ii = 1, 9
        b=b +bsp(ii,isp,jst)*dx(ii)
        majorr=majorr +xsp(ii,isp,jst)*dx(ii)
      enddo
#else
    !> 3D equilibrium coordinates splines
    !> We first calculate the spline indexes. The psi and theta indexesare the same
    !> as calculated from other equilibrium splines.
    !> dzeta needs to be calculated again because of sub-section used in 3D equilibrium
      zetatmp=zfast(3,m)
      ksz=max(1,min(nzsp_sec,ceiling((zetatmp-zeta0)*spdzeta_inv)))
      dzx = zetatmp - (zeta0 + spdzeta*(ksz-1))
  
    !> Now we evaluate b and majorr using spline coefficients
      dx(1)=1.0_lk
      dx(2)=dpx
      dx(3)=dp2
      dx(4:6)=dx(1:3)*dtx
      dx(7:9)=dx(1:3)*dt2
      dx(10:18)=dx(1:9)*dzx
      dx(19:27)=dx(10:18)*dzx
      b = 0.0_lk
      majorr = 0.0_lk
      do ii=1,27
        b = b + bsp(ii,isp,jst,ksz)*dx(ii)
        majorr = majorr + xsp(ii,isp,jst,ksz)*dx(ii)
      enddo
#endif
    !yang, 06/11/2020
    !> particle velocity
      q = spq(psitmp)
      dxdptmp=dxdp(psitmp,thetatmp)
      dxdttmp=dxdt(psitmp,thetatmp)
      dzdptmp=dzdp(psitmp,thetatmp)
      dzdttmp=dzdt(psitmp,thetatmp)
      gpp=dxdptmp**2+dzdptmp**2
      gtt=dxdttmp**2+dzdttmp**2
      gzz=majorr**2
      gpt=dxdptmp*dxdttmp+dzdttmp*dzdptmp
      tmp=sqrt(gtt+q**2*gzz) ! B_0 Jacobian
      tmp1=(dxdttmp**2+dzdttmp**2)/(dxdptmp*dzdttmp-dxdttmp*dzdptmp)**2   ! gupp
      b=tmp/majorr/abs(dxdptmp*dzdttmp-dxdttmp*dzdptmp)
 
      upara=(gpt*zfast(4,m)+gtt*zfast(6,m)+q*gzz*zfast(8,m))/tmp
    !  energy=0.5*afast*(zfast(4,m)**2/tmp1+(zfast(4,m)*q*gpt/gtt+q*zfast(6,m)-zfast(8,m))**2*tmp1/b**2+upara**2)/rho0/rho0
      energy=0.5*afast*(zfast(4,m)**2*gpp+zfast(6,m)**2*gtt+zfast(8,m)**2*gzz+2.0_lk*zfast(4,m)*zfast(6,m)*gpt)/rho0/rho0
      upara=upara*vf_inv

      fullf=zfast(7,m)
      delf=fullf*zfast(5,m)
      pitch=upara/sqrt(2.0_lk*energy)

    !> particles sorted into bins in radius, energy, and pitch
      ip=max(1,min(mpsi,ceiling((rg-rg0)*delr)))! radial bin
      jenergy=max(1,min(nvgrid,ceiling(real(nvgrid)*energy*emax_inv*rho0*rho0/tfpp(1,1))))! energy bin
      kpitch=max(1,min(nvgrid,ceiling(real(nvgrid)*0.5*(pitch+1.0))))! pitch bin

      pdum=(psimesh(ip)-psitmp)/deltap(ip)
      profile(ip-1,1,ns)=profile(ip-1,1,ns)+pdum*fullf
      profile(ip,1,ns)=profile(ip,1,ns)+(1.0_lk-pdum)*fullf
      profile(ip-1,2,ns)=profile(ip-1,2,ns)+pdum*delf
      profile(ip,2,ns)=profile(ip,2,ns)+(1.0_lk-pdum)*delf
      profile(ip-1,3,ns)=profile(ip-1,3,ns)+pdum*fullf*upara
      profile(ip,3,ns)=profile(ip,3,ns)+(1.0_lk-pdum)*fullf*upara
      profile(ip-1,4,ns)=profile(ip-1,4,ns)+pdum*delf*upara
      profile(ip,4,ns)=profile(ip,4,ns)+(1.0_lk-pdum)*delf*upara
      profile(ip-1,5,ns)=profile(ip-1,5,ns)+pdum*fullf*energy
      profile(ip,5,ns)=profile(ip,5,ns)+(1.0_lk-pdum)*fullf*energy
      profile(ip-1,6,ns)=profile(ip-1,6,ns)+pdum*delf*energy
      profile(ip,6,ns)=profile(ip,6,ns)+(1.0_lk-pdum)*delf*energy

      pdf(jenergy,1,ns)=pdf(jenergy,1,ns)+fullf
      pdf(jenergy,2,ns)=pdf(jenergy,2,ns)+delf
      pdf(kpitch,3,ns)=pdf(kpitch,3,ns)+fullf
      pdf(kpitch,4,ns)=pdf(kpitch,4,ns)+delf
    enddo
  elseif(fload>0)then
  ns=ns+1
  ai_inv=1.0/afast
  vi_inv=sqrt(afast)/rho0
  
  do m=1,mf
     psitmp=zfast(1,m)
     isp=max(1,min(lsp-1,ceiling(psitmp*spdpsi_inv)))
     dpx=psitmp-spdpsi*real(isp-1)
! radaial spline of rg & b avoids sigularity near axis: y=y1+y2*sqrt(x)+y3*x for psi<spdpsi
     if(isp==1)dpx=sqrt(dpx)
     dp2=dpx*dpx

! poloidal spline is regular
     thetatmp=zfast(2,m)
     jst=max(1,min(lst-1,ceiling(thetatmp*spdtheta_inv)))
     dtx=thetatmp-spdtheta*real(jst-1)
     dt2=dtx*dtx

     rg= rgpsi(1,isp)  +rgpsi(2,isp)*dpx  +rgpsi(3,isp)*dp2
#ifndef _TOROIDAL3D
     b=    bsp(1,isp,jst)+bsp(2,isp,jst)*dpx+bsp(3,isp,jst)*dp2+ &
          (bsp(4,isp,jst)+bsp(5,isp,jst)*dpx+bsp(6,isp,jst)*dp2)*dtx+ &
          (bsp(7,isp,jst)+bsp(8,isp,jst)*dpx+bsp(9,isp,jst)*dp2)*dt2
#else
     ztmp=zfast(3,m)
     b=  spb(psitmp,thetatmp,ztmp)
#endif

     fullf=zfast(7,m)
     delf=fullf*zfast(5,m)
     rmu=zfast(6,m)
     upara=zfast(4,m)*b*qfast*ai_inv*vi_inv 
     energy=0.5_lk*upara*upara+rmu*rmu*b*ai_inv*vi_inv*vi_inv
     pitch=upara/sqrt(2.0_lk*energy)

! particles sorted into bins in radius, energy, and pitch
     ip=max(1,min(mpsi,ceiling((rg-rg0)*delr)))! radial bin
     jenergy=max(1,min(nvgrid,ceiling(real(nvgrid)*energy*emax_inv*rho0*rho0/tfpp(1,1))))! energy bin
     kpitch=max(1,min(nvgrid,ceiling(real(nvgrid)*0.5_lk*(pitch+1.0))))! pitch bin
     !lambda=mu*B0/E
     lambdabin=max(1,min(nvgrid,ceiling(real(nvgrid)*rmu*rmu*lambmax_inv/(energy*rho0*rho0))))

     pdum=(psimesh(ip)-psitmp)/deltap(ip)
     profile(ip-1,1,ns)=profile(ip-1,1,ns)+pdum*fullf
     profile(ip,1,ns)=profile(ip,1,ns)+(1.0_lk-pdum)*fullf
     profile(ip-1,2,ns)=profile(ip-1,2,ns)+pdum*delf
     profile(ip,2,ns)=profile(ip,2,ns)+(1.0_lk-pdum)*delf
     profile(ip-1,3,ns)=profile(ip-1,3,ns)+pdum*fullf*upara
     profile(ip,3,ns)=profile(ip,3,ns)+(1.0_lk-pdum)*fullf*upara
     profile(ip-1,4,ns)=profile(ip-1,4,ns)+pdum*delf*upara
     profile(ip,4,ns)=profile(ip,4,ns)+(1.0_lk-pdum)*delf*upara
     profile(ip-1,5,ns)=profile(ip-1,5,ns)+pdum*fullf*energy
     profile(ip,5,ns)=profile(ip,5,ns)+(1.0_lk-pdum)*fullf*energy
     profile(ip-1,6,ns)=profile(ip-1,6,ns)+pdum*delf*energy
     profile(ip,6,ns)=profile(ip,6,ns)+(1.0_lk-pdum)*delf*energy

     pdf(jenergy,1,ns)=pdf(jenergy,1,ns)+fullf
     pdf(jenergy,2,ns)=pdf(jenergy,2,ns)+delf
     pdf(kpitch,3,ns)=pdf(kpitch,3,ns)+fullf
     pdf(kpitch,4,ns)=pdf(kpitch,4,ns)+delf

     pdf2d(lambdabin,jenergy,1,ns)=pdf2d(lambdabin,jenergy,1,ns)+fullf
     pdf2d(lambdabin,jenergy,2,ns)=pdf2d(lambdabin,jenergy,2,ns)+delf*delf
  enddo
  endif
  
! species ns=4: fast electron
  if(feload>0)then
  ns=ns+1
  afe_inv=1.0_lk/afaste
  vfe_inv=sqrt(afaste)/rho0

  do m=1,mfe
     psitmp=zfaste(1,m)
     isp=max(1,min(lsp-1,ceiling(psitmp*spdpsi_inv)))
     dpx=psitmp-spdpsi*real(isp-1)
! radaial spline of rg & b avoids sigularity near axis: y=y1+y2*sqrt(x)+y3*x for psi<spdpsi
     if(isp==1)dpx=sqrt(dpx)
     dp2=dpx*dpx

! poloidal spline is regular
     thetatmp=zfaste(2,m)
     jst=max(1,min(lst-1,ceiling(thetatmp*spdtheta_inv)))
     dtx=thetatmp-spdtheta*real(jst-1)
     dt2=dtx*dtx

     rg= rgpsi(1,isp)  +rgpsi(2,isp)*dpx  +rgpsi(3,isp)*dp2
#ifndef _TOROIDAL3D
     b=    bsp(1,isp,jst)+bsp(2,isp,jst)*dpx+bsp(3,isp,jst)*dp2+ &
          (bsp(4,isp,jst)+bsp(5,isp,jst)*dpx+bsp(6,isp,jst)*dp2)*dtx+ &
          (bsp(7,isp,jst)+bsp(8,isp,jst)*dpx+bsp(9,isp,jst)*dp2)*dt2
#else
     ztmp=zfaste(3,m)
     b=  spb(psitmp,thetatmp,ztmp)
#endif

     fullf=zfaste(7,m)
     delf=fullf*zfaste(5,m)
     rmu=zfaste(6,m)
     upara=zfaste(4,m)*b*qfaste*afe_inv*vfe_inv
     energy=0.5_lk*upara*upara+rmu*rmu*b*afe_inv*vfe_inv*vfe_inv
     pitch=upara/sqrt(2.0_lk*energy)

! particles sorted into bins in radius, energy, and pitch
     ip=max(1,min(mpsi,ceiling((rg-rg0)*delr)))! radial bin
     jenergy=max(1,min(nvgrid,ceiling(real(nvgrid)*energy*emax_inv*rho0*rho0/tfepp(1,1))))! energy bin
     kpitch=max(1,min(nvgrid,ceiling(real(nvgrid)*0.5*(pitch+1.0))))! pitch bin
     !lambda=mu*B0/E
     lambdabin=max(1,min(nvgrid,ceiling(real(nvgrid)*rmu*rmu*lambmax_inv/(energy*rho0*rho0))))

     pdum=(psimesh(ip)-psitmp)/deltap(ip)
     profile(ip-1,1,ns)=profile(ip-1,1,ns)+pdum*fullf
     profile(ip,1,ns)=profile(ip,1,ns)+(1.0-pdum)*fullf
     profile(ip-1,2,ns)=profile(ip-1,2,ns)+pdum*delf
     profile(ip,2,ns)=profile(ip,2,ns)+(1.0-pdum)*delf
     profile(ip-1,3,ns)=profile(ip-1,3,ns)+pdum*fullf*upara
     profile(ip,3,ns)=profile(ip,3,ns)+(1.0-pdum)*fullf*upara
     profile(ip-1,4,ns)=profile(ip-1,4,ns)+pdum*delf*upara
     profile(ip,4,ns)=profile(ip,4,ns)+(1.0-pdum)*delf*upara
     profile(ip-1,5,ns)=profile(ip-1,5,ns)+pdum*fullf*energy
     profile(ip,5,ns)=profile(ip,5,ns)+(1.0-pdum)*fullf*energy
     profile(ip-1,6,ns)=profile(ip-1,6,ns)+pdum*delf*energy
     profile(ip,6,ns)=profile(ip,6,ns)+(1.0-pdum)*delf*energy

     pdf(jenergy,1,ns)=pdf(jenergy,1,ns)+fullf
     pdf(jenergy,2,ns)=pdf(jenergy,2,ns)+delf
     pdf(kpitch,3,ns)=pdf(kpitch,3,ns)+fullf
     pdf(kpitch,4,ns)=pdf(kpitch,4,ns)+delf

     pdf2d(lambdabin,jenergy,1,ns)=pdf2d(lambdabin,jenergy,1,ns)+fullf
     pdf2d(lambdabin,jenergy,2,ns)=pdf2d(lambdabin,jenergy,2,ns)+delf*delf
  enddo
  endif

! sum over MPI tasks
  icount=nspecies*6*(mpsi+1)
  call MPI_REDUCE(profile,proftmp,icount,mpi_Rsize,MPI_SUM,0,MPI_COMM_WORLD,ierror)
  profile=proftmp

  icount=nspecies*4*(nvgrid)
  call MPI_REDUCE(pdf,pdftmp,icount,mpi_Rsize,MPI_SUM,0,MPI_COMM_WORLD,ierror)
  pdf=pdftmp

  icount=nvgrid*nvgrid*2*nspecies
  call MPI_REDUCE(pdf2d,pdf2dtmp,icount,mpi_Rsize,MPI_SUM,0,MPI_COMM_WORLD,ierror)
  pdf2d=pdf2dtmp
! normalization by marker #
  do ns=1,nspecies
     do ip=0,mpsi
        profile(ip,2:6,ns)=profile(ip,2:6,ns)/profile(ip,1,ns)
     enddo
     if(ns==1)then
        profile(:,1,ns)=profile(:,1,ns)/pmarki

     elseif(ns==2)then
        if(nhybrid>0)then
           profile(:,1,ns)=profile(:,1,ns)/pmarke
        else
           if(fload>0)then
              profile(:,1,ns)=profile(:,1,ns)/pmarkf
           else
               if(feload>0)then
                 profile(:,1,ns)=profile(:,1,ns)/pmarkfe
              endif
           endif 
        endif
     elseif(ns==3)then
          if(nhybrid>0)then
           if(fload>0)then
              profile(:,1,ns)=profile(:,1,ns)/pmarkf
           else
              if(feload>0)then
                 profile(:,1,ns)=profile(:,1,ns)/pmarkfe
              endif
           endif
        else
           profile(:,1,ns)=profile(:,1,ns)/pmarkfe
        endif
     elseif(nf==4)then
        profile(:,1,ns)=profile(:,1,ns)/pmarkfe
     endif

     do j=1,nvgrid
        pdf(j,2,ns)=pdf(j,2,ns)/max(1.0,pdf(j,1,ns))
        pdf(j,4,ns)=pdf(j,4,ns)/max(1.0,pdf(j,3,ns))
     enddo
  enddo

! poloidal resolution=poloidal grid on diag_flux surface
  mtgrid=mtheta(diag_flux)

  allocate(poloidata(0:mtgrid,0:mpsi,nfield+2),fluxdata(0:mtgrid,mtoroidal,nfield),&
    eachflux(mtgrid),allflux(mtgrid,mtoroidal),bpararms(0:mpsi),bperprms(0:mpsi),bratio(0:mpsi),bperpdiag(mgrid),bparadiag(mgrid), &
    phieff(mgrid),phieffrms(0:mpsi),phirms(0:mpsi),phiratio(0:mpsi),eachrz(0:mpsi),allrz(0:mpsi,mtoroidal),  STAT=mstat)
  call check_allocation(mstat, "poloidata in snapshot.F90")
  poloidata=0.0_lk
  fluxdata=0.0_lk
  eachrz=0.0_lk
  allrz=0.0_lk

  !$acc update host(sfluidne,sdeltapsi,fluidne,sapara,sapara_na,gradapara,bpara)
  !$acc update host(phi,gradphi,gradphieff)


! field quantities: phi, a_para, fluidne. Last two coloumn of poloidal for coordinates
  dkne=0.0_lk
  phirms=0.0_lk
  phieffrms=0.0_lk
  bpararms=0.0_lk
  bperprms=0.0_lk
  bratio=0.0_lk
  phiratio=0.0_lk
  sbratio=0.0_lk
  sphiratio=0.0_lk
  phieff=0.0_lk
  bparadiag=0.0_lk
  bperpdiag=0.0_lk

  do i=0,mpsi
     pdum=psimesh(i)
     do j=0,mtheta(i)
        ij=igrid(i)+j
        if(fielddir==1 .or. fielddir==3)then
          tdum=modulo(deltat(i)*real(j)+(zeta1-pi2)*qtinv(i),pi2)
        else
          tdum=modulo(deltat(i)*real(j)+zeta1*qtinv(i),pi2)
        endif

        sqdet_inv = 1.0_lk/sqrt(gupp(ij)*(gutt(ij)*guzz(ij)-gutz(ij)*gutz(ij)) - gupt(ij)*(gupt(ij)*guzz(ij)-gupz(ij)*gutz(ij)) + &
                                gupz(ij)*(gupt(ij)*gutz(ij)-gutt(ij)*gutz(ij)))
        gradapara_phi = gradapara(3,1,ij)-gradapara(2,1,ij)/qmesh(i)
        bperpdiag(ij)=sqrt(gradapara(1,1,ij)**2*gupp(ij) + gradapara(2,1,ij)**2*gutt(ij) + gradapara_phi**2*guzz(ij) &
        +2.0_lk*(gradapara(1,1,ij)*gradapara(2,1,ij)*gupt(ij) + gradapara(1,1,ij)*gradapara_phi*gupz(ij) + gradapara(2,1,ij)*gradapara_phi*gutz(ij)))
        bperprms(i)=bperprms(i)+(gradapara(1,1,ij)**2*gupp(ij)+gradapara(2,1,ij)**2*gutt(ij) &
                    +(gradapara(3,1,ij)-gradapara(2,1,ij)/qmesh(i))**2*guzz(ij))/mtheta(i)
        phirms(i)=phirms(i)+phi(1,ij)**2/mtheta(i)
        if(nhybrid>0)then
           dkne=densitye(1,ij)
        else
           dkne=0.0_lk
        endif
        if(deltabpara==0)then
          phieff(ij)=sfluidne(1,ij)*meshte(i)/meshne(i)+sdeltapsi(1,ij)*meshte(i)*kapane(i)-dkne*meshte(i)
          phieffrms(i)=phieffrms(i)+(sfluidne(1,ij)*meshte(i)/meshne(i) &
                       +sdeltapsi(1,ij)*meshte(i)*kapane(i)-dkne*meshte(i))**2/mtheta(i)
        else
          bparadiag(ij)=(bpara(1,ij)*bmesh(ij))
          bpararms(i) = bpararms(i)+(bpara(1,ij)*bmesh(ij))**2/mtheta(i)
          phieff(ij) = sfluidne(1,ij)*meshte(i)/meshne(i)+bpara(1,ij)*meshte(i) &
                       + sdeltapsi(1,ij)*meshte(i)*kapane(i)-dkne*meshte(i)
          phieffrms(i)=phieffrms(i)+(sfluidne(1,ij)*meshte(i)/meshne(i) + bpara(1,ij)*meshte(i) &
                       +sdeltapsi(1,ij)*meshte(i)*kapane(i)-dkne*meshte(i))**2/mtheta(i)
        endif
     enddo
     bratio(i)=sqrt(bpararms(i)/bperprms(i))
     phiratio(i)=sqrt(phieffrms(i)/phirms(i))
  enddo

  if(deltabpara==1)then
    sbratio=sqrt(sum(bpararms)/sum(bperprms))
  endif
  sphiratio=sqrt(sum(phieffrms)/sum(phirms))

  do nf=1,nfield
     if(nf==1)then
        dfield=phi(0,:)!bparadiag(:)
     endif
     if(nf==2)then
#ifdef _FRC
             dfield=densityi(1,:)
#else
             dfield=sapara(0,:)!bperpdiag(:)
#endif
     endif
     if(nf==3)then
#ifdef _FRC
        if(feload>0)then
            dfield=densityfe(1,:)
        elseif(nhybrid>0)then
            dfield=densitye(1,:)
        endif
#else   
        if(cs_method==1)then    
            dfield=sapara_na(0,:)
        else
            dfield=sfluidne(0,:)!xi_MHD(:)*ndiag/istep!
        endif
#endif
    endif 

! gather potential on diag_flux surface
     allflux=0.0
     do j=1,mtgrid
        eachflux(j)=dfield(igrid(diag_flux)+j)
     enddo
!fujy gather potential on r-theta surface at theta=0
     do j=0,mpsi
        eachrz(j)=dfield(igrid(j)+1)
     enddo
     icount=mpsi+1
     call MPI_GATHER(eachrz,icount,mpi_Rsize,allrz,icount,mpi_Rsize,0,toroidal_comm,ierror) 
!fujy allrz is the output rz phi    

if(mype==0.and.istep/=0.and.nf==1)then
!   open(41,file='rzdata.out',status='unknown',position='append')
   open(41,file='rzdata.out',status='replace')
     write(41,35)allrz
   close(41)
endif
35 format(122e15.7)


     icount=mtgrid
     call MPI_GATHER(eachflux,icount,mpi_Rsize,allflux,icount,mpi_Rsize,0,toroidal_comm,ierror) 
     fluxdata(1:mtgrid,:,nf)=allflux

! poloidal BC
     fluxdata(0,:,nf)=fluxdata(mtgrid,:,nf)

! potential data on poloidal plain uses polar coordinates
     do j=0,mtgrid
        tdum=pi2*real(j)/real(mtgrid)
        do i=0,mpsi
           if(fielddir==1 .or. fielddir==3)then
             jt=max(0,min(mtheta(i),1+int(modulo(tdum+pi2*qtinv(i),pi2)/deltat(i))))
             wt=modulo(tdum+pi2*qtinv(i),pi2)/deltat(i)-real(jt-1)
           else
             jt=max(0,min(mtheta(i),1+int(tdum/deltat(i))))
             wt=tdum/deltat(i)-real(jt-1)
           endif
           poloidata(j,i,nf)=(wt*dfield(igrid(i)+jt)+(1.0-wt)*dfield(igrid(i)+jt-1)) !/rho0**2
           if(nf==1)then
              poloidata(j,i,nf)=poloidata(j,i,nf)/(rho0*rho0)
           elseif(nf==2)then
! apara is renormalized in such a way that it has the same
! amplitude as phi in ideal shear Alfven waves
#ifndef _FRC
              poloidata(j,i,nf)=poloidata(j,i,nf)/(rho0*sqrt(betae*aion))
#endif
! no need to renormalize fluidne
! renormalize non-adiabatic delta_apara in conservative scheme
           elseif(nf==3 .and. cs_method==1)then
              poloidata(j,i,nf)=poloidata(j,i,nf)/(rho0*sqrt(betae*aion))
           endif
        enddo
     enddo
  enddo

! poloidal grid position in polar coordinates
  do j=0,mtgrid
     tdum=2.0_lk*pi*real(j)/real(mtgrid)
     do i=0,mpsi
        pdum=psimesh(i)
#ifdef _TOROIDAL3D
        poloidata(j,i,nfield+1)=spx(pdum,tdum,0.0_lk)
        poloidata(j,i,nfield+2)=spz(pdum,tdum,0.0_lk)
#else
        poloidata(j,i,nfield+1)=spx(pdum,tdum)
        poloidata(j,i,nfield+2)=spz(pdum,tdum)
#endif
      enddo
  enddo


  niflux=0.0_lk
  do i=0,mpsi
     do j=0,mtheta(i)
        ij=igrid(i)+j
        niflux(i)=niflux(i)+densityi(0,ij)
     enddo
     niflux(i)=niflux(i)/(mtheta(i)+1)
  enddo

! open snapshot output file
  if(mype==0)then
    nsnap=mstepall+istep
#ifdef ADIOS2
    if (.not. adios2_snapshot_io_initialized) then
      call adios2_declare_io (adios2_snapshot_io, adios2_obj, adios2_snapshot_io_id, ierr)
      call error_check (ierr, "adios2 snapshot declare_io")

      ! define scalars
      call adios2_define_variable(var_step_id, adios2_snapshot_io, "step_id", adios2_type_integer4, 1, &
                                  1_8 * (/ 1 /), 1_8 * (/ 0 /), 1_8 * (/ 1 /), .true., ierr)
      call adios2_define_variable(var_nspecies, adios2_snapshot_io, "nspecies", adios2_type_integer4, ierr)
      call adios2_define_variable(var_nfield, adios2_snapshot_io, "nfield", adios2_type_integer4, ierr)
      call adios2_define_variable(var_nvgrid, adios2_snapshot_io, "nvgrid", adios2_type_integer4, ierr)
      call adios2_define_variable(var_mpsi, adios2_snapshot_io, "mpsi", adios2_type_integer4, ierr)
      call adios2_define_variable(var_mtgrid, adios2_snapshot_io, "mtgrid", adios2_type_integer4, ierr)
      call adios2_define_variable(var_mtoroidal, adios2_snapshot_io, "mtoroidal", adios2_type_integer4, ierr)
      call adios2_define_variable(var_emax_inv, adios2_snapshot_io, "tmax", adios2_type_real_lk, ierr)

      ! define arrays
      shape3 = shape(profile)
      call adios2_define_variable(var_profile, adios2_snapshot_io, "profile", adios2_type_real_lk, 3, &
                                  1_8 * shape3, 1_8 * (/ 0,0,0 /), 1_8 * shape3, .false., ierr)

      shape3 = shape(pdf)
      call adios2_define_variable(var_pdf, adios2_snapshot_io, "pdf", adios2_type_real_lk, 3, &
                                  1_8 * shape3, 1_8 * (/ 0,0,0 /), 1_8 * shape3, .false., ierr)

      shape3 = shape(poloidata)
      call adios2_define_variable(var_poloidata, adios2_snapshot_io, "poloidata", adios2_type_real_lk, 3, &
                                  1_8 * shape3, 1_8 * (/ 0,0,0 /), 1_8 * shape3, .false., ierr)

      shape3 = shape(fluxdata)
      call adios2_define_variable(var_fluxdata, adios2_snapshot_io, "fluxdata", adios2_type_real_lk, 3, &
                                  1_8 * shape3, 1_8 * (/ 0,0,0 /), 1_8 * shape3, .false., ierr)

      call adios2_define_variable(var_density, adios2_snapshot_io, "density", adios2_type_real_lk, 1, &
                                  1_8 * (/ 5*mpsi+5 /), 1_8 * (/ 0 /), 1_8 * (/ 5*mpsi+5 /), .false., ierr)
      
      call adios2_open(snapshot_engine, adios2_snapshot_io, "snapshot.bp"//char(0), adios2_mode_append, mpi_comm_self, ierr)
      call error_check(ierr, "adios2_open in snapshot")

      adios2_snapshot_io_initialized = .true.
    endif

    ! The 3D arrays have changing dimensions. Need to set their shape and selection again
    call adios2_set_shape(var_profile, 3, 1_8 * shape(profile), ierr)
    call adios2_set_selection(var_profile, 3, 1_8 * (/ 0,0,0 /), 1_8 * shape(profile), ierr)

    call adios2_set_shape(var_pdf, 3, 1_8 * shape(pdf), ierr)
    call adios2_set_selection(var_pdf, 3, 1_8 * (/ 0,0,0 /), 1_8 * shape(pdf), ierr)

    call adios2_set_shape(var_poloidata, 3, 1_8 * shape(poloidata), ierr)
    call adios2_set_selection(var_poloidata, 3, 1_8 * (/ 0,0,0 /), 1_8 * shape(poloidata), ierr)

    call adios2_set_shape(var_fluxdata, 3, 1_8 * shape(fluxdata), ierr)
    call adios2_set_selection(var_fluxdata, 3, 1_8 * (/ 0,0,0 /), 1_8 * shape(fluxdata), ierr)

    ! Can't pass 'mpsi+1', 'mtgrid+1', '1.0/emax_inv' to the adios2_put call.
    ! It doesn't get the correct values. Maybe something to do with the bindings.
    ! Put them in a tmp variable and pass them on.
    mpsi_tmp = mpsi+1
    mtgrid_tmp = mtgrid+1
    emax_inv_tmp = 1.0/emax_inv
    
    ! Write density output
    allocate(density_adios(0:5*mpsi+5), stat=mstat)
    call check_allocation(mstat, "density_adios in snapshot.F90")

    density_adios = 0.0

    indx = 0
    do i=0, mpsi
       ij=igrid(i)
       if (nhybrid<1) then
          density_adios(indx) = (densityi(1,ij+mtheta(i)/2)+1.0)*eqmeshni(i)
          indx = indx+1
          density_adios(indx) = eqmeshni(i)*(densityi(1,ij)+1.0)
          indx = indx+1
       else
          density_adios(indx) = (densityi(1,ij+mtheta(i)/2)+1.0)*eqmeshni(i)
          indx = indx+1
          density_adios(indx) = eqmeshni(i)*(densityi(1,ij)+1.0)
          indx = indx+1
          density_adios(indx) = (densitye(1,ij)+1)*eqmeshni(i)
          indx = indx+1
          density_adios(indx) = (densitye(1,ij+mtheta(i)/2)+1)*eqmeshni(i)
          indx = indx+1
          density_adios(indx) = eqmeshni(i)*(1.0+pressureepara(1,ij))
          indx = indx+1
       endif
    enddo

    ! Start step and write data
    call adios2_begin_step(snapshot_engine, adios2_step_mode_append, ierr)

    call adios2_put(snapshot_engine, var_step_id, nsnap, ierr)
    call adios2_put(snapshot_engine, var_nspecies, nspecies, ierr)
    call adios2_put(snapshot_engine, var_nfield, nfield, ierr)
    call adios2_put(snapshot_engine, var_nvgrid, nvgrid, ierr)
    call adios2_put(snapshot_engine, var_mpsi, mpsi_tmp, ierr)
    call adios2_put(snapshot_engine, var_mtgrid, mtgrid_tmp, ierr)
    call adios2_put(snapshot_engine, var_mtoroidal, mtoroidal, ierr)
    call adios2_put(snapshot_engine, var_emax_inv, emax_inv_tmp, ierr)
    call adios2_put(snapshot_engine, var_profile, profile, ierr)
    call adios2_put(snapshot_engine, var_pdf, pdf, ierr)
    call adios2_put(snapshot_engine, var_poloidata, poloidata, ierr)
    call adios2_put(snapshot_engine, var_fluxdata, fluxdata, ierr)
    call adios2_put(snapshot_engine, var_density, density_adios, ierr)

    ! End step
    call adios2_end_step(snapshot_engine, ierr)
    deallocate(density_adios)

#endif
! #else
     write(cdum0,'(i7.7,".out")')nsnap
     cdum1='snap'//trim(cdum0)
     open(iosnap,file=cdum1,status='replace')

! parameters: # of species, fields, and grids in velocity, radius, poloidal, toroidal; T_up
     write(iosnap,101)nspecies,nfield,nvgrid,mpsi+1,mtgrid+1,mtoroidal
     write(iosnap,102)1.0/emax_inv

! write out particle radial profile and pdf, and 2D field
     write(iosnap,102)profile,pdf,poloidata,fluxdata
!javierhn: write the phase space of lambda vs Energy. The values in pdf2d(:,:,2,:) correponds to delta_f**2 distribution
     write(iosnap,102) pdf2d

! close snapshot file
     close(iosnap)

     if(magnetic==1)then
        cdum3='iob'//trim(cdum0)
        open(iob,file=cdum3,status='replace')
        write(iob,*)'polarization !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(iob,102)phiratio,sphiratio
        if(deltabpara==1)then
           write(iob,*)'bpara !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
           write(iob,102)sqrt(bpararms),sqrt(bperprms),bratio,sbratio
        endif
        close(iob)
     endif
     
     if(island==1)then
     cdum2='density'//trim(cdum0)
     open(6570,file=cdum2,status='replace') 
     do i=0, mpsi
        ij=igrid(i)
        !write(6570,*) i,eqmeshni(i),densityi(0,ij+mtheta(i)/2-3),densityi(0,ij+mtheta(i)/2-2),densityi(0,ij+mtheta(i)/2-1),densityi(0,ij+mtheta      (i)/2),densityi(0,ij+mtheta(i)/2+1),densityi(0,ij+mtheta(i)/2+2),densityi(0,ij+mtheta(i)/2+3)
!        write(6570,*) i,eqmeshni(i),eqmeshni(i)+densityi(0,ij+mtheta(i)/2),eqmeshni(i)+densityi(0,ij),niflux(i)

!> @bug
!>The 'write' here is a bug for poincare situation,because we don't load electron,
!> and the output of densitye would kill the process of running
        if (nhybrid<1) then
           write(6570,*) i,(densityi(1,ij+mtheta(i)/2)+1.0)*eqmeshni(i),eqmeshni(i)*(densityi(1,ij)+1.0)
        else
           write(6570,*) i,(densityi(1,ij+mtheta(i)/2)+1.0)*eqmeshni(i),eqmeshni(i)*(densityi(1,ij)+1.0),(densitye(1,ij)+1)*eqmeshni(i),(densitye(1,ij+mtheta(i)/2)+1)*eqmeshni(i),eqmeshni(i)*(1.0+pressureepara(1,ij))
        endif
     enddo
     close(6570)
     endif

! End ifdef for ADIOS2 
! #endif
  endif ! if(mype==0)

#ifdef _SDP
! Write out files for Synthetic Diagnostics Platform
  if(mype==0) then
    call sdp_perturbation_output
  endif
#endif

101 format(i6)
102 format(e13.6)
  
end subroutine snapshot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine phishot

use global_parameters
use field_array
use particle_array, only:zonali
#ifdef ADIOS2
use adios2_common
#endif
!use simple_error_handler
implicit none

  integer,parameter :: iophi=223
  integer i,nsnap,ierr
  character(len=64) cdum
  real(lk),dimension(:),allocatable :: phi_out, phi_phi00

#ifdef ADIOS2
  type(adios2_variable), save :: var_mpsi, var_phi, var_phi00, var_zonali, var_phi_phi00

  if (mype == 0) then
    allocate (phi_out(0:mpsi), stat=ierr)
    call error_check (ierr, "allocate in phishot for phi_out")
    allocate (phi_phi00(0:mpsi), stat=ierr)
    call error_check (ierr, "allocate in phishot for phi_phi00")

    do i=0,mpsi
      phi_out(i) = phi(1,igrid(i))
      phi_phi00(i) = phi(1,igrid(i))+phi00(i)
    enddo

    if (.not. adios2_phi_io_initialized) then

      ! Declare io and define variables
      call adios2_declare_io (adios2_phi_io, adios2_obj, adios2_phi_io_id, ierr)

      call adios2_define_variable (var_mpsi, adios2_phi_io, "mpsi", adios2_type_integer4, ierr)

      call adios2_define_variable (var_phi, adios2_phi_io, "phi", adios2_type_real_lk, 1, &
        1_8 * shape(phi_out), (/ 0_8 /), 1_8 * shape(phi_out), .true., ierr)

      call adios2_define_variable (var_phi00, adios2_phi_io, "phi00", adios2_type_real_lk, 1, &
        1_8 * shape(phi00), (/ 0_8 /), 1_8 * shape(phi00), .true., ierr)

      call adios2_define_variable (var_zonali, adios2_phi_io, "zonali", adios2_type_real_lk, 1, &
        1_8 * shape(zonali), (/ 0_8 /), 1_8 * shape(zonali), .true., ierr)

      call adios2_define_variable (var_phi_phi00, adios2_phi_io, "phi_plus_phi00", adios2_type_real_lk, 1, &
        1_8 * shape(phi_phi00), (/ 0_8 /), 1_8 * shape(phi_phi00), .true., ierr)

      call adios2_open (adios2_phi_e, adios2_phi_io, "phi.bp", adios2_mode_write, mpi_comm_self, ierr)

      adios2_phi_io_initialized = .true.
    endif

    ! Write step
    call adios2_begin_step (adios2_phi_e, adios2_step_mode_append, ierr)
    call adios2_put(adios2_phi_e, var_mpsi, mpsi, ierr)
    call adios2_put(adios2_phi_e, var_phi, phi_out, ierr)
    call adios2_put(adios2_phi_e, var_phi00, phi00, ierr)
    call adios2_put(adios2_phi_e, var_zonali, zonali, ierr)
    call adios2_put(adios2_phi_e, var_phi_phi00, phi_phi00, ierr)
    call adios2_end_step (adios2_phi_e, ierr)

    deallocate (phi_out)
    deallocate (phi_phi00)
  endif
#endif
! #else
  if (mype==0) then
      nsnap=mstepall+istep
      write(cdum,'(i7.7,".out")')nsnap
      cdum='phi'//trim(cdum)
      open(iophi,file=cdum,status='replace') 
      do i=0,mpsi
        write(iophi,*)i,phi(1,igrid(i)),phi00(i),zonali(i),phi(1,igrid(i))+phi00(i)
      enddo
      close(iophi)
   endif
! #endif
end subroutine phishot

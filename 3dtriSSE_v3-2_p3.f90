!------------------------------------------------
! 3D Subduction fault earthquake sequence model
! MPI
! Yajing Liu
! Last modified Jul. 4, 2012
!
! D. Li
! Last modified Jan. 2019
!------------------------------------------------
! using varying plate convergence rate. pay attention to theta !!!!
!  yt(2*i) should not be negative!!!!
!!!  with two buffer zones both sides of 20 km wide each
!!  this is a restart code with locked seismogenic zone where slip==0
! cleaned up for higher efficiency

!--------------------------------------------------
! A.Lu 
! 
! modified outputs


program main
  USE mpi
  USE phy3d_module_non
  implicit none
  integer, parameter :: DP=kind(1.d0)

  logical cyclecont
  integer :: Nt_all,Nt,ndt,ndtnext,kk,ii,n,l,ndt_v,ndt_inter,Itout,&
       ndtmax,i,j,k,nm,nrec,ihrestart,Ifileout, &
       Iperb,record,Isnapshot,record_cor,s1,s2,s3,s4,s5,s6,s7,s8

  real (DP) :: accuracy,areatot, epsv,dt_try, dt,dtmin,dt_did,dt_next,dip, &
       hnucl, sigmadiff,sefffix,Lffix,xdepth,xlength,&
       t,tprint_inter, tint_out,tout,tmin_out,tint_cos,&
       tslip_ave,tslip_aveint, tmax, &
       tslipsse,tslipcos,tstart1,tend1,tstart2,tend2,tstart3,tend3,tssestart,tsseend, &
       factor1,factor2,factor3,factor4,vtmp,fr_dummy,&
       xilock1,xilock2,xi_dummy,x_dummy,z_dummy,stiff_dummy,help

  real (DP) ::  tmbegin,tmrun

  real (DP), DIMENSION(:), ALLOCATABLE :: x,z,xi,yt,dydt,yt_scale,tau,tau0, &
       slip,slipinc

  !Arrays only defined at master cpu
  real (DP), DIMENSION(:), ALLOCATABLE :: x_all,xi_all,z_all,&
       yt_all,dydt_all,yt_scale_all,tau_all, &
       slip_all,slipinc_all,cca_all,ccb_all,ccab_all,xLf_all,seff_all

!ED2 !ED5
  !output related parameters   
  !ALu!ED2!ED3 !ED5
  integer :: islip,issere,indtemp,itau,iseff
  integer :: slipout, sseout, cosout, asout,reout,nsouth,nnorth,neqzone,nssezone, nnshallow,nrecenter
  real (DP) :: tcsse
  real (DP), DIMENSION(:), ALLOCATABLE :: tslip,tssere,recenind
  integer, DIMENSION(:), ALLOCATABLE :: indsouth, indnorth, indeq, indsse, indnshallow !ED4
  real (DP), DIMENSION(:,:), ALLOCATABLE :: slip_tau, slip_v, slip_slip,re_slip,re_v,re_tau, re_seff
  integer :: imv,ias,icos,isse,Ioutput,inul,i_nul
  real (DP) :: vcos,vsse1,vsse2
  real (DP), DIMENSION(:), ALLOCATABLE :: maxv,maxv_s1,maxv_s2,maxv_s3,maxv_s4,maxv_s5,&
       maxnum,tmv,tas,tcos,tsse, maxv_s6,maxv_s7,maxv_s8, &
       ttau, mtau, mtau_s1, mtau_s2, mtau_s3, mtau_s4, mtau_s5, mtau_s6, mtau_s7, mtau_s8, &
       mtaunum
!ED5
   real (DP), DIMENSION(:), ALLOCATABLE :: mseff, tseff, mseff_s1, mseff_s2, mseff_s3, &
      mseff_s4, mseff_s5, mseff_s6, mseff_s7, mseff_s8, mseffnum

!ED4 !ED5
   real (DP), DIMENSION(:), ALLOCATABLE :: avev_s, avev_n, maxv_s, maxv_n, maxnum_s, maxnum_n,&
      avetau_s, avetau_n, mtau_s, mtau_n, mtaunum_s, mtaunum_n,&
      vtemp_s, vtemp_n, tautemp_n, tautemp_s,&
      aveslip_s, aveslip_n, sliptemp_s, sliptemp_n,&
      avev_eq, avev_sse, maxv_eq, maxv_sse, maxnum_eq, maxnum_sse,&
      avetau_eq, avetau_sse, mtau_eq, mtau_sse, mtaunum_eq, mtaunum_sse,&
      aveslip_eq, aveslip_sse, sliptemp_eq, sliptemp_sse,&
      vtemp_eq, vtemp_sse, tautemp_eq, tautemp_sse, &
      avev_ns, maxv_ns, maxnum_ns, avetau_ns, mtau_ns, mtaunum_ns, &
      aveslip_ns, sliptemp_ns, vtemp_ns, tautemp_ns, &
      aveseff_s, aveseff_n, aveseff_ns, mseff_s, mseff_n, mseff_ns, &
      sefftemp_s, sefftemp_n, sefftemp_ns, mseffnum_s, mseffnum_n, mseffnum_ns, &
      aveseff_sse, aveseff_eq, mseff_sse, mseff_eq,  &
      sefftemp_sse, sefftemp_eq, mseffnum_sse, mseffnum_eq

   real (DP) :: vsum_s, vsum_n, tausum_s, tausum_n, slipsum_s, slipsum_n,&
               vsum_eq, vsum_sse, tausum_eq, tausum_sse, slipsum_eq, slipsum_sse, &
               vsum_ns, tausum_ns, slipsum_ns, &
               seffsum_s, seffsum_n, seffsum_ns, seffsum_sse, seffsum_eq

  real (DP), DIMENSION(:,:), ALLOCATABLE :: slipz1_inter, slipz1_v,slipz1_cos,slipz1_tau,slipz1_sse

  character(len=40) :: cTemp,filename,ct

  !MPI RELATED DEFINITIONS
  integer :: ierr,size,myid,master
  !input related
   character*128 :: input_file

   call get_command_argument(1,input_file)
  call MPI_Init(ierr)
  CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
  CALL MPI_COMM_SIZE( MPI_COMM_WORLD, size, ierr )

  master = 0

  !read in intiliazation parameters

  open(12,file=input_file,form='formatted',status='old')

  read(12,'(a)')jobname
  read(12,'(a)')foldername
  read(12,'(a)')sstiffname
  read(12,'(a)')nstiffname
  read(12,'(a)')fricname
  read(12,'(a)')patchname  
  read(12,'(a)')restartname
  read(12,'(a)')profile
  read(12,*)Nab,Nt_all,Nt,nprocs,hnucl,nrecenter,nsouth,nnorth,neqzone,nssezone,nnshallow !ED4
  read(12,*)Idin,Idout,Iprofile,Iperb,Isnapshot
  read(12,*)factor1,factor2,factor3,factor4
  read(12,*)Vpl
  read(12,*)xilock1,xilock2
  read(12,*)sigmadiff,sefffix,Lffix
  read(12,*)tmax
  read(12,*)tslip_ave,tslip_aveint
  read(12,*)tint_out,tmin_out,tint_cos
  read(12,*)vcos,vsse1,vsse2
  read(12,*) tssestart,tsseend
  read(12,*) nmv,nas,ncos,nsse,nslip,nssere
  read(12,*) s1,s2,s3,s4,s5,s6,s7,s8

110 format(A)
  close(12)

   !ALu!ED2
   slipout=0
   sseout=0
   cosout=0
   asout=0
   reout=0

  if(mod(Nt_all,nprocs)/=0)then
     write(*,*)'Nd_all must be integer*nprocs. Change nprocs!'
     STOP
  else
     write(*,*)'Each cpu calculates',Nt_all/nprocs,'cells'
  end if

  ! hnucl = hstarfactor*hd            !h* (how about along strike?)

  if(myid == master)then
     ALLOCATE(x_all(Nt_all),xi_all(Nt_all), &
          cca_all(Nt_all),ccab_all(Nt_all),ccb_all(Nt_all),seff_all(Nt_all),xLf_all(Nt_all), &
          tau_all(Nt_all),slip_all(Nt_all),slipinc_all(Nt_all),yt_all(3*Nt_all), &
          dydt_all(3*Nt_all),yt_scale_all(3*Nt_all))

     ALLOCATE(maxv_s1(nmv),maxv_s2(nmv),maxv_s3(nmv),maxv_s4(nmv),maxv_s5(nmv), &
          maxv_s6(nmv),maxv_s7(nmv),maxv_s8(nmv), &
          maxv(nmv),maxnum(nmv), &
          tmv(nmv),tas(nas),tcos(ncos),tsse(nsse))

!ED2
     ALLOCATE(ttau(nmv),mtau(nmv),mtau_s1(nmv),mtau_s2(nmv),mtau_s3(nmv),mtau_s4(nmv), &
     mtau_s5(nmv),mtau_s6(nmv),mtau_s7(nmv),mtau_s8(nmv),mtaunum(nmv))

!ED5
     ALLOCATE(tseff(nmv),mseff(nmv),mseff_s1(nmv),mseff_s2(nmv),mseff_s3(nmv),mseff_s4(nmv), &
     mseff_s5(nmv),mseff_s6(nmv),mseff_s7(nmv),mseff_s8(nmv),mseffnum(nmv))

!ED4
      ALLOCATE(avev_s(nmv),avev_n(nmv),maxv_s(nmv),maxv_n(nmv),maxnum_s(nmv),maxnum_n(nmv),&
         avetau_s(nmv),avetau_n(nmv),mtau_s(nmv),mtau_n(nmv),mtaunum_s(nmv),mtaunum_n(nmv),&
         vtemp_s(nsouth),vtemp_n(nnorth),tautemp_s(nsouth),tautemp_n(nnorth),&
         sliptemp_s(nsouth),sliptemp_n(nnorth),aveslip_s(nmv),aveslip_n(nmv))
      ALLOCATE(avev_eq(nmv),avev_sse(nmv),maxv_eq(nmv),maxv_sse(nmv),maxnum_eq(nmv),maxnum_sse(nmv),&
         avetau_eq(nmv),avetau_sse(nmv),mtau_eq(nmv),mtau_sse(nmv),mtaunum_eq(nmv),mtaunum_sse(nmv),&
         vtemp_eq(neqzone),vtemp_sse(nssezone),tautemp_eq(neqzone),tautemp_sse(nssezone),&
         sliptemp_eq(neqzone),sliptemp_sse(nssezone),aveslip_eq(nmv),aveslip_sse(nmv),&
         avev_ns(nmv),aveslip_ns(nmv),avetau_ns(nmv),maxv_ns(nmv),mtau_ns(nmv),maxnum_ns(nmv),mtaunum_ns(nmv),&
         vtemp_ns(nnshallow),tautemp_ns(nnshallow),sliptemp_ns(nnshallow))
!ED5
      ALLOCATE(aveseff_s(nmv),aveseff_n(nmv),aveseff_ns(nmv),aveseff_sse(nmv),aveseff_eq(nmv),&
         mseff_s(nmv),mseff_n(nmv),mseff_ns(nmv),mseff_sse(nmv),mseff_eq(nmv), &
         mseffnum_s(nmv),mseffnum_n(nmv),mseffnum_ns(nmv),mseffnum_sse(nmv),mseffnum_eq(nmv), &
         sefftemp_s(nsouth),sefftemp_n(nnorth),sefftemp_ns(nnshallow),sefftemp_sse(nssezone),&
         sefftemp_eq(neqzone))
!ED4
     ALLOCATE(indsouth(nsouth),indnorth(nnorth), indeq(neqzone), indsse(nssezone), indnshallow(nnshallow))
!ALu
     ALLOCATE(tslip(nslip))
!ED3
     ALLOCATE(recenind(nrecenter))
     ALLOCATE(tssere(nssere))
!ED5
     ALLOCATE(re_slip(nrecenter,nssere),re_tau(nrecenter,nssere),re_v(nrecenter,nssere), &
      re_seff(nrecenter, nssere))

!!! modify output number
     ALLOCATE(slipz1_inter(Nt_all,nas),slipz1_cos(Nt_all,ncos), &
          slipz1_tau(Nt_all,nsse),slipz1_sse(Nt_all,nsse))
     ALLOCATE(slipz1_v(Nt_all,ncos))

!ALu !ED2
     ALLOCATE(slip_tau(Nt_all,nslip),slip_v(Nt_all,nslip),slip_slip(Nt_all,nslip))  

  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ALLOCATE (x(Nt),z(Nt),z_all(Nt_all),xi(Nt),cca(Nt),ccb(Nt),seff(Nt),&
       xLf(Nt),tau(Nt),tau0(Nt),slip(Nt),slipinc(Nt), &
       yt(3*Nt),dydt(3*Nt),yt_scale(3*Nt))
  ALLOCATE (stiff(Nt,Nt_all))   !!! stiffness of Stuart green calculation
  ALLOCATE (nstiff(Nt,Nt_all))

  !Read in stiffness matrix, in nprocs segments

  write(cTemp,*) myid
  write(*,*) cTemp

  open(5, file=trim(sstiffname)//'trigreen_'//trim(adjustl(cTemp))//'.bin',form='unformatted',access='stream')
  open(55,file=trim(sstiffname)//'position.bin',form='unformatted',access='stream')
  open(70, file=trim(nstiffname)//'trigreen_'//trim(adjustl(cTemp))//'.bin',form='unformatted',access='stream')
  open(444,file=trim(fricname)//'sigma.txt') !folder
  open(445,file=trim(fricname)//'xLf.txt')
  open(446,file=trim(fricname)//'a.txt')
  open(447,file=trim(fricname)//'ab.txt')
  open(448,file=trim(fricname)//'b.txt')
  record_cor=Nt*(myid)

  !-------------------------------------------------------------------------------------------
  ! read stiffness from Stuart green calculation.
  !-----------------------------------------------------------------------------------------
  if(myid==master)then
     do k=1,Nt_all
        read(55) xi_all(k),x_all(k),z_all(k) !xi is along the dip while x is along the strike
        xi_all(k)=xi_all(k)/1000
        x_all(k)=x_all(k)/1000
        z_all(k)=z_all(k)/1000

        read(444,*) seff_all(k)
        read(445,*) xLf_all(k)
        read(446,*) cca_all(k)
        read(447,*) ccab_all(k)
        read(448,*) ccb_all(k)
     end do
      
      open(2,file=trim(foldername)//'vardep'//jobname,status='unknown')
      !	write(2,300)'z','seff','Lf','ccab','cca'
      do i=1,Nt_all
         write(2,'(5(1x,e20.13))')z_all(i),seff_all(i),xLf_all(i), &
              ccab_all(i),cca_all(i)
      end do
      close(2)

      !ED3
      open(13,file=trim(patchname)//'recenter.txt',form='formatted',status='old')
      do i=1,nrecenter
         read(13,*) recenind(i)
      end do
      close(13)

      !ED4
      open(14,file=trim(patchname)//'indsd.txt',form='formatted',status='old')
      do i=1,nsouth
         read(14,*) indsouth(i)
      end do
      close(14)     

      open(15,file=trim(patchname)//'indnd.txt',form='formatted',status='old')
      do i=1,nnorth
         read(15,*) indnorth(i)
      end do
      close(15)

      open(16,file=trim(patchname)//'indeq.txt',form='formatted',status='old')
      do i=1,neqzone
         read(16,*) indeq(i)
      end do        
      close(16)   

      open(17,file=trim(patchname)//'indsse.txt',form='formatted',status='old')
      do i=1,nssezone
         read(17,*) indsse(i)
      end do        
      close(17)  

      open(18,file=trim(patchname)//'indns.txt',form='formatted',status='old')
      do i=1,nnshallow
         read(18,*) indnshallow(i)
      end do        
      close(18)



  end if

  do i=1,Nt !! observe
     do j=1,Nt_all !! source
        read(5) stiff(i,j)

        if(stiff(i,j).lt.-0.5d0.or.stiff(i,j).gt.0.99d0)then
           stiff(i,j) = 0.d0
           write(*,*) myid,i,'shear stiffness extreme'
	end if
        if(isNaN(stiff(i,j)))then
           stiff(i,j)=0.d0
           write(*,*) myid,i,'shear stiffness extreme'
        end if

     end do
  end do
  close(5)

  do i=1,Nt !! observe
     do j=1,Nt_all !! source
        read(70) nstiff(i,j)
        if(nstiff(i,j).lt.-0.5d0.or.nstiff(i,j).gt.0.99d0)then
           nstiff(i,j) = 0.d0
           write(*,*) myid,i,j,'normal stiffness extreme'
	     end if
        if(isNaN(nstiff(i,j)))then
           nstiff(i,j)=0.d0
           write(*,*) myid,i,'normal stiffness NaN'
        end if
     end do
  end do


  !ALu 
  close(70)

  close(55) 
  close(444)
  close(445)
  close(446)
  close(447)
  close(448)
  !!-----------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------
  
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  call MPI_Scatter(cca_all,Nt,MPI_Real8,cca,Nt,MPI_Real8,master,MPI_COMM_WORLD,ierr)
  call MPI_Scatter(ccb_all,Nt,MPI_Real8,ccb,Nt,MPI_Real8,master,MPI_COMM_WORLD,ierr)
  call MPI_Scatter(xLf_all,Nt,MPI_Real8,xLf,Nt,MPI_Real8,master,MPI_COMM_WORLD,ierr)
  call MPI_Scatter(seff_all,Nt,MPI_Real8,seff,Nt,MPI_Real8,master,MPI_COMM_WORLD,ierr)
  call MPI_Scatter(x_all,Nt,MPI_Real8,x,Nt,MPI_Real8,master,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(z_all,Nt_all,MPI_Real8,master,MPI_COMM_WORLD,ierr)

  call CPU_TIME(tmbegin)

  tm1=tmbegin
  tmday=86400.0
  tmelse=0.0
  tmmult=0.0
  tmmidn=0.0

  t=0.d0
  tprint_inter = 0.d0
  tslipcos = 0.d0
  tslipsse = 0.d0

  !ALu 
  tout =0.d0

  ndtnext=0
  ndt_v=0
  ndt_inter=0
  ndt =0
  nrec=0
  i_nul = 0
  inul = 0
  isse = 0

  Ioutput = 0    !initial always 0 (output)

  imv=0   !counter for maxv, sliptot, slipthresh1, slipthresh2,slipthresh3 output
  ias=0 !counter for slip at iz3 and s.z. average slip output
  icos = 0
  
  !ALu !ED3
  islip=0
  issere=0
  itau=0 !ED5S
  iseff=0

  accuracy = 1.d-3
  epsv = accuracy
  dtmin = 1.d-6
  dt_try=1.d-6


  if(myid==master)then

     open(1,file=trim(foldername)//'phypara'//jobname, status='unknown')
     write(1,*)'procs num = ', nprocs

  end if

  Ifileout = 60   !file index, after 47

  !----Initial values of velocity, state variable, shear stress and slip--
  !--SET INITIAL VPL FOR THE LOCKED PART TO BE 0
  ! ! set plate convergence
  if(IDin.eq.0)then
     do j=1,Nt
!ALu      
        yt(3*j-2)=seff(j) !inital normal stress
        yt(3*j-1)=Vpl !try different !initial velocity
        yt(3*j)=xLf(j)/(1.1*yt(3*j-1)) !initial state variable
        slip(j)=0.d0 
     end do
  end if

  !------------------------------------------------------------------
  if(IDin.eq.1) then               !if this is a restart job
     if(myid==master)then
        call restart(0,'out',4,Nt_all,t,dt,dt_try,ndt,nrec,yt_all,slip_all)
        write(1,*)'This is a restart job. Start time ',t,' yr'
     end if
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     call MPI_Bcast(t,1,MPI_Real8,master,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(dt,1,MPI_Real8,master,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(dt_try,1,MPI_Real8,master,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(ndt,1,MPI_integer,master,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(nrec,1,MPI_integer,master,MPI_COMM_WORLD,ierr)
     call MPI_Scatter(yt_all,3*Nt,MPI_Real8,yt,3*Nt,MPI_Real8,master,MPI_COMM_WORLD,ierr)
     call MPI_Scatter(slip_all,Nt,MPI_Real8,slip,Nt,MPI_Real8,master,MPI_COMM_WORLD,ierr)

     ndtnext = ndt
     tprint_inter = t
     tslip_ave=t
     tout = t

  else
     if(myid==master)then
        write(1,*)'Start time ',t,' yr'
     end if
  end if
  if(myid==master)then
     close(1)
  end if

  !----------------------------------------------
  !      Start of Basic Cycle:
  !----------------------------------------------
  cyclecont=.true.

  do while(cyclecont)

!ALu
     call derivs(myid,dydt,3*Nt,Nt_all,Nt,t,yt,z_all,x)
   
     do j=1,3*Nt
        yt_scale(j)=dabs(yt(j))+dabs(dt_try*dydt(j))
     end do
!ALu 
     CALL rkqs(myid,yt,dydt,3*Nt,Nt_all,Nt,t,dt_try,accuracy,yt_scale, &
          dt_did,dt_next,z_all,x) 
     dt = dt_did
     dt_try = dt_next

!ALu 
     do i=1,Nt
        if(yt(3*i-1).lt.epsv)then
           help=(yt(3*i-1)/(2*V0))*dexp((f0+ccb(i)*dlog(V0*yt(3*i)))/cca(i))
           tau(i)=yt(3*i-2)*cca(i)*dlog(help+dsqrt(1+help**2))
        else
           tau(i)=yt(3*i-2)*(f0+cca(i)*log(yt(3*i-1)/v0)+ccb(i)*log(V0*yt(3*i)/xLf(i)))
        end if
        slipinc(i) = yt(3*i-1)*dt
        slip(i) = slip(i) + slipinc(i)
     end do

     ndt = ndt + 1

     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     call MPI_Gather(yt,3*Nt,MPI_Real8,yt_all,3*Nt,MPI_Real8,master,MPI_COMM_WORLD,ierr)
     call MPI_Gather(slipinc,Nt,MPI_Real8,slipinc_all,Nt,MPI_Real8,master,MPI_COMM_WORLD,ierr)
     call MPI_Gather(slip,Nt,MPI_Real8,slip_all,Nt,MPI_Real8,master,MPI_COMM_WORLD,ierr)
     call MPI_Gather(tau,Nt,MPI_Real8,tau_all,Nt,MPI_Real8,master,MPI_COMM_WORLD,ierr)

     !-------------------
     !      Output:     (a single thread will do the writing while others
     !      proceed
     !-------------------

     if(myid==master)then
        imv=imv+1
        tmv(imv)=t

!!! Initialize max V !ED5
        maxv(imv) = 0.d0
        maxv_s1(imv)=0.d0
        maxv_s2(imv)=0.d0
        maxv_s3(imv)=0.d0
        maxv_s4(imv)=0.d0
        maxv_s5(imv)=0.d0
        maxv_s6(imv)=0.d0
        maxv_s7(imv)=0.d0
        maxv_s8(imv)=0.d0

        maxv_s(imv)=0.d0 !ED4 max V for different patches
        maxv_n(imv)=0.d0
        maxv_ns(imv)=0.d0
        maxv_eq(imv)=0.d0
        maxv_sse(imv)=0.d0

!ED2
!!! Initialize max shear stress !ED5
        itau=itau+1
        ttau(itau)=t
        mtau(itau)=0.d0
        mtau_s1(itau)=0.d0
        mtau_s2(itau)=0.d0
        mtau_s3(itau)=0.d0
        mtau_s4(itau)=0.d0
        mtau_s5(itau)=0.d0
        mtau_s6(itau)=0.d0
        mtau_s7(itau)=0.d0
        mtau_s8(itau)=0.d0

        mtau_s(itau)=0.d0 !ED4 
        mtau_n(itau)=0.d0
        mtau_ns(itau)=0.d0
        mtau_eq(itau)=0.d0
        mtau_sse(itau)=0.d0

!ED5
!!! Initialize max normal stress

         iseff=iseff+1
         tseff(iseff)=t
         mseff(iseff)=0.d0
         mseff_s1(iseff)=0.d0
         mseff_s2(iseff)=0.d0
         mseff_s3(iseff)=0.d0
         mseff_s4(iseff)=0.d0
         mseff_s5(iseff)=0.d0
         mseff_s6(iseff)=0.d0
         mseff_s7(iseff)=0.d0
         mseff_s8(iseff)=0.d0

         mseff_s(iseff)=0.d0 !ED4 
         mseff_n(iseff)=0.d0
         mseff_ns(iseff)=0.d0
         mseff_eq(iseff)=0.d0
         mseff_sse(iseff)=0.d0

! Find Max values    
!ALu
        do i = 1,Nt_all      !!!!!!!!!!!!!!!!! find max velocity
           if(yt_all(3*i-1).ge.maxv(imv))then
              maxv(imv) = yt_all(3*i-1)
              maxnum(imv)=i
           end if
        end do
        write(*,*) ndt,t,maxv(imv)

        !! write fault observations
              maxv_s1(imv) = yt_all(3*s1-1)
              maxv_s2(imv) = yt_all(3*s2-1)
              maxv_s3(imv) = yt_all(3*s3-1)
              maxv_s4(imv) = yt_all(3*s4-1)
              maxv_s5(imv) = yt_all(3*s5-1)
              maxv_s6(imv) = yt_all(3*s6-1)
              maxv_s7(imv) = yt_all(3*s7-1)
              maxv_s8(imv) = yt_all(3*s8-1)         

        do i = 1,Nt_all      !!!!!!!!!!!!!!!!! find max shear stress
           if(tau_all(i).ge.mtau(itau))then
              mtau(itau) = tau_all(i)
              mtaunum(itau)=i
           end if
        end do
        write(*,*) ndt,t,mtau(itau)

        mtau_s1(itau) = tau_all(s1)
        mtau_s2(itau) = tau_all(s2)
        mtau_s3(itau) = tau_all(s3)
        mtau_s4(itau) = tau_all(s4)
        mtau_s5(itau) = tau_all(s5)
        mtau_s6(itau) = tau_all(s6)
        mtau_s7(itau) = tau_all(s7)
        mtau_s8(itau) = tau_all(s8)

         do i = 1,Nt_all      !!!!!!!!!!!!!!!!! find max effective normal stress
            if(yt_all(3*i-2).ge.mseff(iseff))then
               mseff(iseff) = yt_all(3*i-2)
               mseffnum(iseff)=i
            end if
         end do
         write(*,*) ndt,t,mseff(iseff)

        !! write fault observations
         mseff_s1(iseff) = yt_all(3*s1-2)
         mseff_s2(iseff) = yt_all(3*s2-2)
         mseff_s3(iseff) = yt_all(3*s3-2)
         mseff_s4(iseff) = yt_all(3*s4-2)
         mseff_s5(iseff) = yt_all(3*s5-2)
         mseff_s6(iseff) = yt_all(3*s6-2)
         mseff_s7(iseff) = yt_all(3*s7-2)
         mseff_s8(iseff) = yt_all(3*s8-2)

!!! Find max Velocity for different patches
        do i=1,nsouth !South
            if(yt_all(3*indsouth(i)-1).ge.maxv_s(imv))then
               maxv_s(imv)=yt_all(3*indsouth(i)-1)
               maxnum_s(imv)=indsouth(i)
            end if
        end do  

        do i=1,nnorth !North
            if(yt_all(3*indnorth(i)-1).ge.maxv_n(imv))then
               maxv_n(imv)=yt_all(3*indnorth(i)-1)
               maxnum_n(imv)=indnorth(i)
            end if
        end do      

        do i=1,neqzone !Earthquake zone
            if(yt_all(3*indeq(i)-1).ge.maxv_eq(imv))then
               maxv_eq(imv)=yt_all(3*indeq(i)-1)
               maxnum_eq(imv)=indeq(i)
            end if
        end do             

        do i=1,nssezone !SSE zone
            if(yt_all(3*indsse(i)-1).ge.maxv_sse(imv))then
               maxv_sse(imv)=yt_all(3*indsse(i)-1)
               maxnum_sse(imv)=indsse(i)
            end if
        end do 

        do i=1,nnshallow !North Shallow
            if(yt_all(3*indnshallow(i)-1).ge.maxv_ns(imv))then
               maxv_ns(imv)=yt_all(3*indnshallow(i)-1)
               maxnum_ns(imv)=indnshallow(i)
            end if
        end do 

!!! Find max normal stress for different patches
        do i=1,nsouth !South
            if(yt_all(3*indsouth(i)-2).ge.mseff_s(imv))then
               mseff_s(imv)=yt_all(3*indsouth(i)-2)
               mseffnum_s(imv)=indsouth(i)
            end if
        end do  

        do i=1,nnorth !North
            if(yt_all(3*indnorth(i)-2).ge.mseff_n(imv))then
               mseff_n(imv)=yt_all(3*indnorth(i)-2)
               mseffnum_n(imv)=indnorth(i)
            end if
        end do      

        do i=1,neqzone !Earthquake zone
            if(yt_all(3*indeq(i)-2).ge.mseff_eq(imv))then
               mseff_eq(imv)=yt_all(3*indeq(i)-2)
               mseffnum_eq(imv)=indeq(i)
            end if
        end do             

        do i=1,nssezone !SSE zone
            if(yt_all(3*indsse(i)-2).ge.mseff_sse(imv))then
               mseff_sse(imv)=yt_all(3*indsse(i)-2)
               mseffnum_sse(imv)=indsse(i)
            end if
        end do 

        do i=1,nnshallow !North Shallow
            if(yt_all(3*indnshallow(i)-2).ge.mseff_ns(imv))then
               mseff_ns(imv)=yt_all(3*indnshallow(i)-2)
               mseffnum_ns(imv)=indnshallow(i)
            end if
        end do 

!!! Find max shear stress for different patches
        do i=1,nsouth
           if(tau_all(indsouth(i)).ge.mtau_s(itau))then
               mtau_s(itau)=tau_all(indsouth(i))
               mtaunum_s(itau)=indsouth(i)
           end if
        end do 
        do i=1,nnorth
           if(tau_all(indnorth(i)).ge.mtau_n(itau))then
               mtau_n(itau)=tau_all(indnorth(i))
               mtaunum_n(itau)=indnorth(i)
           end if
        end do 

        do i=1,neqzone
           if(tau_all(indeq(i)).ge.mtau_eq(itau))then
               mtau_eq(itau)=tau_all(indeq(i))
               mtaunum_eq(itau)=indeq(i)
           end if
        end do 

        do i=1,nssezone
           if(tau_all(indsse(i)).ge.mtau_sse(itau))then
               mtau_sse(itau)=tau_all(indsse(i))
               mtaunum_sse(itau)=indsse(i)
           end if
        end do 

        do i=1,nnshallow
           if(tau_all(indnshallow(i)).ge.mtau_ns(itau))then
               mtau_ns(itau)=tau_all(indnshallow(i))
               mtaunum_ns(itau)=indnshallow(i)
           end if
        end do 

!!! Find the Average Velocity for different patches
         vsum_s=0.d0 !South
         do i=1,nsouth
            vtemp_s(i)=yt_all(3*indsouth(i)-1)
            vsum_s=vsum_s+vtemp_s(i)
         end do 
         avev_s(imv)=vsum_s/nsouth

!         write(*,*) 'check ave'
!         write(*,*) sum(vtemp_s), nsouth
!         write(*,*) 'ave',avev_s(imv)

         vsum_n=0.d0 !North
         do i=1,nnorth
            vtemp_n(i)=yt_all(3*indnorth(i)-1)
            vsum_n=vsum_n+vtemp_n(i)
         end do 
         avev_n(imv)=vsum_n/nnorth

         vsum_eq=0.d0 !Earthquake zone
         do i=1,neqzone
            vtemp_eq(i)=yt_all(3*indeq(i)-1)
            vsum_eq=vsum_eq+vtemp_eq(i)
         end do 
         avev_eq(imv)=vsum_eq/neqzone         

         vsum_sse=0.d0 !SSE zone
         do i=1,nssezone
            vtemp_sse(i)=yt_all(3*indsse(i)-1)
            vsum_sse=vsum_sse+vtemp_sse(i)
         end do 
         avev_sse(imv)=vsum_sse/nssezone  

         vsum_ns=0.d0 !North Shallow
         do i=1,nnshallow
            vtemp_ns(i)=yt_all(3*indnshallow(i)-1)
            vsum_ns=vsum_ns+vtemp_ns(i)
         end do 
         avev_ns(imv)=vsum_ns/nnshallow  

!!! Find the average normal stress for different patches
         seffsum_s=0.d0 !South
         do i=1,nsouth
            sefftemp_s(i)=yt_all(3*indsouth(i)-2)
            seffsum_s=seffsum_s+sefftemp_s(i)
         end do 
         aveseff_s(imv)=seffsum_s/nsouth

         write(*,*) 'check normal ave'
         write(*,*) sum(sefftemp_s), nsouth
         write(*,*) 'ave',aveseff_s(imv)

         seffsum_n=0.d0 !North
         do i=1,nnorth
            sefftemp_n(i)=yt_all(3*indnorth(i)-2)
            seffsum_n=seffsum_n+sefftemp_n(i)
         end do 
         aveseff_n(imv)=seffsum_n/nnorth

         seffsum_eq=0.d0 !Earthquake zone
         do i=1,neqzone
            sefftemp_eq(i)=yt_all(3*indeq(i)-2)
            seffsum_eq=seffsum_eq+sefftemp_eq(i)
         end do 
         aveseff_eq(imv)=seffsum_eq/neqzone         

         seffsum_sse=0.d0 !SSE zone
         do i=1,nssezone
            sefftemp_sse(i)=yt_all(3*indsse(i)-2)
            seffsum_sse=seffsum_sse+sefftemp_sse(i)
         end do 
         aveseff_sse(imv)=seffsum_sse/nssezone  

         seffsum_ns=0.d0 !North Shallow
         do i=1,nnshallow
            sefftemp_ns(i)=yt_all(3*indnshallow(i)-2)
            seffsum_ns=seffsum_ns+sefftemp_ns(i)
         end do 
         aveseff_ns(imv)=seffsum_ns/nnshallow  


!!! Find the Average Shear Stress for North and South Patches
         tausum_s=0.d0
         do i=1,nsouth
            tautemp_s(i)=tau_all(indsouth(i))
            tausum_s=tausum_s+tautemp_s(i)
         end do 
         avetau_s(imv)=tausum_s/nsouth

         tausum_n=0.d0
         do i=1,nnorth
            tautemp_n(i)=tau_all(indnorth(i))
            tausum_n=tausum_n+tautemp_n(i)
         end do  
         avetau_n(imv)=tausum_n/nnorth

         tausum_eq=0.d0
         do i=1,neqzone
            tautemp_eq(i)=tau_all(indeq(i))
            tausum_eq=tausum_eq+tautemp_eq(i)
         end do  
         avetau_eq(imv)=tausum_eq/neqzone

         tausum_sse=0.d0
         do i=1,nssezone
            tautemp_sse(i)=tau_all(indsse(i))
            tausum_sse=tausum_sse+tautemp_sse(i)
         end do  
         avetau_sse(imv)=tausum_sse/nssezone     

         tausum_ns=0.d0
         do i=1,nnshallow
            tautemp_ns(i)=tau_all(indnshallow(i))
            tausum_ns=tausum_ns+tautemp_ns(i)
         end do  
         avetau_ns(imv)=tausum_ns/nnshallow      

!!! Find the Average Slip for North and South Patches
         slipsum_s=0.d0
         do i=1,nsouth
            sliptemp_s(i)=slip_all(indsouth(i))
            slipsum_s=slipsum_s+sliptemp_s(i)
         end do
         aveslip_s(imv)=slipsum_s/nsouth

         slipsum_n=0.d0
         do i=1,nnorth
            sliptemp_n(i)=slip_all(indnorth(i))
            slipsum_n=slipsum_n+sliptemp_n(i)
         end do
         aveslip_n(imv)=slipsum_n/nnorth         

         slipsum_eq=0.d0
         do i=1,neqzone
            sliptemp_eq(i)=slip_all(indeq(i))
            slipsum_eq=slipsum_eq+sliptemp_eq(i)
         end do
         aveslip_eq(imv)=slipsum_eq/neqzone

         slipsum_sse=0.d0
         do i=1,nssezone
            sliptemp_sse(i)=slip_all(indsse(i))
            slipsum_sse=slipsum_sse+sliptemp_sse(i)
         end do
         aveslip_sse(imv)=slipsum_sse/nssezone   

         slipsum_ns=0.d0
         do i=1,nnshallow
            sliptemp_ns(i)=slip_all(indnshallow(i))
            slipsum_ns=slipsum_ns+sliptemp_ns(i)
         end do
         aveslip_ns(imv)=slipsum_ns/nnshallow  

!!!!! checkout point  - not necessary now !!!!!!!!!!!!!!!

        !-----Interseismic slip every ? years----

        if (t.ge.tslip_ave)then
           ias = ias + 1
           tas(ias)=t

           do i=1,Nt_all
              slipz1_inter(i,ias) = slip_all(i)*(1.d-3)
           end do

           tslip_ave = tslip_ave + tslip_aveint
        end if


!ALu
       !!!!!!!!!!!!!!!!---Find SSE and Coseismic Slip------------
         if(t.ge.tssestart.and.t.le.tsseend)then
            if((maxv(imv)/yrs).le.vcos)then
               tslipsse = tslipsse + dt
               tcsse=tcsse+dt
               if(tslipsse.ge.0.002739726)then   !! output step = 1 day
 		            write(*,*)'SSE',t,maxv(imv)
                  isse = isse + 1
                  islip=islip+1
                  tslip(islip)=t
                  tsse(isse) = t
                  
                  do i=1,Nt_all
                     slipz1_tau(i,isse)= tau_all(i)
                     slipz1_sse(i,isse)= dlog10(yt_all(3*i-1)*(1.d-3)/yrs)

                     slip_tau(i,islip)=tau_all(i)
                     slip_v(i,islip)= dlog10(yt_all(3*i-1)*(1.d-3)/yrs)
                     slip_slip(i,islip)=slip_all(i)*(1.d-3)
                  end do
                  tslipsse = 0.d0

               end if
               
!ED3

               if(tcsse.ge.0.0137)then !! output step = 5 days
                  write(*,*)'reduced SSE output',t,maxv(imv)
                  issere=issere+1
                  tssere(issere)=t
                  do i=1,nrecenter
                     indtemp=recenind(i)
                     re_slip(i,issere)=slip_all(indtemp)*(1.d-3)
                     re_v(i,issere)=dlog10(yt_all(3*indtemp-1)*(1.d-3)/yrs)
                     re_tau(i,issere)=tau_all(indtemp)
                     re_seff(i,issere)=yt_all(3*indtemp-2)
                  end do
                  tcsse=0.d0
               end if
            end if 
         end if
              
            if((maxv(imv)/yrs).ge.vcos)then
               tslipcos = tslipcos+dt
               if(tslipcos.ge.tint_cos)then    !!!!tint_cos: every 5s output
                  write(*,*) 'coseis',t,maxv(imv)
                  icos = icos +1
                  islip=islip+1
                  tslip(islip)=t
                  tcos(icos) = t

                  do i=1,Nt_all
                     slipz1_cos(i,icos) = slip_all(i)*(1.d-3)
                     slipz1_v(i,icos) = dlog10(yt_all(3*i-1)*(1.d-3)/yrs)

                     slip_tau(i,islip)=tau_all(i)
                     slip_v(i,islip)=dlog10(yt_all(3*i-1)*(1.d-3)/yrs)
                     slip_slip(i,islip)=slip_all(i)*(1.d-3)
                  end do 
!ED3
                  write(*,*)'reduced SSE output',t,maxv(imv)
                  issere=issere+1
                  tssere(issere)=t
                  do i=1,nrecenter
                     indtemp=recenind(i)
                     re_slip(i,issere)=slip_all(indtemp)*(1.d-3)
                     re_v(i,issere)=dlog10(yt_all(3*indtemp-1)*(1.d-3)/yrs)
                     re_tau(i,issere)=tau_all(indtemp)
                     re_seff(i,issere)=yt_all(3*indtemp-2)
                  end do


                  tslipcos = 0.d0
               end if
            end if 
     end if


!!! find the total number of output elements and cell number of specific depth.

!              do i=1,Nt_all
!                 slipz1_cos(i,icos) = slip_all(i)*(1.d-3)
!                 slipz1_v(i,icos) = dlog10(yt_all(3*i-1)*(1.d-3)/yrs)
!              end do

!              tslipcos = 0.d0
!           end if
!        end if
!     end if

     !----Output restart files -------------
     if(myid==master)then
        if(mod(ndt,1000).eq.0)ihrestart=1
        if(IDout.eq.1.and.ihrestart.eq.1)then
	   filename='out0'
           call restart(1,filename,4,Nt_all,t,dt,dt_try,ndt,nrec,yt_all,slip_all)
           ihrestart=0
        end if
        if(abs(t-tout).le.tmin_out)then
           ihrestart = 1
           Itout=int(tout)
           write(ct,*)Itout
           ct=adjustl(ct)
           filename='out'//trim(ct)
           call restart(1,filename,Ifileout,Nt_all,t,dt,dt_try,ndt,nrec,yt_all,slip_all)
           ihrestart = 0
           tout = tout+tint_out
        end if
     end if

     !----Output velocity and slip records ---
     !----velocity in mm/yr ; slip in meter, moment in 10^{14} Nm --
     if(myid==master)then
        Ioutput = 0


        !call output(Ioutput,Isnapshot,Nt_all,Nt,inul,imv,ias,icos,isse,x,&
        !     tmv,tas,tcos,tsse,maxv,maxv_s1,maxv_s2,maxv_s3,maxv_s4,maxv_s5,maxv_s6,&
        !     maxv_s7,maxv_s8,maxnum, &
        !     slipz1_inter,slipz1_tau,slipz1_sse,slipz1_cos, xi_all,x_all,slipz1_v)

        !ALu
        !call output(Ioutput,Isnapshot,Nt_all,Nt,inul,imv,ias,icos,isse,islip,x,&
        !     tmv,tas,tcos,tsse,tslip,maxv,maxv_s1,maxv_s2,maxv_s3,maxv_s4,maxv_s5,maxv_s6,&
        !     maxv_s7,maxv_s8,maxnum, &
        !     slipz1_inter,slipz1_tau,slipz1_sse,slipz1_cos,&
        !     slip_tau,slip_v,slip_slip,xi_all,x_all,slipz1_v,&
        !     slipout, sseout, cosout, asout)
        !ED3 !ED4
        call output(Ioutput,Isnapshot,Nt_all,Nt,inul,imv,ias,icos,isse,islip,x,&
             tmv,tas,tcos,tsse,tslip,maxv,maxv_s1,maxv_s2,maxv_s3,maxv_s4,maxv_s5,maxv_s6,&
             maxv_s7,maxv_s8,maxnum, &
             slipz1_inter,slipz1_tau,slipz1_sse,slipz1_cos,&
             slip_tau,slip_v,slip_slip,xi_all,x_all,slipz1_v,&
             slipout, sseout, cosout, asout, &
             itau, ttau,mtau, mtau_s1,mtau_s2,mtau_s3,mtau_s4,mtau_s5,mtau_s6,mtau_s7,mtau_s8,mtaunum,&
             issere,tssere,re_slip,re_tau,re_v,nrecenter,reout, re_seff, & 
             maxnum_s, maxnum_n, maxv_s, maxv_n, avev_s, avev_n, aveslip_s, aveslip_n,&
             mtaunum_s, mtaunum_n, mtau_s, mtau_n, avetau_s, avetau_n,&
             maxnum_eq, maxnum_sse, maxv_eq, maxv_sse, avev_eq, avev_sse, aveslip_eq, aveslip_sse,&
             mtaunum_eq, mtaunum_sse, mtau_eq, mtau_sse, avetau_eq, avetau_sse,&
             maxnum_ns, maxv_ns, avev_ns, aveslip_ns, mtaunum_ns, mtau_ns, avetau_ns, &
             iseff, tseff, mseff, mseff_s1, mseff_s2, mseff_s3, mseff_s4, mseff_s5, mseff_s6,&
             mseff_s7, mseff_s8, mseffnum, &
             mseff_s, mseff_n, mseff_ns, mseff_sse, mseff_eq, &
             mseffnum_s, mseffnum_n, mseffnum_ns, mseffnum_sse, mseffnum_eq, &
             aveseff_s, aveseff_n, aveseff_ns, aveseff_sse, aveseff_eq )


     end if

     !      end of output
     !----------------------------------------------------
     !   End of output conmmands.
     !----------------------------------------------------

     !----------------------------------------------------
     !   Go to next cycle if t < tmax or ndt > ndtmax
     !---------------------------------------------------

     if (t>tmax)cyclecont = .false.

  end do

  !--- Final output -------

  if(myid==master)then
     filename='outlast'
     call restart(1,filename,Ifileout,Nt_all,t,dt,dt_try,ndt,nrec,yt_all,slip_all)
     Ioutput = 1
!     call output(Ioutput,Isnapshot,Nt_all,Nt,inul,imv,ias,icos,isse,x,&
!          tmv,tas,tcos,tsse,maxv,maxv_s1,maxv_s2,maxv_s3,maxv_s4,maxv_s5,maxv_s6,maxv_s7,maxv_s8, &
!          maxnum, &
!          slipz1_inter,slipz1_tau,slipz1_sse, &
!          slipz1_cos,&
!          xi_all,x_all,slipz1_v)

      call output(Ioutput,Isnapshot,Nt_all,Nt,inul,imv,ias,icos,isse,islip,x,&
         tmv,tas,tcos,tsse,tslip,maxv,maxv_s1,maxv_s2,maxv_s3,maxv_s4,maxv_s5,maxv_s6,&
         maxv_s7,maxv_s8,maxnum, &
         slipz1_inter,slipz1_tau,slipz1_sse,slipz1_cos,&
         slip_tau,slip_v,slip_slip,xi_all,x_all,slipz1_v,&
         slipout, sseout, cosout, asout, &
         itau, ttau,mtau, mtau_s1,mtau_s2,mtau_s3,mtau_s4,mtau_s5,mtau_s6,mtau_s7,mtau_s8,mtaunum,&
         issere,tssere,re_slip,re_tau,re_v,nrecenter,reout, re_seff, & 
         maxnum_s, maxnum_n, maxv_s, maxv_n, avev_s, avev_n, aveslip_s, aveslip_n,&
         mtaunum_s, mtaunum_n, mtau_s, mtau_n, avetau_s, avetau_n,&
         maxnum_eq, maxnum_sse, maxv_eq, maxv_sse, avev_eq, avev_sse, aveslip_eq, aveslip_sse,&
         mtaunum_eq, mtaunum_sse, mtau_eq, mtau_sse, avetau_eq, avetau_sse,&
         maxnum_ns, maxv_ns, avev_ns, aveslip_ns, mtaunum_ns, mtau_ns, avetau_ns, &
         iseff, tseff, mseff, mseff_s1, mseff_s2, mseff_s3, mseff_s4, mseff_s5, mseff_s6,&
         mseff_s7, mseff_s8, mseffnum, &
         mseff_s, mseff_n, mseff_ns, mseff_sse, mseff_eq, &
         mseffnum_s, mseffnum_n, mseffnum_ns, mseffnum_sse, mseffnum_eq, &
         aveseff_s, aveseff_n, aveseff_ns, aveseff_sse, aveseff_eq )

!ED2
!     call output(Ioutput,Isnapshot,Nt_all,Nt,inul,imv,ias,icos,isse,islip,x,&
!          tmv,tas,tcos,tsse,tslip,maxv,maxv_s1,maxv_s2,maxv_s3,maxv_s4,maxv_s5,maxv_s6,maxv_s7,maxv_s8, &
!          maxnum, &
!          slipz1_inter,slipz1_tau,slipz1_sse, &
!          slipz1_cos,slip_tau,slip_v,slip_slip, &
!          xi_all,x_all,slipz1_v,&
!         slipout, sseout, cosout, asout, &
!             itau, ttau,mtau, mtau_s1,mtau_s2,mtau_s3,mtau_s4,mtau_s5,mtau_s6,mtau_s7,mtau_s8,mtaunum,&
!             issere,tssere,re_slip,re_tau,re_v,nrecenter,reout)
  end if

  !---End of final output ----

  call CPU_TIME(tmrun)
  tmrun = tmmidn*tmday + tmrun - tmbegin
  tmmult = tmmult/tmrun*100.
  tmelse = tmelse/tmrun*100.
  if(myid==master)then
     open(10,file=trim(foldername)//'summary'//jobname,status='unknown')

     write(10,*)'processor',myid
     write(10,*)'      PARTS OF RUNNING TIME     '
     write(10,'(T9,A,T40,F10.5)')'message passing (percent)', tmmult
     write(10,'(T9,A,T40,F10.5)')'everything else (percent)', tmelse
     write(10,*)
     write(10,'(T9,A,T45,F10.5)')'total running time(min)', tmrun/60.0
     write(10,*)
     write(10,*)'       INFORMATION ABOUT THE END OF THE RUN      '
     write(10,'(T9,A,I7)')'ndtend = ', ndt
     write(10,'(T9,A,D20.13)')'tend = ',t
     write(10,'(T9,A,I7)')'nrec = ',nrec
     write(10,*)
     write(10,*)'Nprocs=', nprocs
     close(10)
  end if


  if(myid==master)then
     DEALLOCATE (x_all,xi_all,yt_all,dydt_all,yt_scale_all,tau_all, &
          slip_all,slipinc_all,cca_all,ccb_all,xLf_all,seff_all, &
          maxnum,maxv,maxv_s1,maxv_s2,maxv_s3,maxv_s4,maxv_s5,maxv_s6,maxv_s7,maxv_s8,&
          tmv,tcos,tas,tsse)
!ED2
     DEALLOCATE (ttau,mtau,mtau_s1,mtau_s2,mtau_s3,mtau_s4,mtau_s5,mtau_s6,mtau_s7,mtau_s8, &
          mtaunum)

!ED5 
      DEALLOCATE (tseff,mseff,mseff_s1,mseff_s2,mseff_s3,mseff_s4,mseff_s5,mseff_s6, &
      mseff_s7,mseff_s8,mseffnum)

!ED4
      DEALLOCATE(indsouth,indnorth,indeq,indsse,indnshallow)
      DEALLOCATE(avev_s,avev_n,maxv_s,maxv_n,maxnum_s,maxnum_n,&
      avetau_s,avetau_n,mtau_s,mtau_n,mtaunum_s,mtaunum_n,aveslip_s,aveslip_n)
      DEALLOCATE(vtemp_s,vtemp_n,tautemp_s,tautemp_n,sliptemp_s,sliptemp_n)
      DEALLOCATE(avev_eq,avev_sse,maxv_eq,maxv_sse,maxnum_eq,maxnum_sse,&
      avetau_eq,avetau_sse,mtau_eq,mtau_sse,mtaunum_eq,mtaunum_sse,aveslip_eq,aveslip_sse)
      DEALLOCATE(vtemp_eq,vtemp_sse,tautemp_eq,tautemp_sse,sliptemp_eq,sliptemp_sse)
      DEALLOCATE(avev_ns,maxv_ns,maxnum_ns,&
      avetau_ns,mtau_ns,mtaunum_ns,aveslip_ns)
      DEALLOCATE(vtemp_ns,tautemp_ns,sliptemp_ns)
   !ED5
      DEALLOCATE(aveseff_s,aveseff_n,aveseff_ns,aveseff_sse,aveseff_eq, &
      mseff_s,mseff_n,mseff_ns,mseff_sse,mseff_eq, &
      mseffnum_s,mseffnum_n,mseffnum_ns,mseffnum_sse,mseffnum_eq, &
      sefftemp_s,sefftemp_n,sefftemp_ns,sefftemp_sse,sefftemp_eq)


     DEALLOCATE (slipz1_inter,slipz1_tau,slipz1_sse, &
          slipz1_cos,slipz1_v)
!     DEALLOCATE (intdepz1,intdepz2,intdepz3,ssetime,slipz1_v)
!ALu
     DEALLOCATE (tslip)
     DEALLOCATE (slip_tau,slip_v,slip_slip)
     !ED3
     DEALLOCATE(recenind)
     DEALLOCATE(tssere)
!ED5
     DEALLOCATE(re_slip,re_tau,re_v,re_seff)

  end if

  DEALLOCATE (nstiff)
  DEALLOCATE (stiff)
  DEALLOCATE (x,z_all,xi,yt,dydt,yt_scale,tau,tau0,slip,slipinc)
  DEALLOCATE (cca,ccb,xLf,seff)
  call MPI_finalize(ierr)
END program main

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine rkqs(myid,y,dydx,n,Nt_all,Nt,x,htry,eps,yscal,hdid,hnext,z_all,p)
  Use mpi
  USE phy3d_module_non, only : nprocs
  implicit none
  integer, parameter :: DP = kind(1.0d0)
  integer :: n,i,j,k,NMAX,Nt,Nt_all
  real (DP) :: eps,hdid,hnext,htry,x
  real (DP) :: dydx(n),y(n),yscal(n),z_all(Nt_all),p(Nt) !p is position
  external derivs
  real (DP) :: errmax,errmax1,h,htemp,xnew,errmax_all(nprocs)
  real (DP), dimension(:), allocatable :: yerr,ytemp
  real (DP), parameter :: SAFETY=0.9, PGROW=-.2,PSHRNK=-.25,ERRCON=1.89e-4

  !MPI RELATED DEFINITIONS
  integer :: ierr,myid,master
  master = 0

  nmax=3*Nt
  h=htry
  allocate (yerr(nmax),ytemp(nmax))
  !        write(*,*)'in rkqs, bf rkck, myid=',myid,y(n-1)
1 call rkck(myid,dydx,h,n,Nt_all,Nt,y,yerr,ytemp,x,derivs,z_all,p)
  !        write(*,*)'in rkqs, af rkck, myid=',myid, ytemp(n-1)
  errmax=0.
  do  i=1,nmax
     !j = int(ceiling(real(i)/2)) ! position within central part
     !   if(p(j).gt.-700.and.p(j).lt.180)then
     errmax = max(errmax,dabs(yerr(i)/yscal(i)))
     !   end if
  end do
  errmax=errmax/eps
!  write(*,*) errmax

  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  call MPI_Gather(errmax,1,MPI_Real8,errmax_all,1,MPI_Real8,master,MPI_COMM_WORLD,ierr)

  if(myid==master)then
     errmax1=maxval(errmax_all)
  end if
  CALL MPI_BCAST(errmax1,1,MPI_REAL8,master,MPI_COMM_WORLD, ierr)

  if(errmax1.gt.1.)then
     htemp = SAFETY*h*(errmax1**PSHRNK)
     h = dsign(max(dabs(htemp),0.1*dabs(h)),h)
     xnew = x+h
     if(xnew.eq.x)pause 'stepsize underflow in rkqs'
     goto 1
  else
     if(errmax1.gt.ERRCON)then
        hnext=SAFETY*h*(errmax1**PGROW)
     else
        hnext=5.*h
     end if
     hdid=h
     x=x+h
     do i=1,nmax
        y(i)=ytemp(i)
     end do
  end if

  deallocate (yerr,ytemp)

  RETURN
end subroutine rkqs
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
     subroutine rkck(myid,dydx,h,n,Nt_all,Nt,y,yerr,yout,x,derivs,z_all,p)
       USE phy3d_module_non, only :nprocs
       implicit none
       integer, parameter :: DP = kind(1.0d0)
       integer :: n,i,NMAX,myid,Nt_all,Nt
       external derivs
       real (DP) :: h,x,dydx(n),y(n),yerr(n),yout(n),z_all(Nt_all),p(Nt)
       real (DP), dimension(:), ALLOCATABLE :: ak2,ak3,ak4,ak5,ak6,ytemp
       REAL (DP),  parameter :: A2=.2,A3=.3,A4=.6,A5=1.,A6=.875, &
            B21=.2,B31=3./40.,B32=9./40.,B41=.3,&
            B42=-.9,B43=1.2,B51=-11./54.,B52=2.5, &
            B53=-70./27.,B54=35./27., B61=1631./55296., &
            B62=175./512.,B63=575./13824.,B64=44275./110592., &
            B65=253./4096.,C1=37./378., C3=250./621.,  &
            C4=125./594.,C6=512./1771.,DC1=C1-2825./27648., &
            DC3=C3-18575./48384.,DC4=C4-13525./55296.,  &
            DC5=-277./14336.,DC6=C6-.25


       nmax = 3*Nt
       ALLOCATE (ak2(nmax),ak3(NMAX),ak4(NMAX),ak5(NMAX),ak6(NMAX),ytemp(NMAX))

       do  i=1,n
          ytemp(i)=y(i)+B21*h*dydx(i)
       end do
       call derivs(myid,ak2,n,Nt_all,Nt,x+A2*h,ytemp,z_all,p)
       do  i=1,n
          ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
       end do
       call derivs(myid,ak3,n,Nt_all,Nt,x+A3*h,ytemp,z_all,p)
       do  i=1,n
          ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
       end do
       call derivs(myid,ak4,n,Nt_all,Nt,x+A4*h,ytemp,z_all,p)
       do  i=1,n
          ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+  &
               B54*ak4(i))
       end do
       call derivs(myid,ak5,n,Nt_all,Nt,x+A5*h,ytemp,z_all,p)
       do  i=1,n
          ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+  &
               B64*ak4(i)+B65*ak5(i))
       end do
       call derivs(myid,ak6,n,Nt_all,Nt,x+A6*h,ytemp,z_all,p)
       do  i=1,n
          yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+   &
               C6*ak6(i))
       end do
       do  i=1,n
          yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+     &
               DC5*ak5(i)+DC6*ak6(i))
       end do
       DEALLOCATE (ak2,ak3,ak4,ak5,ak6,ytemp)
       return
     end subroutine rkck
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
     subroutine derivs(myid,dydt,nv,Nt_all,Nt,t,yt,z_all,x)
       USE mpi
       USE phy3d_module_non, only: stiff,nstiff,cca,ccb,seff,xLf,eta,f0,Vpl,V0,Lratio,nprocs,&
            tm1,tm2,tmday,tmelse,tmmidn,tmmult
       implicit none
       integer, parameter :: DP = kind(1.0d0)
       integer :: nv,n,i,j,k,kk,l,ii,Nt,Nt_all
       real (DP) :: t,yt(nv),dydt(nv)
       real (DP) :: deriv3,deriv2,deriv1,small
       real (DP) :: psi,help1,help2,help
       real (DP) :: SECNDS
       real (DP) :: z_all(Nt_all),x(Nt),zz(Nt),zzfric(Nt),zz_all(Nt_all)
       intrinsic imag,real
   !ALu
       real (DP) :: deriv4
       real (DP) :: nzzfric(Nt)

       !MPI RELATED DEFINITIONS
       integer :: ierr,myid,master
       master = 0


       small=1.d-3

       do i=1,Nt
          zzfric(i) = 0.0
          nzzfric(i)=0.0
          zz(i)=yt(3*i-1)-Vpl
       end do


       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_Gather(zz,Nt,MPI_Real8,zz_all,Nt,MPI_Real8,master,MPI_COMM_WORLD,ierr)
       call MPI_Bcast(zz_all,Nt_all,MPI_Real8,master,MPI_COMM_WORLD,ierr)

       !----------------------------------------------------------------------
       !    summation of stiffness of all elements in slab
       !----------------------------------------------------------------------
       ! initilize zzfric
       call CPU_TIME(tm2)
       tmelse=tmelse+tm2-tm1
       tm1=tm2

       ! calculate stiffness*(Vkl-Vpl) of one proccessor.
       do i=1,Nt
          do j=1, Nt_all
             ! number in all element, total Nt_all
             zzfric(i) = zzfric(i) + stiff(i,j)*zz_all(j)
          end do
       end do

       do i=1,Nt
          do j=1, Nt_all
             ! number in all element, total Nt_all
             nzzfric(i) = nzzfric(i) + nstiff(i,j)*zz_all(j)
          end do
       end do

 !      do i=1,1
 !      write(*,'(1x,I6,E23.13)') myid*Nt+1,zzfric(1)
 !      end do

       call CPU_TIME(tm2)
       !if ((tm2-tm1) .lt. 0.03)then
          tmmult=tmmult+tm2-tm1
       !else
       !   tmmidn=tmmidn+1
       !   tmmult=tmmult+tm2-tm1+tmday
       !end if
       tm1=tm2

       do i=1,Nt
          !if(yt(3*i).lt.0.00001)then
          !   yt(3*i) = yt(3*i)+0.00001
          !end if

         if(yt(3*i-1).le.small)then
             psi = dlog(V0*yt(3*i)/xLf(i))
             help1 = yt(3*i-1)/(2*V0)
             help2 = (f0+ccb(i)*psi)/cca(i)
             help = dsqrt(1+(help1*dexp(help2))**2)

             deriv1 = (yt(3*i-2)*ccb(i)/yt(3*i))*help1*dexp(help2)/help
             deriv2 = (yt(3*i-2)*cca(i)/(2*V0))*dexp(help2)/help
             deriv3 = 1 - yt(3*i-1)*yt(3*i)/xLf(i)
             deriv4=cca(i)*dasinh(help1*dexp(help2))
             dydt(3*i-2)=-nzzfric(i)
             dydt(3*i-1)=-(zzfric(i)+deriv1*deriv3+deriv4*dydt(3*i-2))/(eta+deriv2)
             dydt(3*i)=deriv3
          else
             deriv1 = yt(3*i-2)*ccb(i)/yt(3*i)
             deriv2 = yt(3*i-2)*cca(i)/yt(3*i-1)
             deriv3 = 1 - yt(3*i-1)*yt(3*i)/xLf(i)
             deriv4=f0+cca(i)*dlog(yt(3*i-1)/v0)+ccb(i)*dlog(V0*yt(3*i)/xLf(i))
             dydt(3*i-2)=-nzzfric(i)
             dydt(3*i-1) = -(zzfric(i)+deriv1*deriv3+deriv4*dydt(3*i-2))/(eta+deriv2)
             dydt(3*i) = deriv3
          end if
       end do

!       if(myid==0) write(*,*) t,dydt(1),dydt(2)

       RETURN
     END subroutine derivs

!-----------------------------------------------------------------------------
!    read parameters: sigma_effective, a,b,D_c
!----------------------------------------------------------------------------

    subroutine resdep(Nt_all,dip,hnucl,sigmadiff,sefffix,Lffix,Iperb,&
         factor1,factor2,factor3,factor4, &
         xilock1,xilock2,cca_all,ccb_all,xLf_all, &
         seff_all,x_all,z_all)
      USE mpi
      USE phy3d_module_non, only: p18,Nab,xmu,xnu,gamma, &
           Iprofile,foldername,jobname,profile
      implicit none
      integer, parameter :: DP = kind(1.0d0)
      integer, parameter :: DN=9
      integer :: k,i,j,kk,Iperb,record,l,m,nn,Nt,Nt_all

      real (DP) :: temp(DN),dep(DN),dist(DN),ptemp(Nt_all), &
           ccabmin(Nt_all),xLfmin(Nt_all),xilock1,xilock2, &
           hnucl,sigmadiff,sefffix,Lffix,dip, &
           factor,factor1,factor2,factor3,factor4
      real (DP) :: cca_all(Nt_all),ccb_all(Nt_all),ccab_all(Nt_all), &
           xLf_all(Nt_all),seff_all(Nt_all),x_all(Nt_all),z_all(Nt_all)

      real (DP) ::a(Nab),tpr(Nab),zp(Nab),b(nab),ab(nab)
      open(447,file="ab.txt",status='old')
       do i=1,Nt_all
        read(447,*) ccab_all(i)
       end do
      close(447)

      open(444,file="sigma.txt",status='old')
       do i=1,Nt_all
        read(444,*) seff_all(i)
       end do
      close(444)

      open(445,file="xLf.txt",status='old')
       do i=1,Nt_all
        read(445,*) xLf_all(i)
       end do
      close(445)

 
      open(446,file="a.txt",status='old')
       do i=1,Nt_all
        read(446,*) cca_all(i)
       end do
      close(446)

      
      open(448,file="b.txt",status='old')
       do i=1,Nt_all
        read(448,*) ccb_all(i)
       end do
      close(447)

      !     To save info about some of the quantities
      open(2,file=trim(foldername)//'vardep'//jobname,status='unknown')
      !	write(2,300)'z','seff','Lf','ccab','cca'
      do i=1,Nt_all
         write(2,'(5(1x,e20.13))')z_all(i),seff_all(i),xLf_all(i), &
              ccab_all(i),cca_all(i)
      end do
      close(2)
300   format(5(1x,A20))
      RETURN
    END subroutine resdep
!
!------------------------------------------------------------------------------
! restart file
!------------------------------------------------------------------------------

      subroutine restart(inout,filename,Ifileout,Nt_all,t,dt,dt_try,ndt,nrec,yt,slip)
USE phy3d_module_non, ONLY : jobname,foldername,restartname, &
                        tm1,tm2,tmday,tmelse,tmmidn,tmmult
      implicit none
      integer, parameter :: DP = kind(1.0d0)
      integer :: inout,i,ndt,nrec,Ifileout,Nt,Nt_all
      real (DP) :: t,dt,dt_try
      real (DP) ::  yt(3*Nt_all),slip(Nt_all)
character(len=40) :: filename

      if(inout.eq.0) then
         open(Ifileout,file=trim(restartname),status='old')
          read(Ifileout,*)t,ndt,nrec
          read(Ifileout,*)dt,dt_try
          do i=1,3*Nt_all
             read(Ifileout,*)yt(i)
          end do

          do i=1,Nt_all
             read(Ifileout,*)slip(i)
          end do
	  close(Ifileout)
      else
         open(Ifileout,file=trim(foldername)//trim(filename)//jobname,status='unknown')
         write(Ifileout,*)t,ndt,nrec
         write(Ifileout,*)dt,dt_try
         do i=1,3*Nt_all
              write(Ifileout,*)yt(i)
         end do
         do i=1,Nt_all
              write(Ifileout,*)slip(i)
         end do
         close(Ifileout)
        end if

      RETURN
      END

!------Output -------------------------------------------
!--------------------------------------------------------
!subroutine output(Ioutput,Isnapshot,Nt_all,Nt,inul,imv,ias,icos,isse,x,&
!    tmv,tas,tcos,tsse,maxv,maxv_s1,maxv_s2,maxv_s3,maxv_s4,maxv_s5,maxv_s6,&
!    maxv_s7,maxv_s8,maxnum, &
!     slipz1_inter,slipz1_tau,slipz1_sse,&
!     slipz1_cos,&
!     xi_all,x_all,slipz1_v)

!ED2
!subroutine output(Ioutput,Isnapshot,Nt_all,Nt,inul,imv,ias,icos,isse,islip,x,&
!    tmv,tas,tcos,tsse,tslip,maxv,maxv_s1,maxv_s2,maxv_s3,maxv_s4,maxv_s5,maxv_s6,&
!    maxv_s7,maxv_s8,maxnum, &
!    slipz1_inter,slipz1_tau,slipz1_sse,&
!    slipz1_cos,slip_tau,slip_v, slip_slip, &
!    xi_all,x_all,slipz1_v,&
!    slipout, sseout, cosout, asout,&
!    itau, ttau,mtau, mtau_s1,mtau_s2,mtau_s3,mtau_s4,mtau_s5,&
!    mtau_s6,mtau_s7,mtau_s8,mtaunum,&
!    issere,tssere,re_slip,re_tau,re_v,nrecenter,reout)

!ED3 !ED4
subroutine output(Ioutput,Isnapshot,Nt_all,Nt,inul,imv,ias,icos,isse,islip,x,&
   tmv,tas,tcos,tsse,tslip,maxv,maxv_s1,maxv_s2,maxv_s3,maxv_s4,maxv_s5,maxv_s6,&
   maxv_s7,maxv_s8,maxnum, &
   slipz1_inter,slipz1_tau,slipz1_sse,slipz1_cos,&
   slip_tau,slip_v,slip_slip,xi_all,x_all,slipz1_v,&
   slipout, sseout, cosout, asout, &
   itau, ttau,mtau, mtau_s1,mtau_s2,mtau_s3,mtau_s4,mtau_s5,mtau_s6,mtau_s7,mtau_s8,mtaunum,&
   issere,tssere,re_slip,re_tau,re_v,nrecenter,reout, re_seff, & 
   maxnum_s, maxnum_n, maxv_s, maxv_n, avev_s, avev_n, aveslip_s, aveslip_n,&
   mtaunum_s, mtaunum_n, mtau_s, mtau_n, avetau_s, avetau_n,&
   maxnum_eq, maxnum_sse, maxv_eq, maxv_sse, avev_eq, avev_sse, aveslip_eq, aveslip_sse,&
   mtaunum_eq, mtaunum_sse, mtau_eq, mtau_sse, avetau_eq, avetau_sse,&
   maxnum_ns, maxv_ns, avev_ns, aveslip_ns, mtaunum_ns, mtau_ns, avetau_ns, &
   iseff, tseff, mseff, mseff_s1, mseff_s2, mseff_s3, mseff_s4, mseff_s5, mseff_s6,&
   mseff_s7, mseff_s8, mseffnum, &
   mseff_s, mseff_n, mseff_ns, mseff_sse, mseff_eq, &
   mseffnum_s, mseffnum_n, mseffnum_ns, mseffnum_sse, mseffnum_eq, &
   aveseff_s, aveseff_n, aveseff_ns, aveseff_sse, aveseff_eq )


USE mpi
!USE phy3d_module_non, only: xmu,nmv,nas,ncos,nsse,yrs,Vpl, &
!		foldername,jobname
!ALu!ED3
USE phy3d_module_non, only: xmu,nmv,nas,ncos,nsse,nslip,yrs,Vpl, &
		foldername,jobname,nssere
implicit none
integer, parameter :: DP = kind(1.0d0)
integer :: Nt,Nt_all,i,j,k,l,kk,inul,imv,ias,icos,isse,Ioutput,Isnapshot
!ED2!ED3!ED5
integer :: islip, slipout, sseout, cosout, asout,itau,issere,nrecenter,reout,iseff
real (DP) :: tslip(nslip),tssere(nssere)
!ED5
real (DP) :: slip_tau(Nt_all,nslip),slip_v(Nt_all,nslip),slip_slip(Nt_all,nslip), &
            re_slip(nrecenter,nssere), re_v(nrecenter, nssere), re_tau(nrecenter,nssere), &
            re_seff(nrecenter,nssere)
character(20) ::        cTemp1
!ED2 ! ED5
real (DP) :: x(Nt),maxnum(nmv),maxv(nmv),maxv_s1(nmv),maxv_s2(nmv),maxv_s3(nmv),maxv_s4(nmv),maxv_s5(nmv),&
        maxv_s6(nmv),maxv_s7(nmv),maxv_s8(nmv), &
	tmv(nmv),tas(nas),tcos(ncos),tsse(nsse), &
   ttau(nmv),mtau_s1(nmv),mtau_s2(nmv),mtau_s3(nmv),mtau_s4(nmv),mtau_s5(nmv),mtau_s6(nmv),mtau_s7(nmv),&
   mtau_s8(nmv),mtaunum(nmv),mtau(nmv), &
   tseff(nmv),mseff_s1(nmv),mseff_s2(nmv),mseff_s3(nmv),mseff_s4(nmv),mseff_s5(nmv),mseff_s6(nmv),mseff_s7(nmv),&
   mseff_s8(nmv),mseffnum(nmv),mseff(nmv)

!ED4 !ED5
real(DP) :: avev_s(nmv),avev_n(nmv),maxv_s(nmv),maxv_n(nmv),maxnum_s(nmv),maxnum_n(nmv),&
   avetau_s(nmv),avetau_n(nmv),mtau_s(nmv),mtau_n(nmv),mtaunum_s(nmv),mtaunum_n(nmv),&
   aveslip_s(nmv),aveslip_n(nmv),&
   avev_eq(nmv),avev_sse(nmv),maxv_eq(nmv),maxv_sse(nmv),maxnum_eq(nmv),maxnum_sse(nmv),&
   avetau_eq(nmv),avetau_sse(nmv),mtau_eq(nmv),mtau_sse(nmv),mtaunum_eq(nmv),mtaunum_sse(nmv),&
   aveslip_eq(nmv),aveslip_sse(nmv),&
   avev_ns(nmv),maxv_ns(nmv),maxnum_ns(nmv),&
   avetau_ns(nmv),mtau_ns(nmv),mtaunum_ns(nmv),&
   aveslip_ns(nmv), &
   aveseff_s(nmv),aveseff_n(nmv),aveseff_ns(nmv),aveseff_sse(nmv),aveseff_eq(nmv), &
   mseff_s(nmv),mseff_n(nmv),mseff_ns(nmv),mseff_sse(nmv),mseff_eq(nmv), &
   mseffnum_s(nmv),mseffnum_n(nmv),mseffnum_ns(nmv),mseffnum_sse(nmv),mseffnum_eq(nmv) 

real (DP) :: slipz1_inter(Nt_all,nas),slipz1_cos(Nt_all,ncos),&
        slipz1_tau(Nt_all,nsse),slipz1_sse(Nt_all,nsse), &
     xi_all(Nt_all),x_all(Nt_all),&
      slipz1_v(Nt_all,ncos)

! integer :: n_intz1,n_intz2,n_intz3,n_cosz1,n_cosz2,n_cosz3
!integer :: intdepz1(Nt_all),intdepz2(Nt_all),intdepz3(Nt_all)



if(Ioutput == 0)then    !output during run


   if(imv==nmv)then
      open(30,file=trim(foldername)//'maxv_all'//jobname,access='append',status='unknown')
      open(311,file=trim(foldername)//'maxv_s1'//jobname,access='append',status='unknown')
      open(312,file=trim(foldername)//'maxv_s2'//jobname,access='append',status='unknown')
      open(313,file=trim(foldername)//'maxv_s3'//jobname,access='append',status='unknown')
      open(314,file=trim(foldername)//'maxv_s4'//jobname,access='append',status='unknown')
      open(315,file=trim(foldername)//'maxv_s5'//jobname,access='append',status='unknown')
      open(316,file=trim(foldername)//'maxv_s6'//jobname,access='append',status='unknown')
      open(317,file=trim(foldername)//'maxv_s7'//jobname,access='append',status='unknown')
      open(318,file=trim(foldername)//'maxv_s8'//jobname,access='append',status='unknown')
      open(319,file=trim(foldername)//'maxnum'//jobname,access='append',status='unknown')
      open(320,file=trim(foldername)//'maxnum_south'//jobname,access='append',status='unknown')
      open(321,file=trim(foldername)//'maxnum_north'//jobname,access='append',status='unknown')
      open(322,file=trim(foldername)//'maxv_s'//jobname,access='append',status='unknown')
      open(323,file=trim(foldername)//'maxv_n'//jobname,access='append',status='unknown')
      open(324,file=trim(foldername)//'avev_s'//jobname,access='append',status='unknown')
      open(325,file=trim(foldername)//'avev_n'//jobname,access='append',status='unknown')
      open(326,file=trim(foldername)//'aveslip_s'//jobname,access='append',status='unknown')
      open(327,file=trim(foldername)//'aveslip_n'//jobname,access='append',status='unknown')
      open(328,file=trim(foldername)//'maxnum_eq'//jobname,access='append',status='unknown')
      open(329,file=trim(foldername)//'maxnum_sse'//jobname,access='append',status='unknown')
      open(330,file=trim(foldername)//'maxv_eq'//jobname,access='append',status='unknown')
      open(331,file=trim(foldername)//'maxv_sse'//jobname,access='append',status='unknown')
      open(332,file=trim(foldername)//'avev_eq'//jobname,access='append',status='unknown')
      open(333,file=trim(foldername)//'avev_sse'//jobname,access='append',status='unknown')
      open(334,file=trim(foldername)//'aveslip_eq'//jobname,access='append',status='unknown')
      open(335,file=trim(foldername)//'aveslip_sse'//jobname,access='append',status='unknown') 
      open(336,file=trim(foldername)//'maxnum_ns'//jobname,access='append',status='unknown')  
      open(337,file=trim(foldername)//'maxv_ns'//jobname,access='append',status='unknown')
      open(338,file=trim(foldername)//'avev_ns'//jobname,access='append',status='unknown')
      open(339,file=trim(foldername)//'aveslip_ns'//jobname,access='append',status='unknown')
      open(340,file=trim(foldername)//'aveseff_s'//jobname,access='append',status='unknown')
      open(341,file=trim(foldername)//'aveseff_n'//jobname,access='append',status='unknown')
      open(342,file=trim(foldername)//'aveseff_ns'//jobname,access='append',status='unknown')
      open(343,file=trim(foldername)//'aveseff_sse'//jobname,access='append',status='unknown')
      open(344,file=trim(foldername)//'aveseff_eq'//jobname,access='append',status='unknown')
      open(345,file=trim(foldername)//'mseff_s'//jobname,access='append',status='unknown')
      open(346,file=trim(foldername)//'mseff_n'//jobname,access='append',status='unknown')
      open(347,file=trim(foldername)//'mseff_ns'//jobname,access='append',status='unknown')
      open(348,file=trim(foldername)//'mseff_sse'//jobname,access='append',status='unknown')
      open(349,file=trim(foldername)//'mseff_eq'//jobname,access='append',status='unknown')
      open(350,file=trim(foldername)//'mseffnum_s'//jobname,access='append',status='unknown')
      open(351,file=trim(foldername)//'mseffnum_n'//jobname,access='append',status='unknown')
      open(352,file=trim(foldername)//'mseffnum_ns'//jobname,access='append',status='unknown')
      open(353,file=trim(foldername)//'mseffnum_sse'//jobname,access='append',status='unknown')
      open(354,file=trim(foldername)//'mseffnum_eq'//jobname,access='append',status='unknown')      
      open(355,file=trim(foldername)//'avetau_s'//jobname,access='append',status='unknown')
      open(356,file=trim(foldername)//'avetau_n'//jobname,access='append',status='unknown')
      open(357,file=trim(foldername)//'avetau_ns'//jobname,access='append',status='unknown')
      open(358,file=trim(foldername)//'avetau_sse'//jobname,access='append',status='unknown')
      open(359,file=trim(foldername)//'avetau_eq'//jobname,access='append',status='unknown')
      open(360,file=trim(foldername)//'mtau_s'//jobname,access='append',status='unknown')
      open(361,file=trim(foldername)//'mtau_n'//jobname,access='append',status='unknown')
      open(362,file=trim(foldername)//'mtau_ns'//jobname,access='append',status='unknown')
      open(363,file=trim(foldername)//'mtau_sse'//jobname,access='append',status='unknown')
      open(364,file=trim(foldername)//'mtau_eq'//jobname,access='append',status='unknown')
      open(365,file=trim(foldername)//'mtaunum_s'//jobname,access='append',status='unknown')
      open(366,file=trim(foldername)//'mtaunum_n'//jobname,access='append',status='unknown')
      open(367,file=trim(foldername)//'mtaunum_ns'//jobname,access='append',status='unknown')
      open(368,file=trim(foldername)//'mtaunum_sse'//jobname,access='append',status='unknown')
      open(369,file=trim(foldername)//'mtaunum_eq'//jobname,access='append',status='unknown') 



      do i=1,nmv
         write(30,130)tmv(i),dlog10(maxv(i)*1d-3/yrs)
         write(311,130) dlog10(maxv_s1(i)*1d-3/yrs)
         write(312,130) dlog10(maxv_s2(i)*1d-3/yrs)
         write(313,130) dlog10(maxv_s3(i)*1d-3/yrs)
         write(314,130) dlog10(maxv_s4(i)*1d-3/yrs)
         write(315,130) dlog10(maxv_s5(i)*1d-3/yrs)
         write(316,130) dlog10(maxv_s6(i)*1d-3/yrs)
         write(317,130) dlog10(maxv_s7(i)*1d-3/yrs)
         write(318,130) dlog10(maxv_s8(i)*1d-3/yrs)
         write(319,*) int(maxnum(i))
         write(320,*) int(maxnum_s(i))
         write(321,*) int(maxnum_n(i))
         write(322,*) dlog10(maxv_s(i)*1d-3/yrs)
         write(323,*) dlog10(maxv_n(i)*1d-3/yrs)
         write(324,*) dlog10(avev_s(i)*1d-3/yrs)
         write(325,*) dlog10(avev_n(i)*1d-3/yrs)
         write(326,*) aveslip_s(i)
         write(327,*) aveslip_n(i)
         write(328,*) int(maxnum_eq(i))
         write(329,*) int(maxnum_sse(i))
         write(330,*) dlog10(maxv_eq(i)*1d-3/yrs)
         write(331,*) dlog10(maxv_sse(i)*1d-3/yrs)
         write(332,*) dlog10(avev_eq(i)*1d-3/yrs)
         write(333,*) dlog10(avev_sse(i)*1d-3/yrs)
         write(334,*) aveslip_eq(i)
         write(335,*) aveslip_sse(i)
         write(336,*) int(maxnum_ns(i))  
         write(337,*) dlog10(maxv_ns(i)*1d-3/yrs)
         write(338,*) dlog10(avev_ns(i)*1d-3/yrs)
         write(339,*) aveslip_ns(i)
         write(340,*) aveseff_s(i)
         write(341,*) aveseff_n(i)
         write(342,*) aveseff_ns(i)
         write(343,*) aveseff_sse(i)
         write(344,*) aveseff_eq(i)
         write(345,*) mseff_s(i)
         write(346,*) mseff_n(i)
         write(347,*) mseff_ns(i)
         write(348,*) mseff_sse(i)
         write(349,*) mseff_eq(i)
         write(350,*) mseffnum_s(i)
         write(351,*) mseffnum_n(i)
         write(352,*) mseffnum_ns(i)
         write(353,*) mseffnum_sse(i)
         write(354,*) mseffnum_eq(i)
         write(355,*) avetau_s(i)
         write(356,*) avetau_n(i)
         write(357,*) avetau_ns(i)
         write(358,*) avetau_sse(i)
         write(359,*) avetau_eq(i)
         write(360,*) mtau_s(i)
         write(361,*) mtau_n(i)
         write(362,*) mtau_ns(i)
         write(363,*) mtau_sse(i)
         write(364,*) mtau_eq(i)
         write(365,*) mtaunum_s(i)
         write(366,*) mtaunum_n(i)
         write(367,*) mtaunum_ns(i)
         write(368,*) mtaunum_sse(i)
         write(369,*) mtaunum_eq(i)


      end do
      close(30)
      do i=311,369
         close(i)
      end do
      imv = 0
   end if

   !ED2  
   if(itau==nmv)then
       
      open(200,file=trim(foldername)//'maxtau_all'//jobname,access='append',status='unknown')
      open(201,file=trim(foldername)//'mtau_s1'//jobname,access='append',status='unknown')
      open(202,file=trim(foldername)//'mtau_s2'//jobname,access='append',status='unknown')
      open(203,file=trim(foldername)//'mtau_s3'//jobname,access='append',status='unknown')
      open(204,file=trim(foldername)//'mtau_s4'//jobname,access='append',status='unknown')
      open(205,file=trim(foldername)//'mtau_s5'//jobname,access='append',status='unknown')
      open(206,file=trim(foldername)//'mtau_s6'//jobname,access='append',status='unknown')
      open(207,file=trim(foldername)//'mtau_s7'//jobname,access='append',status='unknown')
      open(208,file=trim(foldername)//'mtau_s8'//jobname,access='append',status='unknown')
      open(209,file=trim(foldername)//'mtaunum'//jobname,access='append',status='unknown')



      do i=1,nmv
         write(200,*) ttau(i), mtau(i)
         write(201,*) mtau_s1(i)
         write(202,*) mtau_s2(i)
         write(203,*) mtau_s3(i)
         write(204,*) mtau_s4(i)
         write(205,*) mtau_s5(i)
         write(206,*) mtau_s6(i)
         write(207,*) mtau_s7(i)
         write(208,*) mtau_s8(i)
         write(209,*) mtaunum(i)

      end do
      do i=200,209
         close(i)
      end do
      itau = 0
   end if

   !ED5
   if(iseff==nmv)then
       
      open(400,file=trim(foldername)//'mseff_all'//jobname,access='append',status='unknown')
      open(401,file=trim(foldername)//'mseff_s1'//jobname,access='append',status='unknown')
      open(402,file=trim(foldername)//'mseff_s2'//jobname,access='append',status='unknown')
      open(403,file=trim(foldername)//'mseff_s3'//jobname,access='append',status='unknown')
      open(404,file=trim(foldername)//'mseff_s4'//jobname,access='append',status='unknown')
      open(405,file=trim(foldername)//'mseff_s5'//jobname,access='append',status='unknown')
      open(406,file=trim(foldername)//'mseff_s6'//jobname,access='append',status='unknown')
      open(407,file=trim(foldername)//'mseff_s7'//jobname,access='append',status='unknown')
      open(408,file=trim(foldername)//'mseff_s8'//jobname,access='append',status='unknown')

      do i=1,nmv
         write(400,*) tseff(i), mseff(i)
         write(401,*) mseff_s1(i)
         write(402,*) mseff_s2(i)
         write(403,*) mseff_s3(i)
         write(404,*) mseff_s4(i)
         write(405,*) mseff_s5(i)
         write(406,*) mseff_s6(i)
         write(407,*) mseff_s7(i)
         write(408,*) mseff_s8(i)
      end do
      do i=400,408
         close(i)
      end do
      iseff = 0
   end if   

   !ED3
   if(issere==nssere)then
      reout=reout+1
      write(ctemp1,*) reout
      write(*,*) ctemp1
      open(53,file=trim(foldername)//'re-v-'//trim(adjustl(cTemp1))//jobname,access='append',status='unknown')
      open(54,file=trim(foldername)//'re-slip-'//trim(adjustl(cTemp1))//jobname,access='append',status='unknown')
      open(56,file=trim(foldername)//'re-t-'//trim(adjustl(cTemp1))//jobname,access='append',status='unknown')
      open(57,file=trim(foldername)//'re-tau-'//trim(adjustl(cTemp1))//jobname,access='append',status='unknown')
      open(58,file=trim(foldername)//'re-seff-'//trim(adjustl(cTemp1))//jobname,access='append',status='unknown')


      do j=1,nssere
         do i=1,nrecenter
            write(53,*) re_v(i,j)
            write(54,*) re_slip(i,j)
            write(57,*) re_tau(i,j)
            write(58,*) re_seff(i,j)
         end do
      end do

      do i=1,nssere
         write(56,*) tssere(i)
      end do

      close(53)
      close(54)
      close(56)
      close(57)
      close(58)

      issere=0
   end if

    if(ias==nas)then
      asout=asout+1
      write(cTemp1,*) asout
      write(*,*) cTemp1
       open(31,file=trim(foldername)//'slipz1-inter'//trim(adjustl(cTemp1))//jobname,access='append',status='unknown')
       open(34,file=trim(foldername)//'t-inter'//trim(adjustl(cTemp1))//jobname,access='append',status='unknown')

       do j=1,nas
          do i=1,Nt_all
             write(31,*) slipz1_inter(i,j)
          end do
       end do
       do i=1,nas
          write(34,*) tas(i)
       end do
       close(31)
       close(34)
       ias = 0
    end if

    if(icos==ncos)then
       cosout=cosout+1
       write(cTemp1,*) cosout
       write(*,*) cTemp1
       open(42,file=trim(foldername)//'coseis-v-'//trim(adjustl(cTemp1))//jobname,access='append',status='unknown')
       open(45,file=trim(foldername)//'coseis-slip-'//trim(adjustl(cTemp1))//jobname,access='append',status='unknown')
       open(48,file=trim(foldername)//'coseis-t-'//trim(adjustl(cTemp1))//jobname,access='append',status='unknown')
!       open(unit=50,file=trim(foldername)//'slipz1-cos.dat',access='append',status='unknown')
!       open(unit=51,file=trim(foldername)//'slipz1-v.dat',access='append',status='unknown')

       do j=1,ncos
          do i=1,Nt_all
             write(42,*) slipz1_v(i,j)
             write(45,*) slipz1_cos(i,j)
!             write(50,*) slipz1_cos(i,j)
!             write(51,*) slipz1_v(i,j)
          end do
       end do

       do i=1,ncos
          write(48,*) tcos(i)
       end do
       close(42)
       close(45)
       close(48)
!       close(50)
!       close(51)
       icos = 0
    end if

   if(isse==nsse)then
      sseout=sseout+1
      write(cTemp1,*) sseout
      write(*,*) cTemp1
      open(25,file=trim(foldername)//'sse-v-'//trim(adjustl(cTemp1))//jobname,access='append',status='unknown')
      open(26,file=trim(foldername)//'sse-tau-'//trim(adjustl(cTemp1))//jobname,access='append',status='unknown')
      open(28,file=trim(foldername)//'sse-t-'//trim(adjustl(cTemp1))//jobname, access='append', status='unknown')
!      open(unit=52,file=trim(foldername)//'slipz1-sse.dat', access='append', status='unknown')
      do j = 1,nsse
         do i=1,Nt_all
            write(25,*) slipz1_sse(i,j)
            write(26,*) slipz1_tau(i,j)
!            write(52,*) slipz1_sse(i,j)
         end do
      end do

      do i=1,nsse
         write(28,*) tsse(i)
      end do

      isse = 0
!      close(52)
      close(25)
      close(26)
      close(28)
   end if

 !ALu
   if(islip==nslip)then
      slipout=slipout+1
      write(cTemp1,*) slipout
      write(*,*) cTemp1

!      open(60,file=trim(foldername)//'tslip'//jobname, access='append', status='unknown')
!      open(61,file=trim(foldername)//'slip_v'//jobname, access='append', status='unknown')
      open(60,file=trim(foldername)//'slip-t-'//trim(adjustl(cTemp1))//jobname,access='append',status='unknown')
      open(61,file=trim(foldername)//'slip-v-'//trim(adjustl(cTemp1))//jobname,access='append',status='unknown')
      open(62,file=trim(foldername)//'slip-tau-'//trim(adjustl(cTemp1))//jobname,access='append',status='unknown')
      open(63,file=trim(foldername)//'slip-slip-'//trim(adjustl(cTemp1))//jobname,access='append',status='unknown')
      do j=1,nslip
         do i=1,Nt_all
            write(61,*) slip_v(i,j)
            write(62,*) slip_tau(i,j)
            write(63,*) slip_slip(i,j)
         end do
      end do

      do i=1,nslip
         write(60,*) tslip(i)
      end do

      islip=0
      close(60)
      close(61)
      close(62)
      close(63)
   end if

else

   if((imv>0).and.(imv<nmv))then
      open(30,file=trim(foldername)//'maxv_all'//jobname,access='append',status='unknown')
      open(311,file=trim(foldername)//'maxv_s1'//jobname,access='append',status='unknown')
      open(312,file=trim(foldername)//'maxv_s2'//jobname,access='append',status='unknown')
      open(313,file=trim(foldername)//'maxv_s3'//jobname,access='append',status='unknown')
      open(314,file=trim(foldername)//'maxv_s4'//jobname,access='append',status='unknown')
      open(315,file=trim(foldername)//'maxv_s5'//jobname,access='append',status='unknown')
      open(316,file=trim(foldername)//'maxv_s6'//jobname,access='append',status='unknown')
      open(317,file=trim(foldername)//'maxv_s7'//jobname,access='append',status='unknown')
      open(318,file=trim(foldername)//'maxv_s8'//jobname,access='append',status='unknown')
 !      open(319,file=trim(foldername)//'maxnum'//jobname,access='append',status='unknown')
     do i=1,imv
         write(30,130)tmv(i),dlog10(maxv(i)*1d-3/yrs)
         write(311,130) dlog10(maxv_s1(i)*1d-3/yrs)
         write(312,130) dlog10(maxv_s2(i)*1d-3/yrs)
         write(313,130) dlog10(maxv_s3(i)*1d-3/yrs)
         write(314,130) dlog10(maxv_s4(i)*1d-3/yrs)
         write(315,130) dlog10(maxv_s5(i)*1d-3/yrs)
         write(316,130) dlog10(maxv_s6(i)*1d-3/yrs)
         write(317,130) dlog10(maxv_s7(i)*1d-3/yrs)
         write(318,130) dlog10(maxv_s8(i)*1d-3/yrs)
 !         write(319,170) maxnum(i)
      end do
      close(30)
      do i=311,315
         close(i)
      end do
      imv = 0
   end if

    if((ias>0).and.(ias<nas))then
       open(31,file=trim(foldername)//'slipz1-inter.bin',form='unformatted',access='append',status='unknown')
       open(34,file=trim(foldername)//'t-inter'//jobname,access='append',status='unknown')

       do j=1,ias
          do i=1,Nt_all
             write(31) slipz1_inter(i,j)
          end do
       end do
       do i=1,ias
          write(34,*) tas(i)
       end do

       close(31)
       close(34)

       ias = 0
     end if


   if(isse<nsse.and.isse>0)then
      open(25,file=trim(foldername)//'slipz1_sse.bin',form='unformatted',access='append',status='unknown')
      open(26,file=trim(foldername)//'slipz1_tau.bin',form='unformatted',access='append',status='unknown')
      open(28,file=trim(foldername)//'t_sse'//jobname, access='append',status='unknown')
      do j=1,isse
         do i=1,Nt_all
            write(25) slipz1_sse(i,j)
            write(26) slipz1_tau(i,j)
         end do
      end do
      do i=1,isse
         write(28,*) tsse(i)
      end do

      isse = 0
      close(25)
      close(26)
      close(28)
  end if


     if((icos>0).and.(icos<ncos))then
       open(45,file=trim(foldername)//'slipz1-cos.bin',form='unformatted',access='append',status='unknown')
       open(42,file=trim(foldername)//'slipz1-v.bin',form='unformatted',access='append',status='unknown')
       open(48,file=trim(foldername)//'t-cos'//jobname,access='append',status='unknown')

       do j=1,icos
          do i=1,Nt_all
             write(45) slipz1_cos(i,j)
             write(42) slipz1_v(i,j)
          end do
       end do
       do i=1,icos
          write(48,*) tcos(i)
       end do
       close(42)
       close(45)
       close(48)
       icos = 0
     end if

!      if((islip>0).and.(islip<nslip))then
!         open(60,file=trim(foldername)//'tslip'//jobname, access='append', status='unknown')
!         open(61,file=trim(foldername)//'slip_v'//jobname, access='append', status='unknown')
!         do j=1,nslip
!            do i=1,Nt_all
!               write(61,*) slipz1_v(i,j)
!            end do
!         end do
!
!         do i=1,nslip
!            write(60,*) tslip(i)
!         end do
!         close(60)
!         close(61)
!      end if
end if


 120    format(E20.13,4X,E20.13,4X,I6)
 130    format(E20.13,2(1X,E13.6))
 140    format(E20.13)
 150    format(E15.8,3(1X,E15.8))
 160    format(E20.13,1x,E20.13)
 500    format(E15.8,1X,E20.13)
 600    format(E15.4,1X,E13.6,1X,E15.8)
 700    format(E13.6)
 900    format(E15.8)
 170    format(F6.0)


RETURN
END subroutine output

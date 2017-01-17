!# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  #
!# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  #
!# Copyright 2011                                                                  #
!# Xevi Biarnés, Fabio Pietrucci, Fabrizio Marinelli, Alessandro Laio.             #
!#                                                                                 # 
!# This file is part of METAGUI 2.1.                                               #
!#                                                                                 #
!# METAGUI is free software: you can redistribute it and/or modify it under the    #
!# terms of the GNU General Public License as published by the Free Software       #
!# Foundation, either version 3 of the License, or (at your option) any later      #
!# version.                                                                        #
!#                                                                                 #
!# METAGUI is distributed in the hope that it will be useful,  but WITHOUT ANY     # 
!#                                                                                 # 
!# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR   #
!# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.     #
!#                                                                                 #
!# You should have received a copy of the GNU General Public License along with    #
!#                                                                                 #
!# METAGUI.  If not, see <http://www.gnu.org/licenses/>.                           #
!#                                                                                 #
!# If you use this plugin for your work, please cite the following paper in your   #
!# publications and communications:                                                #
!#                                                                                 #
!# Xevi Biarnés, Fabio Pietrucci, Fabrizio Marinelli and Alessandro Laio (2011)    #
!# Computer Physics Communications , XXX, XXX-XXX                                  #
!# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  #
!# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  #
!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
!  This program computes the free-energy of the bins in collective-variable 
!  space (clusters of trajectory frames) built by metagui.tcl, following the 
!  weighted-histogram analysis methodology in:
!  Marinelli, Pietrucci, Laio, and Piana, PLoS Comput. Biol. 5, e100045 (2009).
! 
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
subroutine wham_on_microstates(inputfile)

  !
  implicit none
  character*100 inputfile
  integer :: NCV,NH,NCV_ACTIVE, NCL
  real*8 :: T_FILL,T_SMOOTH, KT, DELTA1,G_CORR
  integer :: NT_TOT 
  real*8, allocatable :: TIME(:),CV_T(:,:)
  character*50 :: FILE_HILL(100),HLABEL(100)
  integer, allocatable :: IACTIVE_CVH(:,:),NACTIVE_CVH(:)
  integer, allocatable :: NTH(:),HLABEL_T(:),CL_T(:)
  integer :: NTH_MAX,MIN_OVERLAP,NGRID
  real*8, allocatable :: RANGE_CV(:,:),RANGE_H(:,:,:)
  logical, allocatable :: is_periodic(:)
  real*8, allocatable :: VG(:,:),VG_ERR(:,:),NEFF(:,:),VEFF(:,:),FSHIFT(:),N_EXP(:),MAXVEFF(:)
  real*8, allocatable :: F_C(:),F_CH(:,:),ERR_F_C(:),ERR_F_CH(:,:)
  integer, allocatable :: P_CL(:)
  !
  real*8 ppt,num,den,TR_N_EXP,tmp
  real*8 minGRIDX,minGRIDY,dGX,dGY
  real*8, allocatable :: minGRID(:),maxGRID(:),dG(:)
  real*8, allocatable :: vggrid1(:,:),vggrid2(:,:)
  real*8 sh(2),ds(2),ss(2),sg(2),dds(2),BOX(2)
  real*8  ds2,ww,ww1,ww2,gt,rr,CC,CCbest,errVG
  real*8 dx,dy,V1,V2,V3,V4,maxvggrid1,maxvggrid2
  real*8 :: smin,smax
  real*8, allocatable :: pp(:),num_c(:),den_c(:),neff_h(:)
  real*8, allocatable :: num_c_err(:),den_c_err(:)
  real*8 :: tmp2,ttt(10)
  real*8, allocatable :: nbiased_ch(:,:)
  real*8, allocatable :: vbiased_ch(:,:),vbiased2_ch(:,:)
  real*8, allocatable :: CHI(:,:),VG_ERR_CL(:,:)
  real*8 :: Fm
  !
  integer :: nlines,i,j,ios,it,ih,nass,ic,nact,COUNVEFF 
  integer :: i0,ixmin,ixmax
  character*100 :: cdum,line
  integer :: id,ig,ix,iy,ixp,iyp,ixt,iyt,ixx,iyy,icv,iact
  integer :: nt_h,nhills,nhills1,nhills2
  integer NGRIDX,NGRIDY
  integer, allocatable :: GRIDBINS(:),cv_act(:)
  integer, allocatable :: within_range(:,:),within_range_t(:,:),connected(:,:),connected_t(:,:),visited(:,:)
  integer :: ngood,ntot,nbest,im(2)
  integer :: nbad
  integer, allocatable :: cl_th(:,:)
  !
  logical :: interval
  logical, allocatable :: ok(:,:),ok2(:)
  !
  !........................................................................ title
  write(6,*)
  write(6,*)' ##################################################################'
  write(6,*)' ##                                                              ##'
  write(6,*)' ##  WHAM-BASED CLUSTER ANALYSIS FOR BIAS-EXCHANGE METADYNAMICS  ##'
  write(6,*)' ##                                                              ##'
  write(6,*)' ##################################################################'
  write(6,*)
  !......................................................................... main
  !
  !
  !
  !
  write(6,'(a)')'!!!!!!!!!!!!!!!!!!  I N P U T   T E M P L A T E  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  write(6,'(a)')'KT     2.4                     ! temperature          (units of HILLS files)        !'
  write(6,'(a)')'T_FILL  40000. 40000.          ! fill and smooth time                               !'
  write(6,'(a)')'FILE_HILL  HILLS.A  1  1       ! Hill_file, with one biased variable, the first one !'
  write(6,'(a)')'FILE_HILL  HILLS.B  1  2       ! Hill_file, with one biased variable, the second one!'
  write(6,'(a)')'DELTA  2.                      ! align VG within 2 KT.                              !'
  write(6,'(a)')'NCV_ACTIVE  3                  ! number of active CVs                               !'
  write(6,'(a)')'RANGE_CV  1  0.0  1.0          ! restrict active CV 1 to range 0.0:1.0              !'
  write(6,'(a)')'RANGE_CV  2  0.5  2.0          ! restrict active CV 2 to range 0.5:2.0              !'
  write(6,'(a)')'RANGE_CV  3  -3.14 3.14        ! restrict active CV 2 to range -3.14:3.14           !'
  write(6,'(a)')'PERIODIC  3                    ! active CV 3 is periodic                            !'
  write(6,'(a)')'MIN_OVERLAP  2                 ! retain cluster if at least 2 walkers have NEFF>=1  !'
  write(6,'(a)')'NGRID        50                ! compute VG aligment on a grid of size 50           !'
  write(6,'(a)')'GCORR        10.               ! correlation time of n                              !'
  write(6,'(a)')'TR_N_EXP     10.               ! Threshold on the number of explored clusters       !'
  write(6,'(a)')'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  !

  KT=0.d0
  T_FILL=0.d0
  T_SMOOTH=0.d0
  NH=0  ! number of HILLS files
  NCV_ACTIVE=0

  MIN_OVERLAP=2
  DELTA1=2.d0
  NGRID=50
  G_CORR=10.
  TR_N_EXP=10.
  allocate(IACTIVE_CVH(10,100),NACTIVE_CVH(100))
  open(55,file=inputfile,status='old')
  !
  do 
     read(55,'(a100)',end=666)line
     write(6,'("INPUT: ",a70)') line
     read(line,*) cdum
     if(trim(cdum).eq."KT")            read(line,*) cdum,KT
     if(trim(cdum).eq."NCV_ACTIVE")    read(line,*) cdum,NCV_ACTIVE
     if(trim(cdum).eq."DELTA")         read(line,*) cdum,DELTA1
     if(trim(cdum).eq."T_FILL")        read(line,*) cdum,T_FILL,T_SMOOTH
     if(trim(cdum).eq."MIN_OVERLAP")   read(line,*) cdum,MIN_OVERLAP
     if(trim(cdum).eq."GCORR")         read(line,*) cdum,G_CORR
     if(trim(cdum).eq."TR_N_EXP")      read(line,*) cdum,TR_N_EXP 
     if (trim(cdum)=='NGRID')          read(line,*) cdum,NGRID 
     if (trim(cdum).eq."FILE_HILL") then
        NH=NH+1
        read(line,*) cdum,FILE_HILL(NH),NACTIVE_CVH(NH),IACTIVE_CVH(1:NACTIVE_CVH(NH),NH)
     endif
  enddo
666  continue

  if(KT==0.) stop 'ERROR: provide KT'

  if(NCV_ACTIVE /= 0) then
     ! BIASED case
     NCV=NCV_ACTIVE

     allocate(RANGE_CV(NCV_ACTIVE,2),is_periodic(NCV_ACTIVE))
     allocate(RANGE_H(NH,NCV_ACTIVE,2))
     RANGE_CV(:,1)=-1.d9
     RANGE_CV(:,2)=+1.d9
     is_periodic(:)=.false.
     RANGE_H(:,:,1)=-1.d9
     RANGE_H(:,:,2)=+1.d9

     rewind(55)
     do 
        read(55,'(a100)',end=667)line
        write(6,'("INPUT: ",a70)') line
        read(line,*)cdum
        if (trim(cdum).eq."RANGE_CV")  read(line,*) cdum,j,RANGE_CV(j,1:2)
        if (trim(cdum).eq."PERIODIC") then
           read(line,*) cdum,j
           is_periodic(j)=.true.
        endif
     enddo
667  continue
     close(55)
     !
     ! end read input
     !
     !  Now we read the labels associated to each  HILLS file, and reads the cluster
     !  labels along the trajectory from the LABELS file (written by metagui.tcl).
     !
     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     !
     do ih=1,NH
        open(21,file=FILE_HILL(ih),status='old')
        ! parsing header
        read(21,'(a100)') line
        write(6,*) trim(FILE_HILL(ih)),":",trim(line)
        ! header format: #! ACTIVE 2  1 2  ABC
        read(line,*) cdum,cdum,nact,ttt(1:nact),HLABEL(ih) ! header of hill file
        write(6,*) "nactive=",nact," label=",HLABEL(ih)
        close(21)
     enddo

     ! read cluster labels from file LABELS
     write(6,'(2x,a)') 'reading cluster labels from file LABELS'

     open(20,file='LABELS',status='old')
     NT_TOT=0
     do
        read(20,*,end=668)
        NT_TOT=NT_TOT+1
     enddo
668  continue
     rewind(20)

     allocate(TIME(NT_TOT),CV_T(NCV,NT_TOT),HLABEL_T(NT_TOT),CL_T(NT_TOT),NTH(NH)) ! big allocation
     NTH(:)=0
     do it=1,NT_TOT
        read(20,*) TIME(it),CL_T(it),cdum,cdum,CV_T(1:NCV,it) ! cluster labels on CV frames
        do ih=1,NH
           if (trim(cdum)==trim(HLABEL(ih))) then
              HLABEL_T(it)=ih
              NTH(ih)=NTH(ih)+1
           endif
        enddo
     enddo
     NTH_MAX=maxval(NTH(1:NH),1)
     close(20)
     write(6,'(2x,a,i12,a)') 'read',NT_TOT,' cluster labels from file LABELS'
     nass=0
     do it=1,NT_TOT
        if (CL_T(it).ne.-999) nass=nass+1
     enddo
     write(6,'(2x,a,i9)') 'assigned cluster labels =',nass
     !
     NCL=maxval(CL_T(:),1)
     !
     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     !
     !  we now read the Gaussians from HILLS files, computes the bias
     !  potential along the trajectory averaging in different time windows,
     !  and discards the regions in collective variable space which do not have
     !  a reliable estimate of the bias potential or which are disconnected from 
     !  the global minimum. Files VG_HILLS* are written.
     !
     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     !
     !  
     write(6,*) 
     write(6,*) '::::::::::::::::: bias potentials VG :::::::::::::::::::::::'
     allocate(VG(NTH_MAX,NH))
     allocate(VG_ERR(NTH_MAX,NH))
     VG(:,:)=-999.d0
     VG_ERR(:,:)=1000.d0
     !
     do ih=1,NH ! loop over hills files
        ! collecting data on this bias
        nt_h=NTH(ih) ! number of biased frames
        if (nt_h==0) then
           write(6,'(2x,a)') 'WARNING: no CV frames for this bias'
           cycle ! nothing to do...
        endif
        !
        nact=NACTIVE_CVH(ih)
        if (nact==0) then
           write(6,'(2x,a)') 'unbiased walker: all frames are good'
           VG(:,ih)=0.d0
           VG_ERR(:,ih)=0.d0          
           cycle ! nothing to do...
        endif
        !
        open(20,file=FILE_HILL(ih),status='old')

        nhills=0
        nhills1=0
        do
           read(20,'(A100)',iostat=ios) line
           if (ios/=0) exit
           if (line(1:1).eq."#") cycle
           nhills=nhills+1
           read(line,*)tmp
           if (nhills1==0 .and. tmp>=T_SMOOTH) then
              nhills1=nhills
           endif
        enddo
        rewind(20)
        nhills2=nhills1+(nhills-nhills1)/2

        read(20,'(A100)') line
        interval=.false.   ! now it is able to read the header of the hills file with interval
        if(index(line,"INTERVAL")>0)then
           if(nact .ne. 1)stop 'INTERVAL is allowed only for 1-dimensional biases'
           i0=index(line,"INTERVAL")+8
           read(line(i0:),*)smin,smax
           interval=.true.
        endif

        write(6,*)
        write(6,'(2x,a,a,i7,i8,100i3)') trim(FILE_HILL(ih)), &
             ' :  nhills, CV frames, active CV =', &
             nhills,nt_h,IACTIVE_CVH(1:nact,ih) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !first find the region on which the bias potential is reliable
        !define the grid
        allocate(minGRID(nact),maxGRID(nact),dG(nact),GRIDBINS(nact),cv_act(nact))
        minGRID(1:nact)=0.d0
        maxGRID(1:nact)=0.d0
        dG(1:nact)=1.d0

        GRIDBINS(1:nact)=NGRID
        minGRID(1:nact)=RANGE_CV(IACTIVE_CVH(1:nact,ih),1)
        maxGRID(1:nact)=RANGE_CV(IACTIVE_CVH(1:nact,ih),2)
        dG(1:nact)=(maxGRID(1:nact)-minGRID(1:nact))/(GRIDBINS(1:nact))

        if(nact==1)then
           NGRIDX=GRIDBINS(1)
           NGRIDY=1
           minGRIDX=minGRID(1)
           minGRIDY=0.d0
           dGX=dG(1)
           dGY=1.d0
        else if(nact==2)then
           NGRIDX=GRIDBINS(1)
           NGRIDY=GRIDBINS(2)
           minGRIDX=minGRID(1)
           minGRIDY=minGRID(2)
           dGX=dG(1)
           dGY=dG(2)
        else
           stop 'bias on more than two CVs is not programmed'
        endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !   reading hills and computing the two VG profiles on the grid
        allocate(vggrid1(NGRIDX,NGRIDY),vggrid2(NGRIDX,NGRIDY))
        vggrid1(:,:)=0.d0
        vggrid2(:,:)=0.d0
        sh(:)=0.d0
        ds(:)=1.d0

        write(6,'(a)') "  computing two VG profiles on the grid"
        write(6,'(a,$)') '  '

        i=0
        do while (i<nhills)
           read(20,'(A100)') line
           if (line(1:1).eq."#") cycle
           i=i+1
           if(mod(i,nhills/50)==0) write(6,'(a,$)') '.'
           read(line,*) tmp,sh(1:nact),ds(1:nact),ww 

           if(i<=nhills1)then
              ww1=ww
              ww2=ww
           else if (i>nhills1 .and. i<=nhills2)then
              ww1=ww*(1.-dble(i-nhills1)/dble(nhills2-nhills1))
              ww2=ww
           else
              ww1=0.
              ww2=ww*(1.-dble(i-nhills2)/dble(nhills-nhills2))
           endif
           ! let's do it only for 2D cases maximum
           !
           do iy=1,NGRIDY
              do ix=1,NGRIDX
                 !               sg(1)=minGRIDX+ix*dGX  AAAA
                 !               sg(2)=minGRIDY+iy*dGY  AAAA
                 sg(1)=minGRIDX+(ix-1)*dGX
                 sg(2)=minGRIDY+(iy-1)*dGY
                 dds(1:nact)=(sg(1:nact)-sh(1:nact))
                 BOX(1)=dGX*(NGRIDX)
                 BOX(2)=dGY*(NGRIDY)
                 do j=1,nact
                    if(is_periodic(IACTIVE_CVH(j,ih)))then
                       dds(j)=dds(j)/BOX(j)
                       dds(j)=dds(j)-anint(dds(j))
                       dds(j)=dds(j)*BOX(j)
                    endif
                 enddo
                 ds2=sum((dds(1:nact)/ds(1:nact))**2)/2. 
                 if(ds2>6.) cycle
                 gt=exp(-ds2)
                 vggrid1(ix,iy)=vggrid1(ix,iy)+ww1*gt ! VG in first window
                 vggrid2(ix,iy)=vggrid2(ix,iy)+ww2*gt ! VG in second window
              enddo
           enddo
           if(interval)then
              ixmin=floor((smin-minGRIDX)/dGX)+1
              ixmax=floor((smax-minGRIDX)/dGX)+2
              do ix=1,ixmin-1
                 vggrid1(ix,1)=vggrid1(ixmin,1)
                 vggrid2(ix,1)=vggrid2(ixmin,1)
              enddo
              do ix=ixmax+1,NGRIDX
                 vggrid1(ix,1)=vggrid1(ixmax,1)
                 vggrid2(ix,1)=vggrid2(ixmax,1)
              enddo
           endif
        enddo
        write(6,*)' done'
        close(20)
        maxvggrid1=maxval(vggrid1)
        maxvggrid2=maxval(vggrid2)
        vggrid1(:,:)=vggrid1(:,:)-maxvggrid1       
        vggrid2(:,:)=vggrid2(:,:)-maxvggrid2       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !find the best possible alignment between grid points that are   visited -> within_range -> connected
        write(6,'(a)') "  aligning two VG profiles by"
        write(6,'(a)') "  maximizing # grid points connected to maximum in VG within DELTA*KT range"
        allocate(within_range(NGRIDX,NGRIDY),within_range_t(NGRIDX,NGRIDY))
        allocate(connected(NGRIDX,NGRIDY),connected_t(NGRIDX,NGRIDY))
        allocate(visited(NGRIDX,NGRIDY))
        ! skip from the counting those grid points not visited by any walker
        visited=0
        do it=1,NT_TOT
           if ((HLABEL_T(it)==ih)) then
              ss(1:nact)=CV_T(IACTIVE_CVH(1:nact,ih),it) ! CV of frames to evaluate
              ix=floor((ss(1)-minGRIDX)/dGX)+1
              iy=floor((ss(2)-minGRIDY)/dGY)+1
              if(ix<1 .or. iy<1 .or. ix>NGRIDX .or. iy>NGRIDY) cycle  
              visited(ix,iy)=1
           endif
        enddo
        write(6,'(a,f6.2,a)') "  grid points visited = ",dble(sum(visited)*100)/dble(NGRIDX*NGRIDY)," %"
        ! this should never happen, because this check has been done before in another way
        if(sum(visited)==0) then
           write(6,*) "  useless. no walker frame is biased by this bias."
           deallocate(within_range,within_range_t,connected,connected_t,vggrid1)
           deallocate(vggrid2,visited,minGRID,maxGRID,dG,GRIDBINS,cv_act)
           cycle
        endif
        !
        nbest=0
        CCbest=0
        connected=0
        within_range=0
        write(6,'(a,$)') '  '
        do ic=-75,75
           if(mod(ic,150/50)==0) write(6,'(a,$)') '.'
           CC=dble(ic)*DELTA1*KT/10.d0
           within_range_t=0
           connected_t=0
           do iy=1,NGRIDY
              do ix=1,NGRIDX
                 if(visited(ix,iy)==0)cycle
                 if(abs(vggrid1(ix,iy)-vggrid2(ix,iy)-CC).lt.DELTA1*KT)within_range_t(ix,iy)=1
              enddo
           enddo
           if(sum(within_range_t)==0) cycle
           im=maxloc(((2*vggrid1(:,:)-vggrid2(:,:)-CC)/2+999)*within_range_t(:,:)-999)
           ix=im(1)
           iy=im(2)
           ! scan connected regions, allowing periodicity on the CVs
           ! the maximum is connected, now we try the others, 
           ! checking many times over the neighbors of each bin:
           connected_t(ix,iy)=1
           ! many iterations = many shells of neighbors...
           do it=1,500
              ! scan each bin on the grid 
              do iyt=1,NGRIDY
                 do ixt=1,NGRIDX
                    if (connected_t(ixt,iyt)==0) cycle
                    ! scan the neighbors of each connected bin
                    do ix=ixt-1,ixt+1
                       ixx=ix
                       if(ix==0 .or. ix==NGRIDX+1) then
                          if (is_periodic(IACTIVE_CVH(1,ih))) then
                             ixx=iabs(ix-NGRIDX) 
                          else
                             cycle
                          endif
                       endif
                       do iy=iyt-1,iyt+1
                          iyy=iy
                          if(iy==0 .or. iy==NGRIDY+1) then
                             if(NGRIDY==1)cycle  ! fahimeh, when it is 1D cycle here 
                             if (is_periodic(IACTIVE_CVH(2,ih))) then
                                iyy=iabs(iy-NGRIDY) 
                             else
                                cycle
                             endif
                          endif
                          ! connect the neighbor
                          if(within_range_t(ixx,iyy)==1) connected_t(ixx,iyy)=1
                       enddo
                    enddo
                 enddo
              enddo
           enddo
           !
           ngood=sum(connected_t)
           !write(6,*)"CC=",CC,"withinrange=",sum(within_range_t),"connected=",ngood
           if(ngood .gt. nbest)then
              nbest=ngood
              CCbest=CC
              connected=connected_t
              within_range=within_range_t
           endif
        enddo
        if(sum(connected)==0) then
           write(6,*) ' impossible to align. given up.'
           stop 'ERROR! VG on grid not aligned'
        endif
        write(6,*)' done'
        !      vggrid2(:,:)=vggrid2(:,:)+CCbest
        vggrid1(:,:)=vggrid1(:,:)-sum(vggrid1(:,:)*connected(:,:))/sum(connected)
        vggrid2(:,:)=vggrid2(:,:)-sum(vggrid2(:,:)*connected(:,:))/sum(connected)
        write(6,'(a,f6.2)')"  CCbest = ",CCbest
        write(6,'(a,f6.2,a,f6.2,a)')"  visited grid points within range = ", &
             &  dble(sum(within_range)*100)/dble(sum(visited))," %   connected = ", &
             &  dble(sum(connected)*100)/dble(sum(visited))," %"
        open(30,file="VG_"//(FILE_HILL(ih)),status='unknown')
        do iy=1,NGRIDY
           do ix=1,NGRIDX
              sg(1)=minGRIDX+(ix-1)*dGX
              sg(2)=minGRIDY+(iy-1)*dGY
              write(30,'(10f10.4)')sg(1:nact),vggrid1(ix,iy),vggrid2(ix,iy), &
                   & (((vggrid1(ix,iy)+vggrid2(ix,iy))/2)+999)*visited(ix,iy)-999, &
                   & (((vggrid1(ix,iy)+vggrid2(ix,iy))/2)+999)*within_range(ix,iy)-999, &
                   & (((vggrid1(ix,iy)+vggrid2(ix,iy))/2)+999)*connected(ix,iy)-999
           enddo
           write(30,*) ""
        enddo
        close(30)
        errVG=0.
        do iy=1,NGRIDY
           do ix=1,NGRIDX
              if(connected(ix,iy)==1) errVG=errVG+(vggrid2(ix,iy)-vggrid1(ix,iy))**2
           enddo
        enddo
        errVG=sqrt(errVG/nbest/2)
        vggrid1(:,:)=(vggrid1(:,:)+vggrid2(:,:))/2    !this is the final VG on the grid
        write(6,'(a,f6.2)')"  average error on the connected grid points:",errVG
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !now we evaluate the value of VG on the CVs by interpolation
        ss=0.d0
        j=0
        ngood=0
        write(6,'(a)') "  assigning VG to each frame"
        do it=1,NT_TOT
           if ((HLABEL_T(it)==ih)) then
              j=j+1
              ss(1:nact)=CV_T(IACTIVE_CVH(1:nact,ih),it) ! CV of frames to evaluate
              ix=floor((ss(1)-minGRIDX)/dGX)+1
              iy=floor((ss(2)-minGRIDY)/dGY)+1
              ixp=ix+1
              iyp=min(NGRIDY,iy+1) 
              if(ix<1 .or. iy<1 .or. ixp>NGRIDX .or. iyp>NGRIDY)cycle
              if((connected(ix,iy)==1) .or. (connected(ixp,iy)==1) .or. &
                   &  (connected(ix,iyp)==1) .or. (connected(ixp,iyp)==1) )then
                 ngood=ngood+1
                 dx=(ss(1)-(dble(ix-1)*dGX+minGRIDX))/dGx
                 dy=(ss(2)-(dble(iy-1)*dGY+minGRIDY))/dGy
                 V1=vggrid1(ix,iy)
                 V2=vggrid1(ixp,iy)
                 V3=vggrid1(ix,iyp)
                 V4=vggrid1(ixp,iyp)
                 VG(j,ih)=V1*(1-dx)*(1-dy)+V2*dx*(1-dy)+V3*(1-dx)*dy+V4*dx*dy
                 VG_ERR(j,ih)=errVG
              endif
           endif
        enddo
        if(j.ne.nt_h)stop 'error in assigning VG'
        write(6,'(a,f6.2,a)') '  good_frames =', dble(ngood*100)/dble(nt_h),' %'
        deallocate(within_range,within_range_t,connected,connected_t,vggrid1,vggrid2,visited,minGRID,maxGRID,dG,GRIDBINS,cv_act)
     enddo ! loop over hill files
     !
     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     !
     !  we now compute the free-energy of each microstate and the 
     !  uncertainty on it by combining the estimates from the different walkers 
     !  using weighted histogram analysis. Files FES and FES_ALLW are written.
     !
     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     !

     write(6,*) 
     write(6,*) '::::::::::::::::: free energy of clusters from WHAM ::::::::'
     write(6,*) 
     allocate(NEFF(NCL,NH),VEFF(NCL,NH),FSHIFT(NH),N_EXP(NCL),MAXVEFF(NH))
     allocate(ok(NTH_MAX,NH),cl_th(NTH_MAX,NH))
     allocate(pp(NCL),num_c(NCL),den_c(NCL),ok2(NCL),neff_h(NH))
     allocate(nbiased_ch(NCL,NH),vbiased_ch(NCL,NH),vbiased2_ch(NCL,NH))
     allocate(num_c_err(NCL),den_c_err(NCL))
     allocate(CHI(NCL,NH),VG_ERR_CL(NCL,NH))
     ! putting cluster labels in a matrix
     cl_th(:,:)=-999
     do ih=1,NH
        i=0
        do it=1,NT_TOT
           if (HLABEL_T(it)==ih) then
              i=i+1
              cl_th(i,ih)=CL_T(it)
           endif
        enddo
     enddo
     !.............................................................. VEFF NEFF
     write(6,'(2x,a)') 'computing VEFF and NEFF'
     ! discarding bad frames
     ok(:,:)=.true.
     nbad=0
     do ih=1,NH
        do it=1,NTH(ih)
           if (VG(it,ih)<-998.) ok(it,ih)=.false.
           if (cl_th(it,ih).eq.-999) ok(it,ih)=.false.
           if (.not.ok(it,ih)) nbad=nbad+1
        enddo
     enddo
     write(6,'(4x,i10,a)') nbad,' frames discarded'
     write(6,'(4x,i10,a)') sum(NTH(1:NH))-nbad,' frames accepted'
     ! compute VEFF and NEFF
     nbiased_ch(:,:)=0.d0  ! fahimeh
     vbiased_ch(:,:)=0.d0
     vbiased2_ch(:,:)=0.d0 ! fahimeh
     do ih=1,NH ! loop over biases
        ! VEFF
        num_c(:)=0.
        den_c(:)=0.
        num_c_err(:)=0.
        den_c_err(:)=0.
        COUNVEFF=0
        do it=1,NTH(ih)
           if (.not.ok(it,ih)) cycle
           ic=cl_th(it,ih)
           num_c(ic)=num_c(ic)+exp(2.*VG(it,ih)/KT)
           den_c(ic)=den_c(ic)+exp(VG(it,ih)/KT)
           num_c_err(ic)=num_c_err(ic)+VG_ERR(it,ih)
           den_c_err(ic)=den_c_err(ic)+1.d0

        enddo
        do ic=1,NCL
           if (num_c(ic)==0..or.den_c(ic)==0.) then 
              VEFF(ic,ih)=-999.
              VG_ERR_CL(ic,ih)=1.d6
           else
              VEFF(ic,ih)=KT*log(num_c(ic)/den_c(ic))
              COUNVEFF=COUNVEFF+1
              if(COUNVEFF.eq.1) then
                 MAXVEFF(ih)=VEFF(ic,ih)
              endif
              if(VEFF(ic,ih).gt.MAXVEFF(ih)) MAXVEFF(ih)=VEFF(ic,ih)
              VG_ERR_CL(ic,ih)=num_c_err(ic)/den_c_err(ic) ! average error on VG(ih) of the cluster (ic)
           endif
        enddo
        ! NEFF
        NEFF(:,ih)=0.
        do it=1,NTH(ih)  
           if (.not.ok(it,ih)) cycle
           ic=cl_th(it,ih)
           if (VEFF(ic,ih)<-998.) cycle
           NEFF(ic,ih)=NEFF(ic,ih)+exp((VG(it,ih)-VEFF(ic,ih))/KT)
           nbiased_ch(ic,ih)=nbiased_ch(ic,ih)+1.
           vbiased_ch(ic,ih)=vbiased_ch(ic,ih)+VG(it,ih)
           vbiased2_ch(ic,ih)=vbiased2_ch(ic,ih)+VG(it,ih)**2
        enddo
        do ic=1,NCL
           CHI(ic,ih)=1.d0/(1.d0+NEFF(ic,ih)/G_CORR*(VG_ERR_CL(ic,ih)/KT))
        enddo
     enddo ! loop on walkers
     !.......................................... cleaning from dirty clusters
     write(6,'(2x,a)') 'cleaning from dirty clusters'
     ok2(1:NCL)=.true.
     nbad=0
     do ic=1,NCL
        j=0
        N_EXP(ic)=0
        do ih=1,NH
           if(VEFF(ic,ih)<-998.) NEFF(ic,ih)=0.
           if (NEFF(ic,ih)<1.) j=j+1  ! important for convergence !!!
           N_EXP(ic)=N_EXP(ic)+NEFF(ic,ih)
        enddo
        if ((NH-j)<MIN_OVERLAP) ok2(ic)=.false.
        if(N_EXP(ic)<TR_N_EXP)  ok2(ic)=.false.
        if (.not.ok2(ic)) nbad=nbad+1
     enddo
     write(6,'(4x,i6,a)') nbad,' clusters discarded (not visited by enough walkers)'
     write(6,'(4x,i6,a)') NCL-nbad,' clusters accepted'
     i=0
     tmp=0.
     tmp2=0.
     do ic=1,NCL
        do ih=1,NH
           if (NEFF(ic,ih)>0.and.nbiased_ch(ic,ih)>0.) then
              vbiased_ch(ic,ih)=vbiased_ch(ic,ih)/nbiased_ch(ic,ih)
              vbiased2_ch(ic,ih)=vbiased2_ch(ic,ih)/nbiased_ch(ic,ih)
              vbiased2_ch(ic,ih)=sqrt(vbiased2_ch(ic,ih)-vbiased_ch(ic,ih)**2)
              tmp=tmp+(VEFF(ic,ih)-vbiased_ch(ic,ih))**2
              tmp2=tmp2+vbiased2_ch(ic,ih)
              i=i+1
           else
              vbiased_ch(ic,ih)=-999.
              vbiased2_ch(ic,ih)=0.
              nbiased_ch(ic,ih)=0.
           endif
        enddo
     enddo
     if (i>0) tmp=sqrt(tmp/dble(i))
     if (i>0) tmp2=tmp2/dble(i)
     write(6,'(2x,a,f10.3)') 'rmsd (VEFF-VG_ave) =',tmp
     write(6,'(2x,a,f10.3)') 'rmsd (VG-VG_ave)   =',tmp2
     write(6,'(4x,a)') '(if they are very different there is a problem)'
     !
     ! for correct normalisation
     do ih=1,NH 
        neff_h(ih)=sum(NEFF(:,ih))
     enddo
     !.................................................... f constants to shift VGs
     write(6,'(2x,a)') 'calculating f constants (shift of VG)'
     write(6,'(4x,a)') 'iteration, f (units of HILLS file) :'
     FSHIFT(:)=0.
     !debug
     !NEFF(:,:)=100.
     !debug
     do it=1,10000
        if (mod(it,1000)==0) write(6,'(i8,100f8.2)') it,FSHIFT(:)
        pp(:)=0.
        do ic=1,NCL
           if (.not.ok2(ic)) cycle
           num=0.
           den=0.
           do ih=1,NH
              if (NEFF(ic,ih)<1. .or. VEFF(ic,ih)<-998.) cycle
              num=num+NEFF(ic,ih)*CHI(ic,ih)
           enddo
           do ih=1,NH
              if (NEFF(ic,ih)<1. .or. VEFF(ic,ih)<-998.) then
                 den=den+neff_h(ih)*exp((FSHIFT(ih)-MAXVEFF(ih))/KT)*CHI(ic,ih)
              else
                 den=den+neff_h(ih)*exp((FSHIFT(ih)-VEFF(ic,ih))/KT)*CHI(ic,ih)
              endif
           enddo
           if(den/=0.)then
              pp(ic)=num/den
           endif
        enddo
        pp(:)=pp(:)/sum(pp(:)) ! normalization of probability
        do ih=1,NH
           ppt=0.
           do ic=1,NCL
              if (.not.ok2(ic)) cycle
              if(NEFF(ic,ih)<1. .or. VEFF(ic,ih)<-998.) then
                 ppt=ppt+exp(-MAXVEFF(ih)/KT)*pp(ic)
              else
                 ppt=ppt+exp(-VEFF(ic,ih)/KT)*pp(ic)
              endif
           enddo
           FSHIFT(ih)=-KT*log(ppt)
        enddo
     enddo
     !.................................................... free energy of clusters
     write(6,'(2x,a)') 'assigning free energy to clusters'
     allocate(F_C(NCL),F_CH(NCL,NH),ERR_F_C(NCL),ERR_F_CH(NCL,NH))
     F_C=1000.
     F_CH=1000.
     ERR_F_C=0.
     ERR_F_CH=0.
     do ic=1,NCL
        if (.not.ok2(ic)) cycle
        F_C(ic)=-KT*log(pp(ic))
        ERR_F_C(ic)=sqrt(G_CORR)*KT/sqrt(sum(NEFF(ic,:)*CHI(ic,:)))
        do ih=1,NH
           F_CH(ic,ih)=-KT*log(NEFF(ic,ih)/neff_h(ih))+FSHIFT(ih)-VEFF(ic,ih)
           ERR_F_CH(ic,ih)=KT/sqrt(NEFF(ic,ih)) ! CHI missing!!!
        enddo
     enddo
     tmp=minval(F_C(:),1)  ! shifting F_C, F_CH to zero of F_C
     do ic=1,NCL
        if (.not.ok2(ic)) cycle
        F_C(ic)=F_C(ic)-tmp
        F_CH(ic,:)=F_CH(ic,:)-tmp
     enddo
     !....................................................... output FES, FES_ALLW
     open(100,file='FES',status='unknown')
     open(101,file='FES_ALLW',status='unknown')
     tmp=-1.d6
     do ic=1,NCL
        ! if(sum(NEFF(ic,:)).eq.0.)cycle
        !      if (F_C(ic)>999.) cycle
        if (F_C(ic)>tmp) tmp=F_C(ic)
        write(100,'(i5,2f12.3)') ic,F_C(ic),ERR_F_C(ic) 
        write(101,'(i5,1000f12.3)') ic,F_C(ic),ERR_F_C(ic), &
             (F_CH(ic,ih),ERR_F_CH(ic,ih),ih=1,NH)
     enddo
     close(100)
     close(101)
     write(6,'(2x,a,f8.2)') 'free-energy dispersion of clusters (units of HILLS files) =',&
          tmp-minval(F_C(:),1)
     write(6,'(2x,a)') 'files FES and FES_ALLW written'
     write(6,'(4x,a)') '(contains: cluster F err)'
     !
     write(6,'(4x,a)') '(contains: cluster S(1:NCV) F F_H(1:NH))'
     deallocate(ok,cl_th,pp,num_c,den_c,ok2,neff_h,nbiased_ch,vbiased_ch,vbiased2_ch,CHI,VG_ERR_CL,num_c_err,den_c_err)


  else
     ! UNBIASED CASE

     ! read cluster labels from file LABELS
     write(6,'(2x,a)') 'unbiased trajectory. Computing the free energy from the histogram'
     write(6,'(2x,a)') 'reading cluster labels from file LABELS'

     open(20,file='LABELS',status='old')
     NT_TOT=0
     do
        read(20,*,end=665)
        NT_TOT=NT_TOT+1
     enddo
665  continue
     rewind(20)

     allocate(TIME(NT_TOT),CL_T(NT_TOT)) 
     do it=1,NT_TOT
        read(20,*) TIME(it),CL_T(it)
     enddo
     close(20)
     write(6,'(2x,a,i12,a)') 'read',NT_TOT,' cluster labels from file LABELS'
     nass=0
     do it=1,NT_TOT
        if (CL_T(it).ne.-999) nass=nass+1
     enddo
     write(6,'(2x,a,i9)') 'assigned cluster labels =',nass
     !
     NCL=maxval(CL_T(:),1)
     allocate(F_C(NCL),ERR_F_C(NCL))
     F_C(:)=0.d0
     do it=1,NT_TOT
        if(TIME(it)>T_FILL .and. CL_T(it)>0)F_C(CL_T(it))=F_C(CL_T(it))+1.d0
     enddo
     open(100,file='FES',status='unknown')
     Fm=9.9D99
     do ic=1,NCL
        ERR_F_C(ic)=KT*sqrt(G_CORR/F_C(ic))
        F_C(ic)=-KT*log(F_C(ic))
        if (F_C(ic).lt.Fm) Fm=F_C(ic)
     enddo
     do ic=1,NCL
        F_C(ic)=F_C(ic)-Fm
        write(100,'(i5,2f12.3)') ic,F_C(ic),ERR_F_C(ic) 
     enddo
     close(100)

  end if
  !

  
end subroutine wham_on_microstates


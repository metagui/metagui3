!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
!#                                                                                 # 
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
!                                                                                  #
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
!                                                                                  # 
! This program computes kinetic basins based on the bins diffusive kinetic model   #
! of ref.                                                                          #
! Marinelli, Pietrucci, Laio, and Piana, PLoS Comput. Biol. 5, e100045 (2009).     #
! Bins in collective-variable space are built by metagui.tcl and the free energies # 
! enetering in the kinetic model are calculated using WHAM (wham_bemeta.f90).      #
! Kinetic basins are build first finding the attractors by analyzing the           #
! rate matrix spectrum (Noe' et al. j.chem.phys. 2007: 126, 155102) and then       #
! assigning the bins to each attractor by a committor analysis.                    #
!                                                                                  #
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#

subroutine kinetic_basins(inputfile)

  implicit none
  ! global variables
  integer NDIM,NC
  real*8, allocatable :: ds(:),ss(:,:),FF(:),KMAT(:,:),prob(:,:)
  integer, allocatable :: nn(:),nnlist(:,:),connected(:),ic_to_i(:),i_to_ic(:)
  integer, allocatable :: old_nn(:),old_nnlist(:,:)
  integer NCconn
  integer ic,jc,idum,it,ir,i,j,icc,jjc,n,ii,jj,icc_min
  real*8 dist,dist_min,r,diag,KT,tauinv,diff,deltaF_max,box,averfes,coefc
  real*8, allocatable :: lambda(:),Z(:,:),WK(:),LAST(:),Zord(:,:),DD(:,:),deltas(:)
  real*8, allocatable :: eigenvec(:,:),disteig(:),diff_CV(:),grid_CV(:),coefa(:) 
  integer, allocatable :: attr(:)
  real*8, allocatable :: base(:,:),vec(:),aver(:),covar(:,:),covarfes(:),committor(:)
  integer, allocatable :: basin(:,:)
  integer NB,b_old,b_new
  logical, allocatable :: periodic(:) 
  logical same_sign,refine
  character*100 inputfile

  integer IER
  integer NVECT
  character*100 tmp
  character*100 CLUSTERS_FES_file
  character*100 evect_file
  refine=.false.
  open (11,file=inputfile)
  read(11,*)CLUSTERS_FES_file
  read(11,*)NDIM,NVECT,KT
  allocate(ds(NDIM),DD(NDIM,NDIM),periodic(NDIM),diff_CV(NDIM))
  allocate(grid_CV(NDIM))
  read(11,*) ds(:)
  read(11,*) grid_CV(:)
  read(11,*) periodic(:)
  do i=1,NDIM
     read(11,*)(DD(i,j),j=1,NDIM)
  end do
  read(11,*) refine
  close (11)


  do i=1,NDIM
     do j=1,NDIM
        DD(i,j)=DD(i,j)/(ds(i)*ds(j))
     enddo
  enddo

  do i=1,NDIM
     do j=1,NDIM
        if(i.ne.j) then
           DD(i,i)=DD(i,i)-dabs(DD(i,j))
        endif
     enddo
  enddo


  deltaF_max=100 !max value of the free energy difference for allowed transitions
  write(*,*) 'deltaF_max=',deltaF_max


  write(6,*) "reading data"

  open(10,file=CLUSTERS_FES_file,status='old')
  NC=0
  do 
     read(10,*,end=111)
     NC=NC+1
  enddo
111 continue
  idum=2*(NDIM**2)
  allocate(ss(NDIM,NC),nn(NC),nnlist(idum,NC),FF(NC),connected(NC),deltas(NDIM))
  allocate(old_nnlist(idum,NC),old_nn(NC),coefa(NDIM))
  allocate(aver(NDIM),covar(NDIM,NDIM),covarfes(NDIM))
  idum=2*NDIM
  rewind(10)
  do ic=1,NC                               
     read(10,*)idum,idum,ss(:,ic),FF(ic)
     ss(:,ic)=ss(:,ic)/ds(:)
  enddo
  close(10)
  
  ! refine free energy by plane fitting neighbours free energies
  if(refine) then 
    call refine_fes(NDIM,NC) 
  endif
  !end refine fes 
  !compute the nnlist
  write(6,*) "reading the connectivity matrix" 
  open(15,file="connectivity",status="old")
  old_nn(:)=0
  old_nnlist(:,:)=0
  do
    read(15,*,end=555)i,j
    if((old_nn(i).gt.2*(NDIM**2)).or.(old_nn(j).gt.2*(NDIM**2))) then
      print*, "Number of neighbours exceeded 2*(NDIM**2)"
      stop
    endif
    old_nn(i)=old_nn(i)+1
    old_nn(j)=old_nn(j)+1
    old_nnlist(old_nn(i),i)=j
    old_nnlist(old_nn(j),j)=i
  enddo
555 continue
  
  do ic=1,NC                                  
     FF(ic)=FF(ic)/KT
  enddo
  FF(:)=FF(:)-minval(FF(:),1)
  do ic=1,NC
     if(old_nn(ic) .eq. 0)cycle
     do jc=1,old_nn(ic)
        jjc=old_nnlist(jc,ic)
        diff=FF(ic)-FF(jjc)
        if(diff>deltaF_max)then
           old_nn(ic)=0
           goto 115
        endif
     enddo
115  continue
  enddo
  do ic=1,NC
     nn(ic)=0
     if(old_nn(ic) .eq. 0)cycle
     do jc=1,old_nn(ic)
        jjc=old_nnlist(jc,ic)
        if(old_nn(jjc)==0)cycle
        diff=FF(jjc)-FF(ic)
        if(diff<deltaF_max) then
           nn(ic)=nn(ic)+1
           nnlist(nn(ic),ic)=jjc
        endif
     enddo
  enddo
  deallocate(old_nn,old_nnlist)
  ! find the clusters that are connected to the minimum in free energy
  write(6,*) "finding the clusters connected to the minimum in free energy"
  ic=minloc(FF(:),1)
  write(*,*) 'location of minimum: ic=',ic
  write(*,*) 'number of beighbors of the minimum: nn(ic)=',nn(ic)
  if(nn(ic)==0)then
    write(*,*) 'ERROR: isolated minimum! You may want to increase deltaF_max...'
    stop
  endif
  connected(:)=0
  connected(ic)=1
  do it=1,2000
    do ic=1,NC
      if(connected(ic)==1)then
        do jc=1,nn(ic)
          connected(nnlist(jc,ic))=1
        enddo
      endif
    enddo
  enddo
  NCconn=sum(connected(:))
  write(6,*) "number of connected clusters: ",NCconn

  open(111,file="CLUSTERS.CONNECTED",status="unknown")
  do ic=1,NC
    if(connected(ic)==1) then
      write(111,'(2i8,100f10.4)') ic,0,ss(1:NDIM,ic)*ds(1:NDIM),FF(ic)*KT
    endif
  enddo
  close(111)
  write(*,*) 'connected clusters written to CLUSTERS.CONNECTED'

  allocate(KMAT(NCconn,NCconn),prob(NCconn,NCconn),ic_to_i(NCconn),i_to_ic(NC))
  ic=0
  do i=1,NC
     if(connected(i)==1)then
        ic=ic+1
        i_to_ic(i)=ic
        ic_to_i(ic)=i
     endif
  enddo
  !now compute the rate matrix
  write(6,*) "computing the rate matrix"
  do ic=1,NCconn
     diag=0.d0
     do ir=1,nn(ic_to_i(ic))
        j=nnlist(ir,ic_to_i(ic))
        ! ic is one cluster in the new list (only connected ones)
        ! ic_to_i(ic) is one cluster in the old list (not only connected ones)
        ! j is its neighbor in the old list
        ! i_to_ic(j) is its neighbor in the new list
        do n=1,NDIM
          diff_CV(n)=(ss(n,ic_to_i(ic))-ss(n,j))
          if (periodic(n)) then
              box=grid_CV(n)
              diff_CV(n)=diff_CV(n)/box
              diff_CV(n)=diff_CV(n)-anint(diff_CV(n))
              diff_CV(n)=diff_CV(n)*box
          endif
        enddo 
        tauinv=0.d0
        do n=1,NDIM
           tauinv=tauinv+DD(n,n)*diff_CV(n)*diff_CV(n)
        enddo
        do ii=1,NDIM-1
           do jj=ii+1,NDIM
              tauinv=tauinv+DD(ii,jj)*diff_CV(ii)*diff_CV(jj)
           enddo
        enddo
        if(tauinv<0.d0) tauinv=0.d0
        if(tauinv>0.d0) then
           KMAT(ic,i_to_ic(j))=tauinv
           diag=diag-tauinv*exp((FF(ic_to_i(ic))-FF(j))/2)
           prob(i_to_ic(j),ic)=tauinv*exp((-FF(j)+FF(ic_to_i(ic)))/2)
        endif
     enddo
     KMAT(ic,ic)=diag
     prob(:,ic)=prob(:,ic)/sum(prob(:,ic))
  enddo

  ! write the rate matrix to file
  write(6,*) "writing the rate matrix and transition probabilities to KINMAT"
  open(997,file="KINMAT",status="UNKNOWN")
  write(997,'(a)') "#  i_conn j_conn  i_old j_old  kinmat_ij  prob_ij"
  write(997,'(a)') "#  where:"
  write(997,'(a)') "#  i_conn j_conn are connected clusters"
  write(997,'(a)') "#  i_old j_old are the original clusters in MICROSTATES"
  write(997,'(a)') "#  kinmat_ij is the symmetrized kinetic matrix, i.e. sqrt(k_ij*k_ji)"
  write(997,'(a)') "#  prob_ij is the normalized (sum_j prob_ij = 1) transition probability i->j"
  do ic=1,NCconn        ! new i
     i=ic_to_i(ic)      ! old i
     ! first write the diagonal element
     write(997,'(2i8,2x,2i8,2x,e20.8)') ic,ic,i,i,KMAT(ic,ic)
     do ir=1,nn(i)
        j=nnlist(ir,i)  ! old j
        jc=i_to_ic(j)   ! new j
        write(997,'(2i8,2x,2i8,2x,e20.8,f20.9)') ic,jc,i,j,KMAT(ic,jc),prob(jc,ic) ! note the swap for prob!
     enddo
  enddo
  close(997)

  write(6,*) "diagonalizing the rate matrix"
  NVECT=NVECT+1
  allocate(lambda(NCconn),Z(NCconn,NCconn))
  call compute_eigenvals(NCconn,KMAT,lambda,Z,NVECT)

  write(6,*) "writing eigenvalues to EIGENVALUES"
  open(998,file="EIGENVALUES",status="UNKNOWN")
  do ic=1,NVECT
     icc=NVECT-ic+1
     write(998,*)ic-1,lambda(icc)
  enddo
  close(998)

  allocate(eigenvec(NCconn,1:NVECT))
  write(6,*) "writing eigenvectors to EIGENVECTOR.*"
  do ic=1,NVECT-1
     icc=NVECT-ic
     eigenvec(:,ic)=Z(:,icc)*Z(:,NVECT) ! eigenvec(:,1) is the eigenvector associated with the slowst mode
  enddo

  do ic=1,NVECT-1
     write(tmp,*) ic
     tmp=trim(adjustl(tmp))
     evect_file="EIGENVECTOR."//tmp
     open(998,file=evect_file,status="UNKNOWN")
     do j=1,NCconn
        i=ic_to_i(j)
        write(998,*) i,eigenvec(j,ic)
     enddo
     close(998)
  enddo

  write(6,*) "writing relaxation times to TAU"
  open(999,file="TAU",status="UNKNOWN")
  do ic=2,NVECT
     icc=NVECT-ic+1
     write(999,*)ic-1,log(1.d0/(-lambda(icc)))
  enddo
  close(999)

  write(tmp,*) NVECT
  tmp=trim(adjustl(tmp))
  call system("(echo set xtics 1 ; echo set grid)> plot_TAU.gpl")
  call system("(echo p [0:"//tmp//"] \'TAU\' w p pt 7 ps 2) >> plot_TAU.gpl")
  call system("gnuplot -persist plot_TAU.gpl")

  
  !now find the attractors
  NVECT=NVECT-1
  allocate(attr(NVECT),disteig(NCconn),base(NVECT,NVECT),vec(NVECT))
  attr(1)=minloc(eigenvec(:,1),1)
  attr(2)=maxloc(eigenvec(:,1),1)
  base(:,1)=eigenvec(attr(2),:)-eigenvec(attr(1),:)
  do ic=2,NVECT-1
    do j=1,NCconn
      do icc=1,ic
        vec(icc)=eigenvec(j,icc)-eigenvec(attr(1),icc)
      enddo
      call compute_dist_from_plane(ic,base(1:ic,1:ic-1),vec(1:ic),disteig(j))
    enddo
    attr(ic+1)=maxloc(disteig(:),1)
    base(:,ic)=eigenvec(attr(ic+1),:)-eigenvec(attr(1),:)
  enddo

  !assign each bin to its attractor
  allocate(basin(NC,NVECT-1),committor(NVECT-1))
  basin(:,:)=0
  committor(:)=0.d0
  
  basin(:,:)=0
  NB=NVECT-1
  do j=1,NCconn
    call compute_committor(j,NB,attr(1:NB),committor(1:NB))
    basin(ic_to_i(j),NB)=attr(maxloc(committor(1:NB),1))
  enddo
   
  do while (NB>2)
    j=attr(NB)
    b_old=basin(ic_to_i(j),NB)
    NB=NB-1
    basin(:,NB)=basin(:,NB+1)
    !compute the committor of j=attr(NB+1)
    call compute_committor(j,NB,attr(1:NB),committor(1:NB))
    b_new=attr(maxloc(committor(1:NB),1)) !this is the attractor to which j is committed
    do j=1,NCconn
      if(basin(ic_to_i(j),NB)==b_old)basin(ic_to_i(j),NB)=b_new !reassigne the state to the new attractor
    enddo
  enddo

  open(60,file='CLUSTERS.BASINS.ALL_EIGENVECTORS',status='unknown')
  do NB=2,NVECT-1
     do  j=1,NC
        if(basin(j,NB)==0)cycle
        do jj=1,NB
           if(attr(jj)==basin(j,NB))basin(j,NB)=jj
           if(attr(jj)==i_to_ic(j))basin(j,NB)=-jj
        enddo
     enddo
  enddo

  do  j=1,NC
    write(60,'(60i8)')j,basin(j,2:)
  enddo
  close(60)
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
contains
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  ! 
  !  This subroutine calculates to which attractor is committed a bin
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine compute_committor(istart,NB,attr,committor)
  implicit none
  integer istart,NB
  integer attr(NB)
  real*8 committor(NB)
  !locals
! integer nt,i,j,state,ib,nn_conn(NCconn),nnlist_conn(NCconn,NCconn),statetry
  integer nt,i,j,state,ib,statetry
  integer, allocatable :: nn_conn(:),nnlist_conn(:,:)
  real*8 rand,thisp

  allocate(nn_conn(NCconn),nnlist_conn(NCconn,NCconn))
  committor(:)=0.d0
  do ib=1,NB
    if(istart==attr(ib))then
       committor(ib)=1.d0
       return
    endif
  enddo

  do ic=1,NCconn
     i=ic_to_i(ic)
     nn_conn(ic)=nn(i)
     do j=1,nn_conn(ic)
        nnlist_conn(j,ic)=i_to_ic(nnlist(j,i))
     enddo
  enddo

  nt=0
  do while (nt<100)
     state=istart
     do
        call random_number(rand)
        thisp=0.d0
        do j = 1,nn_conn(state)
           statetry=nnlist_conn(j,state)
           thisp=thisp+prob(statetry,state)
           if(thisp.gt.rand) then
              state=statetry
              do ib=1,NB
                if(state==attr(ib))then
                   committor(ib)=committor(ib)+1.d0
                   nt=nt+1
                   goto 89
                endif
              enddo
              goto 88
           endif
        enddo
88      continue
     enddo
89   continue
  enddo
  committor(:)=committor(:)/dble(nt)
  deallocate(nn_conn,nnlist_conn)
    
    
  end subroutine compute_committor

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  ! 
  !  This subroutine finds the first NVECT (n_eigen) eigenvalues and eigenvectors
  !  of the rate matrix KMAT (MAT)
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine compute_eigenvals(NDIM,MAT,lambda,eigenvec,n_eigen)
    implicit none
    integer NDIM,n_eigen
    real*8 MAT(NDIM,NDIM),lambda(NDIM),eigenvec(NDIM,NDIM)
    integer LWORK
    real*8, allocatable :: UU(:,:),WORK(:)
    integer, allocatable :: iwork(:),ifail(:)
    integer info,IU,IL,M
    allocate(UU(NDIM,NDIM),WORK(NDIM**2),iwork(5*NDIM),ifail(NDIM))
    LWORK=NDIM**2
    IL=1
    IU=n_eigen
    IL=NDIM-n_eigen+1
    IU=NDIM
    M=IU-IL+1

    !   call DSYEV( 'V', 'U', NDIM, MAT, NDIM, lambda, WORK, LWORK, INFO )
    call DSYEVX( 'V', 'I', 'U', NDIM, MAT, NDIM, 0.d0, 0.d0, IL, IU,  &
         0.d0, M, lambda, eigenvec, NDIM, WORK, LWORK, IWORK,  &
         IFAIL, INFO )

    deallocate(UU,WORK,iwork,ifail)

  end subroutine compute_eigenvals

  subroutine compute_dist_from_plane(ND,base,vec,dist)

  implicit none
  integer ND
  real*8 base(ND,ND-1),vec(ND),dist
  ! locals
  integer i,j
  real*8 vecl(ND),scalp,norm,basel(ND,ND-1)
  !remove from vec the components along the base vectors
  basel(:,1)=base(:,1)/sqrt(sum(base(:,1)*base(:,1)))
  do i=2,ND-1
    basel(:,i)=base(:,i)
    do j=1,i-1
      scalp=sum(basel(:,j)*basel(:,i))
      basel(:,i)=basel(:,i)-basel(:,j)*scalp
    enddo
    basel(:,i)=basel(:,i)/sqrt(sum(basel(:,i)*basel(:,i)))
  enddo

  vecl(:)=vec(:)
  do i=1,ND-1
    scalp=sum(vecl(:)*basel(:,i))
    vecl(:)=vecl(:)-basel(:,i)*scalp
  enddo
  dist=sqrt(sum(vecl(:)**2))
  
  end subroutine compute_dist_from_plane
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  ! 
  ! This subroutine assigns the free energy to bins having FES=1000
  ! by plane fitting neighbouring bin free energies
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX 

  subroutine refine_fes(NDIM,NC)
  integer NDIM,NC
  integer ic,jc,idum,it,ir,i,j,icc,jjc,n
  integer old_nn(NC),old_nnlist(2*(NDIM**2),NC),old_nenn(NC),old_nennlist(2*NDIM,NC) 
  real*8 FF_ref(NC),sigmafes,coefa(NDIM),aver(NDIM),covar(NDIM,NDIM),covarfes(NDIM)
  logical skip
  FF_ref(:)=FF(:)
  old_nn(:)=0
  old_nenn(:)=0 
  do ic=1,NC-1
     do jc=ic+1,NC
        dist=0.d0
        do n=1,NDIM
           diff_CV(n)=ss(n,ic)-ss(n,jc)
           if (periodic(n)) then
              box=grid_CV(n)
              diff_CV(n)=diff_CV(n)/box
              diff_CV(n)=diff_CV(n)-anint(diff_CV(n))
              diff_CV(n)=diff_CV(n)*box
           endif  
           dist=dist+diff_CV(n)**2
        enddo
        if(dist<2.01d0)then
           old_nn(ic)=old_nn(ic)+1
           old_nn(jc)=old_nn(jc)+1
           old_nnlist(old_nn(ic),ic)=jc
           old_nnlist(old_nn(jc),jc)=ic
           if((old_nn(ic).gt.2*(NDIM**2)).or.(old_nn(jc).gt.2*(NDIM**2))) then
             print*, "Number of neighbours exceeded 2*(NDIM**2)"
             stop
           endif
           if(dist<1.01d0)then
              old_nenn(ic)=old_nenn(ic)+1
              old_nenn(jc)=old_nenn(jc)+1
              old_nennlist(old_nenn(ic),ic)=jc 
              old_nennlist(old_nenn(jc),jc)=ic
              if((old_nenn(ic).gt.2*NDIM).or.(old_nenn(jc).gt.2*NDIM)) then
                print*, "Number of near neighbours exceeded 2*NDIM"
                stop
              endif
           endif
        endif
     enddo
  enddo
  do ic=1,NC
     if(FF(ic)>900.0) then
       idum=0
       do jc=1,old_nenn(ic)
          if(FF(old_nennlist(jc,ic))>900.0) cycle
          idum=idum+1
       enddo 
       if(idum>=NDIM+1) then
         aver(:)=0.d0
         covar(:,:)=0.d0
         covarfes(:)=0.d0
         averfes=0.d0
         sigmafes=0.d0 
         it=0
         coefa=0
         coefc=0
         do jc=1,old_nenn(ic)
            if(FF(old_nennlist(jc,ic))>900.0) cycle
            jjc=old_nennlist(jc,ic)
            it=it+1
            averfes=averfes+FF(jjc)
            sigmafes=sigmafes+FF(jjc)*FF(jjc) 
            do i=1,NDIM
              diff_CV(i)=ss(i,jjc)-ss(i,ic)
              if (periodic(i)) then
                 box=grid_CV(i)
                 diff_CV(i)=diff_CV(i)/box
                 diff_CV(i)=diff_CV(i)-anint(diff_CV(i))
                 diff_CV(i)=diff_CV(i)*box
              endif
              aver(i)=aver(i)+diff_CV(i)
              covarfes(i)=covarfes(i)+(FF(jjc)*diff_CV(i))
              covar(i,i)=covar(i,i)+diff_CV(i)*diff_CV(i) 
              if(i<NDIM) then
                do j=i+1,NDIM
                   diff_CV(j)=ss(j,jjc)-ss(j,ic)
                   if (periodic(j)) then
                      box=grid_CV(j)
                      diff_CV(j)=diff_CV(j)/box
                      diff_CV(j)=diff_CV(j)-anint(diff_CV(j))
                      diff_CV(j)=diff_CV(j)*box
                   endif
                   covar(i,j)=covar(i,j)+diff_CV(i)*diff_CV(j)
                enddo
              endif
            enddo
         enddo 
         aver(:)=aver(:)/it
         covarfes(:)=covarfes(:)/it
         averfes=averfes/it
         sigmafes=(sigmafes/it)-(averfes**2)
         sigmafes=sqrt(sigmafes) 
         covarfes(:)=covarfes(:)-averfes*aver(:)
         skip=.false.
         do i=1,NDIM
            covar(i,i)=(covar(i,i)/it)-(aver(i)*aver(i)) 
            if(covar(i,i).eq.0.d0) skip=.true.
         enddo
         do i=1,NDIM-1
            do j=i+1,NDIM
               covar(i,j)=(covar(i,j)/it)-(aver(i)*aver(j))
               covar(j,i)=covar(i,j)
            enddo
         enddo
         if(skip.eqv..false.) then
           call Tline(covar,covarfes,coefa,NDIM)
           do i=1,NDIM
              coefc=coefc-coefa(i)*aver(i) 
           enddo
           coefc=averfes+coefc
           FF_ref(ic)=coefc
           if(FF_ref(ic)<averfes-sigmafes.or.FF_ref(ic)>averfes+sigmafes) FF_ref(ic)=FF(ic)
         endif
       else
         idum=0
         do jc=1,old_nn(ic)
            if(FF(old_nnlist(jc,ic))>900.0) cycle
            idum=idum+1
         enddo
         if(idum>=NDIM+1) then
           aver(:)=0.d0
           covar(:,:)=0.d0
           covarfes(:)=0.d0
           averfes=0.d0
           sigmafes=0.d0
           it=0
           coefa=0
           coefc=0
           do jc=1,old_nn(ic)
              if(FF(old_nnlist(jc,ic))>900.0) cycle
              jjc=old_nnlist(jc,ic)
              it=it+1
              averfes=averfes+FF(jjc)
              sigmafes=sigmafes+FF(jjc)*FF(jjc) 
              do i=1,NDIM
                diff_CV(i)=ss(i,jjc)-ss(i,ic)
                if (periodic(i)) then
                   box=grid_CV(i)
                   diff_CV(i)=diff_CV(i)/box
                   diff_CV(i)=diff_CV(i)-anint(diff_CV(i))
                   diff_CV(i)=diff_CV(i)*box
                endif
                aver(i)=aver(i)+diff_CV(i)
                covarfes(i)=covarfes(i)+(FF(jjc)*diff_CV(i))
                covar(i,i)=covar(i,i)+diff_CV(i)*diff_CV(i)
                if(i<NDIM) then
                  do j=i+1,NDIM
                     diff_CV(j)=ss(j,jjc)-ss(j,ic)
                     if (periodic(j)) then
                        box=grid_CV(j)
                        diff_CV(j)=diff_CV(j)/box
                        diff_CV(j)=diff_CV(j)-anint(diff_CV(j))
                        diff_CV(j)=diff_CV(j)*box
                     endif
                     covar(i,j)=covar(i,j)+diff_CV(i)*diff_CV(j)
                  enddo
                endif
              enddo
           enddo
           aver(:)=aver(:)/it
           covarfes(:)=covarfes(:)/it
           averfes=averfes/it
           sigmafes=(sigmafes/it)-(averfes**2)
           sigmafes=sqrt(sigmafes)
           covarfes(:)=covarfes(:)-averfes*aver(:) 
           skip=.false.
           do i=1,NDIM
              covar(i,i)=(covar(i,i)/it)-(aver(i)*aver(i))
              if(covar(i,i).eq.0.d0) skip=.true.
           enddo
           do i=1,NDIM-1
              do j=i+1,NDIM
                 covar(i,j)=(covar(i,j)/it)-(aver(i)*aver(j))
                 covar(j,i)=covar(i,j)
              enddo
           enddo
           if(skip.eqv..false.) then
             call Tline(covar,covarfes,coefa,NDIM)
             do i=1,NDIM
                coefc=coefc-coefa(i)*aver(i)
             enddo
             coefc=averfes+coefc
             FF_ref(ic)=coefc
             if(FF_ref(ic)<averfes-sigmafes.or.FF_ref(ic)>averfes+sigmafes) FF_ref(ic)=FF(ic)
           endif
         endif
       endif
     endif
  enddo
  FF(:)=FF_ref(:ic)

  end subroutine refine_fes 

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  ! 
  ! This subroutine solve a linear system of equations
  ! 
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine Tline(A,B,X,N)
  INTEGER N,I,J,K
  REAL*8 A(N, N), X(N), B(N)
  REAL*8 S
  
  ! Transform A into triangular matrix
  DO K = 1, N-1
    DO I = K + 1, N
      B(I) = B(I) - A(I, K) / A(K, K) * B(K)
      DO J = K + 1, N
        A(I, J) = A(I, J) - A(I, K) / A(K, K) * A(K, J)
      END DO
    END DO
  END DO

  ! Solve triangular system
  X(N) = B(N) / A(N, N)
  DO I = N - 1, 1, -1
    S = 0.d0
    DO K = I + 1, N
      S = S + A(I, K) * X(K)
    END DO
    X(I) = (B(I) - S) / A(I, I)
  END DO

  ! Print results:

  10 format('  B(',I2,') = ',F6.2)
  11 format('  X(',I2,') = ',F8.4)
  20 format(3(F6.2,'  '))


  end subroutine Tline


end subroutine kinetic_basins

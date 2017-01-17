! 
! program clustering v 1.0
! 
! Program to deal with clustering in METAGUI 
! 
! Nowadays it includes:
! 1) k-medoids (a) clustering with initial membering chosen as in k-means++ (b)
! 2) gromacs like clustering (c) 
! 3) grid clustering
! 4) cluster connectivity to be employed in the free energy surface analysis
! 5) silhouette checking (optative, clustering quality assesment(d))
! 6) Random sieved version (e) for all the algorithms except grid
! 
! TODO list: 
!
! 1.- Improve deallocation (in fact, introduce deallocation of unnecessary
! arrays)
!
! REFERENCES:
!
! (a) Kaufman, L. and Rousseeuw, P.J. (1987), Clustering by means of Medoids, in
! Statistical Data Analysis Based on the L_1–Norm and Related Methods, edited by
! Y. Dodge, North-Holland, 405–416.
! (b) Arthur, D. and Vassilvitskii, S. (2007). "k-means++: the advantages of
! careful seeding". Proceedings of the eighteenth annual ACM-SIAM symposium on
! Discrete algorithms. pp. 1027–1035.
! (c) Daura, X.; Gademann, K.; Jann, B.; Seebach, D.; van Gunsteren, W. F.;
! Mark, A. E. Angew. Chem., Int. Ed. 1999, 38, 236−240.
! (d) Peter J. Rousseeuw (1987). "Silhouettes: a Graphical Aid to the
! Interpretation and Validation of Cluster Analysis". Computational and Applied
! Mathematics 20: 53–65.
! (e)     Clustering Molecular Dynamics Trajectories: 1. Characterizing the
! Performance of Different Clustering Algorithms Jianyin Shao,Stephen W. Tanner,
!  Nephi Thompson, and, and Thomas E.  Cheatham, III Journal of Chemical Theory 
! and Computation 2007 3 (6), 2312-2334
! 
! Alex Rodriguez (alexdepremia@gmail.com)
! 
  subroutine make_microstates(inputfile)
  implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! VARIABLE DECLARATION !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!
  !!! GENERAL VARIABLES !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!
  character*100 inputfile
  integer,parameter:: single= selected_real_kind(7,2)
  integer,parameter:: double= selected_real_kind(15,3)
  integer,parameter:: digit= selected_int_kind(12)
  integer(kind=digit) i,ii,j,jj,k,kk,n,nna,l,nCV,ierr,nframes,ncl,nlines
  !nlines=number of lines in coordinates file
  !ncl= Total number of clusters                       
  !nframes=total number of frames considered
  !ierr=error return
  !nCV=number of collective variables
  !dummy variables: i,ii,j,jj,k,kk,n,nna,l
  real(kind=single) d,dprev,dpost,distance
  !distance=actual value of the weighted euclidean distance between two frames 
  ![The distance between two frames is weighted by the grid employed, in such 
  ! a way that if distance is "1", it means that the distance between two 
  ! frames is one grid hypercube. Thus, the gridding define the metrics]
  ! 
  !d=dummy variable that stores the difference between a given CV
  !dprev & dpost=dummy variables that control the value of d in case of periodicity
  real(kind=single) r,x,y,z
  !r=generated random number
  !x,y,z=dummy real(kind=single) numbers
  integer(kind=digit), dimension(:), allocatable :: ngrid,cluster
  !ngrid(nCV)=number of gridding points per CV
  !cluster(nframes)=cluster assigned to each frame
  integer(kind=digit), dimension(:), allocatable :: ncenter,ncele
  !ncenter(ncl)=element that corresponds to the center of a cluster
  !ncele(ncl)=number of elements that behaves to the cluster
  real(kind=single), dimension(:,:), allocatable :: center
  ! center(ncl,nCV)= coordinates of the center of the cluster
  integer(kind=digit), dimension(:,:), allocatable :: eq
  !eq(nframes,2) = information for each frame readed in COORDINATES and needed to generate MICROSTATES
  character(len=800) :: cdum,line,algorithm
  character(len=500) :: longline
  !cdum, line... character for reading input
  real(kind=single), dimension(:,:), allocatable :: CV
  !CV(nframes,nCV)= CV value for each frame
  real(kind=single), dimension(:), allocatable :: CVdum
  !CVdum(nCV)=dummy variable modified with the aim of accounting of out of
  !limits
  real(kind=single), dimension(:), allocatable :: time
  !time(nframes)= time value per frame
  real(kind=single), dimension(:), allocatable :: grid,period
  !grid(nCV)= size of the grid for each CV
  !period(nCV)= size of the period for each CV
  real(kind=single), dimension(:,:), allocatable :: cgrid
  !cgrid(2,nCV)= minimum and maximum values for the gridding
  logical(1), dimension(:), allocatable :: periodic
  !periodic(nCV)=true if the CV is periodic, false otherwise
  real(kind=single), dimension(:), allocatable :: tdist
  !tdist(nframes)= temporal array for distances (allows parallelism)
  integer(kind=digit) :: ilab
    !!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! SIEVING PARAMETERS !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!
  logical(1) sieved
  !sieved= true if we are performing sieving clustering, false otherwise
  integer(kind=digit) nsieve,nelecl,nmin
  !nmin=cluster at which each frame has minimum distance
  !nelecl=number of elements for clustering
  logical(1), dimension(:), allocatable :: chosen
  !chosen(nframes)=true if the frame has been selected along the random screening for sieving
  integer(kind=digit), dimension(:), allocatable :: chlist
  !chlist(nelecl)=list of chosen elements for sieved clustering
  integer(kind=digit), dimension(:), allocatable :: tempcluster
  !tempcluster(nelecl)=temporaly saves the cluster assigment to undo sieving
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! ALGORITHM SPECIFIC VARIABLES !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!
    !!! k-medoids !!!
    !!!!!!!!!!!!!!!!!

  real(kind=single), dimension (:), allocatable :: distmat
  ! distmat(nframes*(nframes-1)/2)=upper triangle of distance matrix
  real(kind=single), dimension (:), allocatable :: prob
  ! prob(nframes)=probability of a given frame of being chose as center (upper bound)
  !               (k-means++ algorithm)
  real(kind=single) :: dmin, score, minscore
  ! dmin=dummy that allows the computing of minimum distance
  ! score=dummy that allows the computing of the score function
  ! minscore=dummy that allows the computing of minimum scoring member of the cluster
  integer(kind=digit) new_cluster,ndimmat,limit,niter,ndifr,maxiter
  ! ndimmat= dimension of dismat nframes*(nframes-1)/2
  ! new_cluster=new cluster assignation
  ! limit= dummy that allows following the i,j indexing while scanning distmat
  ! niter=counter of the number of iterations
  ! maxiter=maximum number of iterations
  ! ndifr=number of frames with a different cluster asssignation between iterations
  logical finished
  ! finished= controls the iterative algorithm

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! gromacs_like clustering !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer(kind=digit) maxneigh,nmax,nns
  !nsieve=number of elements for sieved clustering
  !nns= dummy that controls the number of assigned neighbours in a neighbour list after each iteration
  !nmax=element with the maximum number of neighbours
  !maxneigh=maximum number of neighbours in the first clustering step (to allocate in memory)
  real(kind=single) cutoff!!!! ,dmin (already declared in k-medoids but with different meaning)
  !dmin=minimum distance between a frame and the cluster center,computed for sieved clustering (not 
  !cutoff=cutoff distance for the generation of the neighbour list
  integer(kind=digit), dimension(:), allocatable :: nneigh,temp
  !nneigh(nframes)=number of neighbours of each frame
  !temp(nframes)=dummy array to store the neighbour list for each frame before write it to file
  integer(kind=digit), dimension(:,:), allocatable :: neighbour
  !neighbour(nframes,maxneigh)=neighbour list for each frame  
  logical(1), dimension(:), allocatable :: assigned
  !assigned(nframes)=true if the frame has already assigned to a cluster, false otherwise

    !!!!!!!!!!!!!!!!!!!!!!!
    !!! grid clustering !!!
    !!!!!!!!!!!!!!!!!!!!!!!
  integer(kind=digit), dimension(:), allocatable :: pos
  !pos(nframes)= gives the position of a frame in the grid
  integer(kind=digit), dimension(:), allocatable :: bincl
  !bincl(ncluster)= gives the equivalence position-cluster
  integer(kind=digit), dimension(:), allocatable :: shift
  !shift(nCV)= help to determine position in a grid
  integer(kind=digit) :: myshift,pp
  !myshift=dummy to compute shift
  !pp=dummy for position
  logical(1) outside
  !outside= control if the CVs are outside the grid
  logical(1), dimension(:), allocatable :: posed
  ! posed(nframes): determine if a given position has been assigned
  real(kind=single) x1,x2
  !dummy for connectivity
  real(kind=single) cent_distance
  ! cent_distance= distance between centers
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! POST-CLUSTERING PROCESS VARIABLES !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!
    !!! silhouette !!!
    !!!!!!!!!!!!!!!!!!
  logical compute_sil
  ! compute_sil= if true computes: 1.-the silhouette for each frame
  !                                2.-the average silhouette for each cluster
  !                                3.-the average silhouette for all the frames
  ! Quite time consumming. Take a look to J. of Computational & Applied Mathematics 20 (1987) 53-65
  !
  real(kind=single), dimension(:), allocatable :: silhouette, average_dissimilarity
  ! silhouette (nframes) = silhouette per each frame
  ! average_dissimilarity (ncl) = average dissimilarity per each frame
  real(kind=single) :: sil_average
  ! sil_average = average value of silhouette
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! CONNECTIVITY MATRIX !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
  real(kind=single),dimension (:),allocatable :: d_t_prevt
  !d_t_prevt=distance between frame at t and frame at previous
  real(kind=single) :: maxjump, minjump, deltajump,xj,yj,cumulative_jump
  ! maxjump, minjump = maximum and minimum lenght of d_t_prevt
  ! deltajump = delta chosen for represent the histogram of the jumps
  ! xj,yj =dummy real(kind=single) variables for helping in the construction of histo_d
  ! cumulative_jump = cumulative probability for the jumps 
  integer(kind=digit) :: nnhisto,nconmin
  !nnhisto = number of elements to be considered in jump histogram representation (parameterized later on)
  !nconmin = minimum number of continous trajectories to consider two clusters
  !          connected
  real(kind=single) :: cutting
  !cutting = distance to generate a new stage in the trajectory (if the jump is bigger than 
  !          cutting, the trajectory will not be considered continuous)
  integer(kind=digit), dimension (:), allocatable :: histo_d
  !histo_d=histogram representation of d_t_prevt
  integer(kind=digit), dimension (:), allocatable :: traj_stage
  ! traj_stage=stage representation of trajectory
  integer(kind=digit), dimension (:,:), allocatable :: connectivity_matrix
  ! traj_stage=stage representation of trajectory
  parameter (nnhisto=50)
  ! number of windows for obtain the jump histogram representation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! ASSIGN DEFAULT VALUES !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  compute_sil=.false.
  sieved=.false.
  nsieve=0
  cutting=1.0
  cutoff=1.0
  nconmin=1
  maxiter=10
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! READ INPUT  file !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  open (11,file=inputfile,status="old",err=110)
  nna=0
  do while (nna.ne.1)
    read (11,'(a800)',end=111)line
    read (line,*) cdum
    write (6,*) line
    if (trim(cdum).eq."NCV") then 
      read (line,*) cdum,nCV
      allocate(cgrid(2,nCV),stat=ierr)
      if (ierr.ne.0) then
        write (6,*) "allocation error cgrid"
        stop
      endif
      allocate(periodic(nCV),stat=ierr)
      if (ierr.ne.0) then
        write (6,*) "allocation error periodic"
        stop
      endif
      allocate(ngrid(nCV),stat=ierr)
      if (ierr.ne.0) then
        write (6,*) "allocation error ngrid"
        stop
      endif
      allocate(grid(nCV),stat=ierr)
      if (ierr.ne.0) then
        write (6,*) "allocation error grid"
        stop
      endif
      allocate(period(nCV),stat=ierr)
      if (ierr.ne.0) then
        write (6,*) "allocation error period"
        stop
      endif
      allocate(CVdum(nCV),stat=ierr)
      if (ierr.ne.0) then
        write (6,*) "allocation error CVdum"
        stop
      endif
    endif
    if (trim(cdum).eq."ARGUMENTS") then
      read (line,*) cdum,algorithm
      select case (trim(algorithm))
        case ("KMEDOIDS") 
          read (line,*) cdum,algorithm,ncl
        case ("GROMACS") 
          read (line,*) cdum,algorithm,cutoff
        case ("GRID") 
          read (line,*) cdum,algorithm
        case default 
          write(6,*) "ALGORITHM ",trim(algorithm),"NOT INCLUDED";stop
      end select  
    endif
    if (trim(cdum).eq."CVGRID") read (line,*) cdum, j,cgrid(1,j),cgrid(2,j),ngrid(j)
    if (trim(cdum).eq."PERIODIC") read (line,*) cdum, (periodic(i),i=1,nCV)
    if (trim(cdum).eq."CONNECTIVITY") read (line,*) cdum,cutting,nconmin
    if (trim(cdum).eq."SILHOUETTE") compute_sil=.true.
    if (trim(cdum).eq."MAXITER") read (line,*) cdum, maxiter
    if (trim(cdum).eq."SIEVING") then 
      sieved=.true.
      read (line,*) cdum,nsieve
    endif
  enddo
111 nna=1
  do j=1,nCV
    period(j)=(cgrid(2,j)-cgrid(1,j))
    grid(j)=period(j)/float(ngrid(j))
  enddo
  open (11,file="COORDINATES",status='old',err=120)
  nna=0
  n=0
  nlines=0
  do while (nna.ne.1)
    n=n+1
    nlines=nlines+1
    read (11,'(a)',err=1111,end=112) longline
    read (longline,*,err=1111) ii,jj,(CVdum(k),k=1,nCV),x,ilab
    if (ilab.eq.0) n=n-1
  enddo
112 nna=1
  nframes=n-1
  nlines=nlines-1
  write (6,*) "NUMBER OF LINES IN COORDINATES FILE=",nlines
  allocate(eq(nframes,2),stat=ierr)
  if (ierr.ne.0) then
    write (6,*) "allocation error eq"
    stop
  endif
  allocate(CV(nframes,nCV),stat=ierr)
  if (ierr.ne.0) then
    write (6,*) "allocation error CV"
    stop
  endif
  allocate(time(nframes),stat=ierr)
  if (ierr.ne.0) then
    write (6,*) "allocation error time"
    stop
  endif
  rewind (11)
  i=0
  do l=1,nlines
    read (11,'(a)') longline
    read (longline,*,err=1111) ii,jj,(CVdum(k),k=1,nCV),x,ilab
    if (ilab.eq.1) then
      i=i+1
      read (longline,*) eq(i,1),eq(i,2),(CV(i,j),j=1,nCV), time(i),ilab
    endif
  enddo
  write (6,*) "nCV=",nCV
  write (6,*) "NUMBER OF ELEMENTS FOR CLUSTERING=",nframes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! STARTING CLUSTERING (ALGORITHM CASE) !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
select case (trim(algorithm))
  
  !!!!!!!!!!!!!!!!!!!!!
  !!! CASE KMEDOIDS !!!
  !!!!!!!!!!!!!!!!!!!!!
  case ("KMEDOIDS")
  ! sieving set up
    nelecl=nframes
    allocate (chosen(nframes),stat=ierr)
    if (ierr.ne.0) then
      write (6,*) "allocation error chosen"
      stop
    endif
    if (nsieve.le.0) sieved=.false.
    if (sieved) then
      if (nsieve.lt.nframes) then
        write (6,*) "Performing random screening for sieved clustering"
        do i=1,nframes
          chosen(i)=.false.
        enddo
        nelecl=nsieve
        allocate(chlist(nelecl),stat=ierr)
        if (ierr.ne.0) then
          write (6,*) "allocation error chlist"
          stop
        endif
        CALL init_random_seed()
        do i=1,nelecl
          nna=0
          do while (nna.ne.1)
            CALL RANDOM_NUMBER(r)
            j=NINT(r*float(nframes-1))+1
            if (.not.chosen(j)) then
              chlist(i)=j
              nna=1
              chosen(j)=.true.
            endif
          enddo
        enddo
    !   deallocate (chosen)
        open (13,file="CHOSEN_LIST.dat")
        write (13,*) chlist(:)
        close (13)
      else
        write (6,*) "The number of elements for sieved clustering is less than the number of frames"
        write (6,*) "PERFORMING STANDARD CLUSTERING"
        sieved=.false.
      endif
    endif
    if (.not.sieved) then
      allocate(chlist(nelecl),stat=ierr)
      if (ierr.ne.0) then
        write (6,*) "allocation error chlist"
        stop
      endif
      do i=1,nelecl
        chlist(i)=i
        chosen(i)=.true.
      enddo
    endif
 !  end sieving set up
    write (6,*) "NUMBER OF SELECTED ELEMENTS=",nelecl
    ndimmat=nelecl*(nelecl-1)/2
    write (6,*) "NUMBER OF DISTANCES TO STORE=",ndimmat
    allocate(cluster(nframes),stat=ierr)
    if (ierr.ne.0) then
      write (6,*) "allocation error cluster"
      stop
    endif
    allocate(tempcluster(nelecl),stat=ierr)
    if (ierr.ne.0) then
      write (6,*) "allocation error tempcluster"
      stop
    endif
    allocate(prob(nelecl),stat=ierr)
    if (ierr.ne.0) then
      write (6,*) "allocation error prob"
      stop
    endif
    allocate(distmat(ndimmat),stat=ierr)
    if (ierr.ne.0) then
      write (6,*) "allocation error distmat"
      stop
    endif
    allocate(ncenter(ncl),stat=ierr)
    if (ierr.ne.0) then
      write (6,*) "allocation error ncenter"
      stop
    endif
    allocate(ncele(ncl),stat=ierr)
    if (ierr.ne.0) then
      write (6,*) "allocation error ncele"
      stop
    endif
    allocate(center(ncl,nCV),stat=ierr)
    if (ierr.ne.0) then
      write (6,*) "allocation error center"
      stop
    endif
    write (6,*) "input reading finished without any error"
    write (6,*) "Starting clustering process"
    write (6,*) "Computing distance matrix:WAIT"
    write (6,*) "storing ",ndimmat," real numbers"
    nna=ndimmat/72
    k=0
    distmat=0.
!$OMP PARALLEL DO SCHEDULE (dynamic,1000) PRIVATE (i,j,k,l,d,limit,dprev,dpost) 
    do k=1,ndimmat
      limit=nelecl-1
      i=1
      do while (k.gt.limit)
        limit=limit+nelecl-i-1
        i=i+1
      enddo
      j=k-(i-1)*nelecl+(i*i+i)/2
      if (i.gt.nelecl) stop "i gt nframes"
      if (j.gt.nelecl) stop "j gt nframes"
      do l=1,nCV
        d=CV(chlist(i),l)-CV(chlist(j),l)
        if (periodic(l)) then
          dprev=d-period(l)
          dpost=d+period(l)
          if (abs(d).gt.abs(dprev)) d=dprev
          if (abs(d).gt.abs(dpost)) d=dpost
        endif
        distmat(k)=distmat(k)+d*d/(grid(l)*grid(l))
      enddo
      distmat(k)=sqrt(distmat(k))
      if (mod(k,nna).eq.0) write (6,'(a1)',advance='no') "."
    enddo
!$OMP END PARALLEL DO
    write (6,*) ""
    write (6,*) "FINISHED, STARTING CLUSTERING PROCESS"
! Randomly assign cluster ncenter (k-means++ algorithm)
    CALL init_random_seed()
    ncenter=0
    CALL RANDOM_NUMBER(r)
    ncenter(1)=NINT(r*float(nelecl-1))+1
    write (6,*) "CLUSTER",1,"CENTER",ncenter(1)
    do l=2,ncl
      nna=1
      do i=1,nelecl
        dmin=9.9e9
        do jj=1,l-1
          j=ncenter(jj)
          if (j.gt.i) then
            k=(i-1)*nelecl-(i*i+i)/2+j
            if (distmat(k).le.dmin) dmin=distmat(k)
          else if (i.gt.j) then
            k=(j-1)*nelecl-(j*j+j)/2+i
            if (distmat(k).le.dmin) dmin=distmat(k)
          else
            dmin=0.0
          endif
        enddo
        prob(i)=dmin*dmin  
      enddo
      do i=2,nelecl
        prob(i)=prob(i)+prob(i-1)
      enddo
      do i=1,nelecl
        prob(i)=prob(i)/prob(nelecl)
      enddo
      CALL RANDOM_NUMBER(r)
      dmin=0.0
      nna=0
      j=1
      if (r.eq.dmin) then
        j=1
        nna=1
      endif
      do while (nna.ne.1)
        if ((r.gt.dmin).and.(r.le.prob(j))) then
          nna=1
        else
          j=j+1
        endif
      enddo
      ncenter(l)=j
    enddo
    write (6,*) ncl,"CENTERS ASSIGNED WITH K-MEANS++ algorithm"
! initialize the cluster array and start the iterative algorithm
    tempcluster(:)=0
    finished=.false.
    niter=0
    do while (.not.finished)
! Compute the new_cluster
      niter=niter+1
      if (niter.gt.maxiter) exit
      finished=.true.
      ndifr=0
!$OMP PARALLEL DO SCHEDULE (dynamic,100) PRIVATE (i,j,k,l,new_cluster,dmin,d,dprev,dpost) 
      do i=1,nelecl
        dmin=9.9e9
        do j=1,ncl
          l=ncenter(j)
          if (l.gt.i) then
            k=(i-1)*nelecl-(i*i+i)/2+l
            if (distmat(k).le.dmin) then
              new_cluster=j
              dmin=distmat(k)
            endif
          else if (i.gt.l) then
            k=(l-1)*nelecl-(l*l+l)/2+i
            if (distmat(k).le.dmin) then
              new_cluster=j
              dmin=distmat(k)
            endif
          else
            dmin=0.0
            new_cluster=j
          endif
        enddo
        if (new_cluster.ne.tempcluster(i)) then
          tempcluster(i)=new_cluster
          ndifr=ndifr+1
          finished=.false.
        endif      
      enddo
!$OMP END PARALLEL DO
      write (6,*) "ITERATION:",niter
      if (.not.finished) then
!$OMP PARALLEL DO SCHEDULE (static) PRIVATE (i,j,k,l,score,minscore)
        do j=1,ncl
          minscore=9.9e9
          do i=1,nelecl
            if (tempcluster(i).eq.j) then
              score=0.0
              do l=1,nelecl
                if (tempcluster(l).eq.j) then
                  if (l.gt.i) then
                    k=(i-1)*nelecl-(i*i+i)/2+l
                    score=score+distmat(k)
                  else if (i.gt.l) then
                    k=(l-1)*nelecl-(l*l+l)/2+i
                    score=score+distmat(k)
                  endif
                endif
              enddo
              if (score.lt.minscore) then
                minscore=score
                ncenter(j)=i
              endif
            endif
          enddo
        enddo
!$OMP END PARALLEL DO
      endif
    enddo
!
    write (6,*) "CLUSTERING FINISHED"
! Undo sieving and compute the cluster assigment of those elements
! that are not chosen
!
    do i=1,nelecl
      cluster(chlist(i))=tempcluster(i)
    enddo
    do i=1,ncl
      ncenter(i)=chlist(ncenter(i))
    enddo
!$OMP PARALLEL DO SCHEDULE(dynamic,100) PRIVATE(i,dmin,j,distance,l,d,dprev,dpost)
    do i=1,nframes
      if (.not.chosen(i)) then
        dmin=9.9e9
        do j=1,ncl
          distance=0.0
          do l=1,nCV
            d=CV(i,l)-CV(ncenter(j),l)
            if (periodic(l)) then 
              dprev=d-period(l)
              dpost=d+period(l)
              if (abs(d).gt.abs(dprev)) d=dprev
              if (abs(d).gt.abs(dpost)) d=dpost
            endif
            distance=distance+d*d/(grid(l)*grid(l))
          enddo
          distance=sqrt(distance)
          if (distance.le.dmin) then
            dmin=distance
            cluster(i)=j
          endif
        enddo
      endif
    enddo
!$OMP END PARALLEL DO
      
       
! Compute the number of element of each cluster (ncele) and
! assign cluster center
    ncele=0
    do j=1,ncl
      l=ncenter(j)
      do k=1,nCV
        center(j,k)=CV(l,k)
      enddo
      do i=1,nframes
        if (cluster(i).eq.j) ncele(j)=ncele(j)+1
      enddo
    enddo
  !!!!!!!!!!!!!!!!!!!!
  !!! CASE GROMACS !!!
  !!!!!!!!!!!!!!!!!!!!
  case ("GROMACS")
    cutoff=cutoff*cutoff
  ! sieving set up
    nelecl=nframes
    if (nsieve.gt.0) sieved=.true.
    if (sieved) then
      if (nsieve.lt.nframes) then
        write (6,*) "Performing random screening for sieved clustering"
        allocate (chosen(nframes),stat=ierr)
        if (ierr.ne.0) then
          write (6,*) "allocation error chosen"
          stop
        endif
        do i=1,nframes
          chosen(i)=.false.
        enddo
        nelecl=nsieve
        allocate(chlist(nelecl),stat=ierr)
        if (ierr.ne.0) then
          write (6,*) "allocation error chlist"
          stop
        endif
        CALL init_random_seed()
        do i=1,nelecl
          nna=0
          do while (nna.ne.1)
            CALL RANDOM_NUMBER(r)
            j=NINT(r*float(nframes-1))+1
            if (.not.chosen(j)) then
              chlist(i)=j
              nna=1
              chosen(j)=.true.
            endif
          enddo
        enddo
        deallocate (chosen)
        open (13,file="CHOSEN_LIST.dat")
        write (13,*) chlist(:)
        close (13)
      else
        write (6,*) "The number of elements for sieved clustering is less than the number of frames"
        write (6,*) "PERFORMING STANDARD CLUSTERING"
        sieved=.false.
      endif
    endif
    if (.not.sieved) then
      allocate(chlist(nelecl),stat=ierr)
      if (ierr.ne.0) then
        write (6,*) "allocation error chlist"
        stop
      endif
      do i=1,nelecl
        chlist(i)=i
      enddo
    endif
 !  end sieving set up
    allocate(cluster(nframes),stat=ierr)
    if (ierr.ne.0) then
      write (6,*) "allocation error cluster"
      stop
    endif
    allocate(assigned(nframes),stat=ierr)
    assigned(:)=.false.
    if (ierr.ne.0) then
      write (6,*) "allocation error assigned"
      stop
    endif
    allocate(nneigh(nframes),stat=ierr)
    if (ierr.ne.0) then
      write (6,*) "allocation error nneigh"
      stop
    endif
    allocate(temp(nframes),stat=ierr)
    if (ierr.ne.0) then
      write (6,*) "allocation error temp"
      stop
    endif
    allocate(tdist(nelecl),stat=ierr)
    if (ierr.ne.0) then
      write (6,*) "allocation error tdist"
      stop
    endif
    write (6,*) "input reading finished without any error"
    write (6,*) "Starting clustering process"
    if (sieved) write (6,*) "SIEVED"
    write (6,*) "Neighbour search: WAIT"
!Generate neighbour list and write it in a file
!do it twice, but it allows us to compute the maximum number of
!neighbours and then allocate it 
    open (12,file="neighbour_list.dat")
    nna=nelecl/72
    maxneigh=0
    do ii=1,nelecl
      nneigh(ii)=0
!$OMP PARALLEL DO PRIVATE (jj,distance,k,d,dprev,dpost) ORDERED
      do jj=1,nelecl
        tdist(jj)=0.0d0
        do k=1,nCV
          d=CV(chlist(ii),k)-CV(chlist(jj),k)
          if (periodic(k)) then
            dprev=d-period(k)
            dpost=d+period(k)
            if (abs(d).gt.abs(dprev)) d=dprev
            if (abs(d).gt.abs(dpost)) d=dpost
          endif
          tdist(jj)=tdist(jj)+d*d/(grid(k)*grid(k))
        enddo
      enddo
!$OMP END PARALLEL DO
      do jj=1,nelecl
        if (tdist(jj).le.cutoff) then
          nneigh(ii)=nneigh(ii)+1
          temp(nneigh(ii))=jj
        endif
      enddo
      if (nneigh(ii).gt.maxneigh) maxneigh=nneigh(ii)
      write(12,*) (temp(j),j=1,nneigh(ii))
      if (mod(ii,nna).eq.0) write (6,'(a1)',advance='no') "."
    enddo
    write (6,*) ""
    write (6,*) "The maximum number of neighbouring frames is ",maxneigh
    deallocate(temp)
    allocate(neighbour(nelecl,maxneigh),stat=ierr)
    if (ierr.ne.0) then
      write (6,*) "allocation error neighbour"
      stop
    endif
!read neighbour list
    rewind (12)
    do i=1,nelecl
      read(12,*) (neighbour(i,j),j=1,nneigh(i))
    enddo
    close(12)
!start iterative clustering
    allocate(ncenter(nframes),stat=ierr)
    if (ierr.ne.0) then
      write (6,*) "allocation error ncenter"
      stop
    endif
    allocate(ncele(nframes),stat=ierr)
    if (ierr.ne.0) then
      write (6,*) "allocation error ncele"
      stop
    endif
    write (6,*) "start iterative clustering"
    nna=1
    ncl=0
    do while (nna.ne.0)
      n=0
! localize frame with the maximum neighbour number
      do ii=1,nelecl
        i=chlist(ii)
        if (.not.assigned(i)) then
          if (nneigh(ii).gt.n) then
            nmax=ii
            n=nneigh(ii)
          endif
        endif
      enddo
! assign the cluster to that frame and to all its neighbours
      ncl=ncl+1
      ncenter(ncl)=chlist(nmax)
      ncele(ncl)=n
      assigned(ncenter(ncl))=.true.
      cluster(ncenter(ncl))=ncl
      do i=1,nneigh(nmax)
        jj=neighbour(nmax,i)
        j=chlist(jj)
        assigned(j)=.true.
        cluster(j)=ncl
      enddo
      write (6,*) "iteration",ncl," number of neighbours:",nneigh(nmax)
! delete the assigned clusters from the neighbour list
!$OMP PARALLEL DO SCHEDULE (dynamic,100) PRIVATE (nns,n,i,l,kk) 
      do kk=1,nelecl
        if (.not. assigned(chlist(kk))) then
          do i=1,nneigh(nmax)
            nns=0
            do n=1,nneigh(kk)
              if(neighbour(kk,n).eq.neighbour(nmax,i)) then
                nns=nns+1
                do l=n+1,nneigh(kk)
                  neighbour(kk,l-1)=neighbour(kk,l)
                enddo
              endif
            enddo
            nneigh(kk)=nneigh(kk)-nns
          enddo
        endif
      enddo
!$OMP END PARALLEL DO
! compute nna and check if clustering loop has finished
      nna=0
      do ii=1,nelecl
        i=chlist(ii)
        if (.not.assigned(i)) nna=nna+1
      enddo
      if (ncl.ge.nelecl) nna=0
    enddo
    write (6,*) "CLUSTERS FROM iterative process:",ncl
!
! Assign the rest of the frames to the cluster for sieved clustering
!
    do i=1,nframes
      if (.not.assigned(i)) then
        if (.not.sieved) then
          write (6,*) "Something goes wrong with clustering"
          write (6,*) "not all the frames has been assigned along the"
          write (6,*) "clustering process, although there is not sieving"
          write (6,*) "STOP"
          stop
        endif
! find dmin
        dmin=9.9E+34
        do l=1,ncl
          j=ncenter(l)
          distance=0.0d0
          do k=1,nCV
            d=CV(i,k)-CV(j,k)
            if (periodic(k)) then
              dprev=d-period(k)
              dpost=d+period(k)
              if (abs(d).gt.abs(dprev)) d=dprev
              if (abs(d).gt.abs(dpost)) d=dpost
            endif
            distance=distance+d*d/(grid(k)*grid(k))
          enddo
          if (distance.le.dmin) then
            dmin=distance
            nmin=l
          endif
        enddo
        if (dmin.lt.cutoff) then
          cluster(i)=nmin
          ncele(nmin)=ncele(nmin)+1
          assigned(i)=.true.
        else
          ncl=ncl+1
          ncenter(ncl)=i
          cluster(i)=ncl
          ncele(ncl)=1
          assigned(i)=.true.
        endif
      endif
    enddo
    allocate(center(ncl,nCV),stat=ierr)
    if (ierr.ne.0) then
      write (6,*) "allocation error center"
      stop
    endif
! assign cluster center
    do j=1,ncl
      l=ncenter(j)
      do k=1,nCV
        center(j,k)=CV(l,k)
      enddo
    enddo
  !!!!!!!!!!!!!!!!!
  !!! CASE GRID !!!
  !!!!!!!!!!!!!!!!!
  case ("GRID")
    write (6,*) "STARTING GRID"
    allocate(cluster(nframes),stat=ierr)
    if (ierr.ne.0) then
      write (6,*) "allocation error cluster"
      stop
    endif
    allocate(pos(nframes),stat=ierr)
    if (ierr.ne.0) then
      write (6,*) "allocation error pos"
      stop
    endif
    allocate(posed(nframes),stat=ierr)
    if (ierr.ne.0) then
      write (6,*) "allocation error posed"
      stop
    endif
    allocate(shift(nCV),stat=ierr)
    if (ierr.ne.0) then
      write (6,*) "allocation error shift"
      stop
    endif
! compute myshift
    myshift=1
    shift(nCV)=1
    do j=nCV,2,-1
      myshift=myshift*ngrid(j)
      shift(j-1)=myshift
    enddo
! Compute position
    do i=1,nframes
      pos(i)=0
      outside=.false.
      do j=1,nCV
        pos(i)=pos(i)+shift(j)*floor((CV(i,j)-cgrid(1,j))/grid(j))
        if ((CV(i,j).lt.cgrid(1,j)).or.(CV(i,j).gt.cgrid(2,j))) outside=.true.
      enddo
      if (outside) pos(i)=-999
    enddo
! compute cluster
    ncl=0
    posed(:)=.false.
    cluster(:)=-999
    do i=1,nframes
      if (pos(i).ne.-999) then
        if (.not.posed(i)) then
          ncl=ncl+1
          cluster(i)=ncl
          do j=i+1,nframes
            if (pos(j).eq.pos(i)) then
              posed(j)=.true.
              cluster(j)=ncl
            endif
          enddo
        endif
      endif
    enddo
!
! Compute cluster center and number of elements in the cluster
!
    allocate(center(ncl,nCV),stat=ierr)
    if (ierr.ne.0) then
      write (6,*) "allocation error center"
      stop
    endif
    allocate(bincl(ncl),stat=ierr)
    if (ierr.ne.0) then
      write (6,*) "allocation error bincl"
      stop
    endif
    allocate(ncele(ncl),stat=ierr)
    if (ierr.ne.0) then
      write (6,*) "allocation error ncele"
      stop
    endif
    ncele(:)=0
    do i=1,nframes
      if (pos(i).ne.-999) then
        bincl(cluster(i))=pos(i)
        ncele(cluster(i))=ncele(cluster(i))+1
      endif
    enddo
    do i=1,ncl
      n=bincl(i)
      do j=1,nCV
        l=floor(float(n)/float(shift(j)))
        center(i,j)=cgrid(1,j)+(l+0.5)*grid(j)
        n=n-l*shift(j)
      enddo
    enddo
    write (6,*) "GRIDDING PERFORMED"  
  !!!!!!!!!!!!!!!!!!!!
  !!! CASE DEFAULT !!!
  !!!!!!!!!!!!!!!!!!!!
  case default 
    write(6,*) "ALGORITHM ",trim(algorithm),"NOT INCLUDED";stop
end select
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! COMPUTE CONNECTIVITY !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
select case (trim(algorithm))
  case ("GRID")
    allocate(connectivity_matrix(ncl,ncl),stat=ierr)
    if (ierr.ne.0) then
      write (6,*) "allocation error connectivity_matrix"
      stop
    endif
    connectivity_matrix=0
    do i=1,ncl-1
      do j=i+1,ncl
        cent_distance=0.0
        do k=1,nCV
          d=center(i,k)-center(j,k)
          if (periodic(k)) then
            dprev=d-period(k)
            dpost=d+period(k)
            if (abs(d).gt.abs(dprev)) d=dprev
            if (abs(d).gt.abs(dpost)) d=dpost
          endif
          cent_distance=cent_distance+d*d/(grid(k)*grid(k))
        enddo
        if (cent_distance.lt.1.000002) then
          connectivity_matrix(i,j)=connectivity_matrix(i,j)+1
          connectivity_matrix(j,i)=connectivity_matrix(j,i)+1
        endif
      enddo
    enddo
    open (16,file="connectivity") 
    do i=1,ncl-1
      do j=i+1,ncl
        if ((connectivity_matrix(i,j).ge.1).or.(connectivity_matrix(j,i).ge.1)) then
          write (16,*) i,j
        endif
      enddo
    enddo
    close(16)
  case default
    allocate(d_t_prevt(nframes),stat=ierr)
    if (ierr.ne.0) then
      write (6,*) "allocation error d_t_prevt"
      stop
    endif
    allocate(traj_stage(nframes),stat=ierr)
    if (ierr.ne.0) then
      write (6,*) "allocation error traj_stage"
      stop
    endif
    allocate(histo_d(nnhisto),stat=ierr)
    if (ierr.ne.0) then
      write (6,*) "allocation error histo_d"
      stop
    endif
! compute the distance between t and t-1 for all frames, 
! and plot it in a histogram (the histogram will help to 
! stablish the continuity conditions to look for the connectivity)
    minjump=9.9d9
    maxjump=-9.9d9
    open (16,file="distance_t_tprev.dat")
    traj_stage(1)=0
    do i=2,nframes
      j=i-1
      d_t_prevt(i)=0.0d0
      do k=1,nCV
        d=CV(i,k)-CV(j,k)
        if (periodic(k)) then
          dprev=d-period(k)
          dpost=d+period(k)
          if (abs(d).gt.abs(dprev)) d=dprev
          if (abs(d).gt.abs(dpost)) d=dpost
        endif
        d_t_prevt(i)=d_t_prevt(i)+d*d/(grid(k)*grid(k))
      enddo
      d_t_prevt(i)=sqrt(d_t_prevt(i))
      if (d_t_prevt(i).gt.maxjump) maxjump=d_t_prevt(i)
      if (d_t_prevt(i).lt.minjump) minjump=d_t_prevt(i)
      if (d_t_prevt(i).gt.cutting) then
        traj_stage(i)=traj_stage(i-1)+1
      else
        traj_stage(i)=traj_stage(i-1)
      endif
      write (16,*) i,d_t_prevt(i),traj_stage(i)
    enddo
    close (16)
    deltajump=(maxjump-minjump)/float(nnhisto)
    cumulative_jump=0.0
    open (16,file="jump_histogram.dat")
    do i=1,nnhisto
      histo_d(i)=0
      xj=minjump+float(i-1)*deltajump
      yj=minjump+float(i)*deltajump
      do j=2,nframes
        if ((d_t_prevt(j).gt.xj).and.(d_t_prevt(j).lt.yj)) histo_d(i)=histo_d(i)+1
      enddo
      cumulative_jump=cumulative_jump+histo_d(i)/float(nframes-1)*100.
      write (16,*) (yj+xj)/2.,histo_d(i)/float(nframes-1)*100.,cumulative_jump
    enddo
    close (16)
    allocate(connectivity_matrix(ncl,ncl),stat=ierr)
    if (ierr.ne.0) then
      write (6,*) "allocation error connectivity_matrix"
      stop
    endif
    connectivity_matrix=0
    do i=1,nframes-1
      j=i+1
      if (traj_stage(i).eq.traj_stage(j)) connectivity_matrix(cluster(i),cluster(j))=connectivity_matrix(cluster(i),cluster(j))+1
    enddo
    open (16,file="connectivity_matrix")
    do i=1,ncl
      write (16,'(100i6)') (connectivity_matrix(i,j),j=1,ncl)
    enddo
    close (16)
    open (16,file="connectivity") 
    do i=1,ncl-1
      do j=i+1,ncl
        if ((connectivity_matrix(i,j).ge.nconmin).or.(connectivity_matrix(j,i).ge.nconmin)) then
          write (16,*) i,j
        endif
      enddo
    enddo
    close (16)
end select
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Compute silhouette if asked for !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (compute_sil) then
    allocate(silhouette(nframes),stat=ierr)
      if (ierr.ne.0) then
      write (6,*) "allocation error silhouette"
      stop
    endif
    allocate(average_dissimilarity(ncl),stat=ierr)
    if (ierr.ne.0) then
      write (6,*) "allocation error average_dissimilarity"
      stop
    endif
    allocate(tdist(nframes),stat=ierr)
    if (ierr.ne.0) then
      write (6,*) "allocation error tdist"
      stop
    endif
    open (16,file="silhouette.dat")
    sil_average=0.0
    do i=1,nframes
      do j=1,ncl
        average_dissimilarity(j)=0.0
      enddo
!$OMP PARALLEL DO PRIVATE (j,k,d,dprev,dpost) ORDERED
      do j=1,nframes
        tdist(j)=0.
        do k=1,nCV
          d=CV(i,k)-CV(j,k)
          if (periodic(k)) then
            dprev=d-period(k)
            dpost=d+period(k)
            if (abs(d).gt.abs(dprev)) d=dprev
            if (abs(d).gt.abs(dpost)) d=dpost
          endif
          tdist(j)=tdist(j)+d*d/(grid(k)*grid(k))
        enddo
        tdist(j)=sqrt(tdist(j))
      enddo
!$OMP END PARALLEL DO 
      do j=1,nframes
        average_dissimilarity(cluster(j))=average_dissimilarity(cluster(j))+tdist(j)
      enddo
!Compute the neighbouring cluster (l)
      average_dissimilarity=average_dissimilarity/float(ncele)
      distance=9.9e9
      do j=1,ncl
        if (j.ne.cluster(i)) then
          if (average_dissimilarity(j).lt.distance) then
            distance=average_dissimilarity(j)
            l=j
          endif
        endif
      enddo 
!Compute silhouette
      if (average_dissimilarity(cluster(i)).lt.average_dissimilarity(l)) then
        silhouette(i)=1.-(average_dissimilarity(cluster(i))/average_dissimilarity(l))
      else 
        silhouette(i)=(average_dissimilarity(l)/average_dissimilarity(cluster(i)))-1.
      endif
      write (16,*) i,silhouette(i),cluster(i),average_dissimilarity(cluster(i)),l,average_dissimilarity(l)
      sil_average=sil_average+silhouette(i)/float(nframes)
    enddo
    write (6,*) "AVERAGE SILHOUETTE=",sil_average
    write (16,*) "AVERAGE SILHOUETTE=",sil_average
    close (16)
    open (16,file="cluster_silhouette.dat")
    do i=1,ncl
      sil_average=0.0
      do j=1,nframes
        if (cluster(j).eq.i) then
          sil_average=sil_average+silhouette(j)
        endif
      enddo
      sil_average=sil_average/float(ncele(i))
      write (16,*) "CLUSTER",i,"SILHOUETTE",sil_average
    enddo
    close(16)
  endif
!!!!!!!!!!!!!!!!!!!!
!!! WRITE OUTPUT !!!
!!!!!!!!!!!!!!!!!!!!

  open (12,file="MICROSTATES")
  do i=1,ncl
    write(12,*) i,ncele(i),(center(i,j),j=1,nCV),1000.,-999
  enddo
  close (12)
  open (13,file="FRAMES")
  rewind (11)
  i=0
  do l=1,nlines
    read (11,'(a)') longline
    read (longline,*) ii,jj,(CVdum(k),k=1,nCV),x,ilab
    if (ilab.eq.1) then
      i=i+1
      write (13,*) l-1, cluster(i),eq(i,2),eq(i,1),time(i)
    else
      write (13,*) l-1, -999 ,jj,ii,x
    endif
  enddo
  close (13)
  stop
!!!!!!!!!!!!!!!!!!!!!!!!
!!! ERROR MANAGEMENT !!!
!!!!!!!!!!!!!!!!!!!!!!!!

110 write (6,*) "ERROR OPENING INPUT FILE",trim(inputfile)
    write (6,*) "--------- STOP ----------"
  stop
1111 write (6,*) "ERROR READING COORDINATES in line ",nlines
     write (6,*) "--------- STOP ----------"
  stop
120 write (6,*) "ERROR OPENING COORDINATES file"
    write (6,*) "--------- STOP ----------"
  stop
  end subroutine make_microstates
!!!!!!!!!!!!!!!!!!!
!!! SUBROUTINES !!!
!!!!!!!!!!!!!!!!!!!

SUBROUTINE init_random_seed()
 INTEGER i,n, clock
 INTEGER seed (100)
 CHARACTER(len=256) :: fixed_seed
 INTEGER length, status
 
 n=100
 CALL RANDOM_SEED(size = n)

 CALL get_environment_variable("MK_MICROSTATES_FIXED_SEED", fixed_seed, length, status )
 if(status==1) then
    ! Env var does not exist
    CALL SYSTEM_CLOCK(COUNT=clock)
 else
    clock=0
    write (6,*) "--------- USING FIXED SEED ----------"
 end if


 do i=1,100
   seed (i) = clock + 37 * (i-1)
 enddo
 CALL RANDOM_SEED(PUT = seed)
END SUBROUTINE

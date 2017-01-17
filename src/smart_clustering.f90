 subroutine smart_cluster(inputfile)
 ! CRD(:,:) Coordinates for all the elements
 ! D(:) Distance matrix, stored as vector
 ! rho(:) density vector == exp(-F/KT)
 ! delta(:) cluster center likelihood parameter
 ! cluster (:) cluster assignation parameter
 ! halo (:) cluster assignation taking into account halo 
 !
 ! "Clustering by fast search-and-find of density peaks"  Alejandro Rodriguez & Alessandro Laio, 2013

! VARIABLE DECLARATION
 implicit none
 character*100 inputfile
 integer,parameter:: digit= selected_int_kind(11) !! digits in the integer numbers considered
 integer,parameter:: nhisto=1000 !! Histogram for deciding the cutoff and for make a quantizied square
 integer :: int_arg
 integer :: nna
 integer :: ncount
 integer(kind=digit) :: nmicro,pmicro,bmicro
 real,allocatable,dimension(:) :: microcoord
 real :: Fmicro
 integer(kind=digit) :: i,j,k,l,m,n,nvar,nele,ndistances,limit,ii,jj,narg,nn
 integer(kind=digit) :: ngood,ncol
 integer(kind=digit),parameter :: zero=0
 integer(kind=digit) :: nlines,ncluster
 real :: dist,x,y,z,r,selection,rr
 real :: dprev,dpost,dact
 real :: cutoff,decay,percent
 real :: rcut,dcut
 integer(kind=digit) :: nmin
 integer(kind=digit) :: aneigh
 integer(kind=digit) :: jmax
 character(len=5000) :: line, arg1,arg2
 character(len=50) :: fCV,distype,dentype,perfile,neightype
 real, allocatable, dimension (:,:) ::CRD
 real, allocatable, dimension (:) :: rho,delta,D
 real, allocatable, dimension (:) :: spacer
 real, allocatable, dimension (:) :: period
 real, allocatable, dimension (:) :: gama
 real :: KT
 real :: rho_border
 real :: disneigh
 logical, allocatable, dimension (:) :: limited
 real :: dmin,dmax,dhisto,dlow,dhigh,ddelta              ! related with autocrowding density
 integer(kind=digit), allocatable, dimension (:) :: equiv ! 
 integer(kind=digit), allocatable, dimension (:) :: sort ! sort(i) 
 integer(kind=digit), allocatable, dimension (:) :: invsort !invsort(i)
 integer(kind=digit), allocatable, dimension (:) :: cluster
 integer(kind=digit), allocatable, dimension (:) :: halo
! cluster characteristics
 integer(kind=digit), allocatable, dimension (:) :: center !element center of the cluster
 integer(kind=digit), allocatable, dimension (:) :: scenter !sorted element center of the cluster
 integer(kind=digit), allocatable, dimension (:) :: cluster_elements ! cluster population
 integer(kind=digit), allocatable, dimension (:) :: halo_elements !  cluster population ignoring halo elements
 real, allocatable, dimension (:) :: center_density !  density at the center of the cluster
 real, allocatable, dimension (:) :: border_density !  density at the border of the cluster
 integer(kind=digit), allocatable, dimension (:) :: border_cluster !  cluster that determines border density
 integer(kind=digit), allocatable, dimension (:) :: border_element !  element that determines border density
 integer(kind=digit), allocatable, dimension (:) :: border_other !  other element that determines border density
 logical ::fmet
 logical ::blanko
 logical ::wrmat ! write the unformatted distance matrix
 logical ::wrdis ! output the reordered distance matrix
 logical ::noclus! allows to just compute the decision graph and gamma log log plot
 integer(kind=digit), allocatable, dimension (:) :: nnsort !sorting for square
! END VARIABLE DECLARATION
! PREVIOUS WORK TO START ALGORITHM
! write usage
!write (6,*) "Usage: "
! set default values for arguments
 disneigh=1.0
 write (6,*) "let's start"

! READ INPUT FILES (coordinates or distances)
 open (99,file=inputfile,status='old',err=9999)
 read (99,*) fCV ! read MICROSTATES file name
 read (99,*) KT  ! KT value to get density from Free energy
 read (99,*) nvar  ! KT value to get density from Free energy
 open (11,file=fCV,status='old',err=1100)
! get number of columns (variables,ncol) from the first line
 ncol=nvar+4
 allocate(microcoord(nvar))
! get number of files (elements,nele) by reading the file and rewind
 nele=0
 nlines=0
 open (14,file='equiv.dat')
 do 
   read (11,'(a5000)',end=110) line
   nlines=nlines+1
   read (line,*) nmicro,pmicro,(microcoord(i),i=1,nvar),Fmicro,bmicro
   if (Fmicro.lt.999.0) then
     nele=nele+1
     write (14,*) nele,nlines
   endif
 enddo
 close (14)
110 continue
 rewind (11)
 write (6,*) nlines," LINES"
 write (6,*) nlines," MICROSTATES"
 write (6,*) nele," ELEMENTS"
 write (6,*) nvar," VARIABLES"
 ndistances=nele*(nele-1)/2
! Allocate memory
 allocate (spacer(nvar))
 allocate (period(nvar))
 allocate (D(ndistances))
 allocate (CRD(nele,nvar)) 
 allocate (rho(nele)) 
! 
! init periodic variables
 period(:)=0.0
 spacer(:)=1.0
! Read periodicity 
 do 
   read (99,*,end=120) i,x,y
   period(i)=x
   spacer(i)=y
 enddo
 120 continue
 close (99)
! Read coodinates
 nn=0
 k=0
 do i=1,nlines
   read (11,'(a5000)',end=110) line
   read (line,*) nmicro,pmicro,(microcoord(j),j=1,nvar),Fmicro,bmicro
   if (Fmicro.lt.999.0) then
     k=k+1
     rho(k)=exp(-Fmicro/KT)
     do j=1,nvar
       CRD(k,j)=microcoord(j)
     enddo
   endif
 enddo
 close (11)
 write (6,*) "Coordinates readed"
 D(:)=0.
 write (6,*) "Computing distances"
 ncount=ndistances/72
!$OMP PARALLEL DO SCHEDULE (dynamic,1000) PRIVATE (i,j,k,l,limit,dact,dprev,dpost) 
 do k=1,ndistances
   limit=nele-1
   i=1
   do while (k.gt.limit)
     limit=limit+nele-i-1
     i=i+1
   enddo
   j=k-(i-1)*nele+(i*i+i)/2
   do l=1,nvar
     dact=CRD(i,l)-CRD(j,l)
     dprev=dact-period(l)
     dpost=dact+period(l)
     if (abs(dact).gt.abs(dprev)) dact=dprev
     if (abs(dact).gt.abs(dpost)) dact=dpost
     D(k)=D(k)+(dact*dact)/spacer(l)
   enddo
   D(k)=sqrt(D(k))
   if (mod(k,ncount).eq.0) write (6,'(a1)',advance="NO") "."
 enddo
!$OMP END PARALLEL DO
 write (6,*) ""
 deallocate (period)
 deallocate (spacer)
 deallocate (CRD)
 deallocate (microcoord)
 allocate (cluster(nele)) 
 allocate (delta(nele)) 
 allocate (sort(nele)) 
 allocate (invsort(nele)) 
! Compute delta, first sort
! sort(i)=real index of the sorted element i
! sort(1)= index of the element with highest density
 do i=1,nele
   sort(i)=i
 enddo
 do i=1,nele
   do j=i,nele
     if (rho(j).ge.rho(i)) then
       r=rho(i)
       rho(i)=rho(j)
       rho(j)=r
       ii=sort(i)
       sort(i)=sort(j)
       sort(j)=ii
     endif 
   enddo
 enddo
! sort(i)=real index of the sorted element i
! invsort(ii)=order in density of the element with the index ii
 do i=1,nele
   ii=sort(i)
   invsort(ii)=i
 enddo
! Compute delta 
 r=-9.9e9
 rr=9.9e30
 do ii=2,nele
   delta(ii)=9.9e9
   i=sort(ii)
   do jj=ii-1,1,-1
     j=sort(jj)
     if (j.gt.i) then
       k=(i-1)*nele-(i*i+i)/2+j
     else 
       k=(j-1)*nele-(j*j+j)/2+i
     endif
     if (D(k).le.delta(ii)) delta(ii)=D(k)
   enddo
   if (delta(ii).ge.r) r=delta(ii)
   if (delta(ii).lt.rr) rr=delta(ii)
 enddo 
! for visualization pourpouses, scale delta (1), previously assigned delta(1)=dmax
 delta(1)=r+(r-rr)*0.5
!write & plot decision graph
 open (14,file="dec.dat")
 do i=1,nele
   write (14,*) i,rho(i),delta(i),sort(i)
 enddo
 close(14)
!open (14,file="dec.gpl")
!write (14,*) "plot 'dec.dat' u 2:3 pt 7 ps 1 t ''"
!close (14)
!call system ("gnuplot -persist dec.gpl")
! read rho and delta cutoff
 fmet=.false.
 do while (.not.fmet)
   call sleep (1)
   INQUIRE(FILE="input.smart", EXIST=fmet)
 enddo
 open (15,file="input.smart")
 read(15,*) rcut
 read(15,*) dcut
! perform cluster assignation... 
 cluster(:)=0
 ncluster=1
 cluster(sort(1))=1
 do i=2,nele
   ii=sort(i)
   if ((rho(i).gt.rcut).and.(delta(i).gt.dcut)) then
     ncluster=ncluster+1
     cluster(ii)=ncluster
   else
     dmin=9.9e9
     do j=1,i-1
       jj=sort(j)
       l=max(ii,jj)
       m=min(ii,jj)
       k=(m-1)*nele-(m*m+m)/2+l
       if (D(k).lt.dmin) then
         dmin=D(k)
         cluster(ii)=cluster(jj)
       endif
     enddo
   endif
 enddo  
 write (6,*) "TOTAL NUMBER OF CLUSTERS:",ncluster
 allocate (halo(nele)) 
 allocate (limited(ncluster))
 allocate (border_other(ncluster)) 
 allocate (cluster_elements(ncluster))
 allocate (halo_elements(ncluster))
 allocate (center_density(ncluster))
 allocate (border_density(ncluster))
 allocate (border_cluster(ncluster))
 allocate (border_element(ncluster))
 halo(:)=cluster(:)
 limited(:)=.false.
 halo_elements(:)=0
 cluster_elements(:)=0
 border_element(:)=0
 border_other(:)=0
! found border element and density for each cluster
 do ii=1,ncluster
   rho_border=-9.9e9
   do i=1,nele
     if (cluster(i).eq.ii) then
       do j=1,nele
         if (cluster(j).ne.ii) then
           if (j.gt.i) then
             k=(i-1)*nele-(i*i+i)/2+j
           else 
             k=(j-1)*nele-(j*j+j)/2+i
           endif
           if (D(k).le.disneigh) then
             if (rho_border.le.rho(invsort(j))) then
               rho_border=rho(invsort(j))
               border_element(ii)=i
               border_other(ii)=j
               border_cluster(ii)=cluster(j)
               limited(ii)=.true.
             endif
           endif
         endif
       enddo
     endif
   enddo
 enddo
 if (.not.wrdis) deallocate (D)
 do ii=1,ncluster
   if (limited(ii)) then
     border_density(ii)=0.5*(rho(invsort(border_element(ii)))+rho(invsort(border_other(ii))))
   else
     border_density(ii)=-9.9e9
   endif
 enddo
 do i=1,nele
   cluster_elements(cluster(i))=cluster_elements(cluster(i))+1
   if (rho(invsort(i)).le.border_density(cluster(i))) then
     halo(i)=0
   else
     halo_elements(cluster(i))=halo_elements(cluster(i))+1
   endif
 enddo
!get center properties
 allocate (center(ncluster)) 
 allocate (scenter(ncluster)) 
 nn=1
 i=0
 do while (nn.le.ncluster)
   i=i+1
   ii=sort(i) 
   if (cluster(ii).eq.nn) then
     center(nn)=ii
     scenter(nn)=i
     write (6,*)"CLUSTER:",nn,"border rho:",border_density(nn),"center rho:",rho(i)
     write (6,*)"cluster elements:",cluster_elements(nn),"halo elements:",halo_elements(nn)
     nn=nn+1
   endif
 enddo
 ngood=0
 do i=1,nele
   if (halo(i).ne.0) ngood=ngood+1
 enddo
 close(14)
 write (6,*) "CORE SIZE=",ngood
 open (14,file="equiv.dat")
 open (15,file="smart.out")
 write (15,*) ncluster
 jmax=nlines
 allocate (equiv(jmax))
 equiv(:)=0
 do
   read (14,*,end=1142) i,j
   equiv(j)=i
 enddo
1142 continue
 do j=1,jmax
   i=equiv(j)
   if (i.gt.0) then
     l=halo(i)
     if (i.ne.center(l)) then
       write (15,*) j,l
     else
       write (15,*) j,-l
     endif
   else
     write (15,*) j,zero
   endif
 enddo
 close(14)
 close(15)
! Final deallocation... it could be still optimized...
 deallocate (cluster) 
 deallocate (sort) 
 deallocate (invsort) 
 deallocate (delta) 
 deallocate (halo) 
 deallocate (limited)
 deallocate (border_other) 
 deallocate (cluster_elements)
 deallocate (halo_elements)
 deallocate (center_density)
 deallocate (border_density)
 deallocate (border_cluster)
 deallocate (border_element)
 deallocate (center) 
 deallocate (scenter) 
 stop
 9999 write (6,*) "FILE: ",trim(inputfile)," does not exists"
      write (6,*) "STOP"
      stop
 1100 write (6,*) "FILE: ",trim(fCV)," does not exists"
      write (6,*) "STOP"
      stop
 1200 write (6,*) "FILE: ",trim(perfile)," does not exists"
      write (6,*) "STOP"
      stop
 2100 write (6,*) "FILE rho.dat does not exists"
      write (6,*) "STOP"
      stop
 1300 write (6,*) "FILE: distmat.dat does not exists"
      write (6,*) "and is needed for using the -distance read option"
      write (6,*) "STOP"
      stop
 1400 format (i3,x,i6,x,i6,x,i6,x,f9.2,x,f9.2,x,i6,x,i6,x,i6,x,f9.2,x,f9.2)
  end subroutine smart_cluster

program arrange 

    !New version of arrange by Leo Patrick Mulholland 27.09.24

    !Merges multiple OMEGA files (only formatted right now)
    !In the structure OMEGAXXXX - where XXXX is a number e.g 0001,0123
    !Writes a formatted OMEGAZ file.

    !The energies are quick sorted (qsort), and the corresponding OMEGA values written 
    !according to a pointer array produced by the subroutine qsort.

    !I did not write the quick sort subroutine. I lifted it from the previous version of
    !the code.
    !Credit:        
        !    NON-RECURSIVE STACK VERSION OF QUICKSORT FROM N.WIRTH'S PASCAL
        !    BOOK, 'ALGORITHMS + DATA STRUCTURES = PROGRAMS'. 


    implicit none

    integer :: ii 
    integer,parameter :: max_iter = 9999 
    integer :: num_points_om, total_num_points,offset
    character*9 :: current_file
    logical :: file_exists
    integer :: num_files 
    integer :: index_exists(max_iter)
    integer :: NZED,NELEC,NAST,NUM_TRAN

    real*8:: t1,t2

    !for each individual omega
    real*8,allocatable :: energies_bound(:)
    real*8,allocatable :: energies_incident_om(:)
    real*8,allocatable :: omega_om(:,:)

    !for total omega
    real*8,allocatable :: energies_incident_all(:)
    real*8,allocatable :: omega_all(:,:)

    integer,allocatable :: angular_momentum_dump(:)
    integer,allocatable :: bigL(:),multiplicity(:),indexpointer(:)

    num_files = 0
    total_num_points = 0
    !scan through all omegas and find the dimensions

    do ii = 0,max_iter

        current_file = file_name(ii)
        inquire(file=current_file,exist=file_exists)
        if (file_exists) then 
            num_files = num_files + 1
            index_exists(num_files) = ii
            call initial_read(current_file,num_points_om,NZED,NELEC,NAST,NUM_TRAN)
            print*,current_file, 'exists',num_points_om,' pts'
            total_num_points = total_num_points + num_points_om
        end if  

    end do

    !house keeping, make this better formatted

    print*,'total number of points      : ',total_num_points
    print*,'NUCLEAR CHARGE              : ',NZED
    print*,'NUMBER ELECTRONS            : ',NELEC
    print*,'number of bound states      : ',nast
    print*,'number of transitions       : ',NUM_TRAN
    print*,'number of files to be sorted: ',num_files

    !now do proper read.

    allocate(omega_all(total_num_points,NUM_TRAN))
    allocate(angular_momentum_dump(2*nast))
    allocate(energies_incident_all(total_num_points))
    allocate(indexpointer(total_num_points))
    offset = 0

    do ii = 1,num_files
        current_file = file_name(index_exists(ii))
        call read_omega(NZED,NELEC,NAST,omega_om,energies_bound,energies_incident_om,bigL,multiplicity,num_points_om,current_file)
        energies_incident_all(1+offset:num_points_om+offset) = energies_incident_om
        omega_all(1+offset:num_points_om+offset,:) = omega_om
        offset = offset + num_points_om
        deallocate(omega_om,energies_incident_om)
        if (ii.ne.num_files) then
            !lazy soln that works
            deallocate(energies_bound,bigL,multiplicity)
        end if 
    end do 
    do ii =1,total_num_points
        indexpointer(ii) = ii
    end do 

    !Quick Sorting Energies.
    write(*,10)
    call cpu_time(t1)
    call qsort(energies_incident_all,total_num_points,indexpointer)
    call cpu_time(t2)
    write(*,11) t2-t1


    !Write out the new OMEGA file, using the sorted indexpointer.
    call write_arranged_omega(NZED,NELEC,NAST,omega_ALL,&
    energies_bound,energies_incident_all,&
    bigL,multiplicity,total_num_points,num_tran,indexpointer)

    10 format ('Quick Sorting Energies.')
    11 format ('Sort time: ',F5.1,' seconds.')

    contains
    
    function file_name(number)
        !finds omega file name from an integer.
        character*9 :: file_name
        integer :: number
        integer :: i1,i2,i3,i4
        character*1 NUM(0:9),filec,filed,fileu,filev
        DATA NUM /'0','1','2','3','4','5','6','7','8','9'/

        i1= number/1000
        i2=(number-1000*(number/1000))/100
        i3=(number-(1000*(number/1000))-i2*100)/10
        i4=(number-(1000*(number/1000))-i2*100)-i3*10

        filec = NUM(i1)
        filed = NUM(i2)
        fileu = NUM(i3)
        filev = NUM(i4)

        file_name = 'OMEGA'//filec//filed//fileu//filev
        
    end function

    subroutine initial_read(file_name,num_points,NZED,NELEC,NAST,NUM_TRAN)
        character*9 :: file_name
        integer :: num_points
        integer :: NZED,NELEC,NAST
        integer :: NUM_TRAN
        open(1,file = file_name,status = 'old')
        READ(1,*)NZED,NELEC
        READ(1,*)NAST,num_points,NUM_TRAN
        close(1)
    end subroutine
    
    subroutine read_omega(NZED, &
        NELEC,&
        NAST ,&
        OMEGA,&
        energies_bound,&
        energies_incident,&
        bigL,&
        mult,&
        num_points,&
        omega_path)

        !reads OMEGA, funny that.
        !returns, charge, num electrons,
        !omega, 
        !Bound and incident energies in Rydbergs/(z^2). 

        !Parameters determined from OMEGA
        integer,intent(inout) :: NZED , NELEC
        integer               :: NAST , NUM_POINTS,NUM_TRAN

        !Allocatables 
        integer,allocatable :: bigL(:),mult(:)
        !integer*8,allocatable :: stat_weight(:)
        real*8   ,allocatable :: energies_bound(:)
        real*8   ,allocatable :: omega(:,:)
        real*8   ,allocatable :: energies_incident(:)

        !Input options
        !logical      :: debug
        character*9 :: omega_path

        !Iterators
        integer*8 :: i
        integer*8 :: j

        !Read in times
        real*8 :: read_time_1,read_time_2
        !debug = .false.
        !print*,'-------------------------------------'
        !print*,'Entering routine read_omega()'

        call cpu_time(read_time_1)

        open(1,file=omega_path,status='old',form='formatted')

        !print*,'reading file,' ,omega_path

        !reading in basic data.
        !this allows allocataion of energies, OMEGA

        read(1,*) NZED,NELEC 
        read(1,*) NAST,NUM_POINTS,NUM_TRAN 

        !allocating arrays.
        allocate(energies_bound(NAST))  
        allocate(bigL(NAST))  
        allocate(mult(NAST))  
        allocate(energies_incident(NUM_POINTS))  
        allocate(omega(NUM_POINTS,NUM_TRAN))  

        !reading in multiplicity, big L, bound energies.
        read(1,*) (mult(i),bigL(i),i=1,nast)
        read(1,*) (energies_bound(i),i=1,nast)

        !regular omega format
        !should probably add an option for unformatted
        !but i think it's as simple as changing the open statement
        do i = 1, num_points
            read(1,*) energies_incident(i),(omega(i,j),j=1,NUM_TRAN)
        end do     

        close(1)

        call cpu_time(read_time_2)

        write(*,100) omega_path,read_time_2-read_time_1

        100 format (A9,' read time: ',F6.2,' seconds.')

        !110 format ('Number of transitions unexpected. \n & 
        !      &This code assumes no elastic transitions. & 
        !      &Found transitions ',i6, 'expected (NAST*(NAST-1)/2) = ',i6)

    end subroutine

    subroutine write_arranged_omega(NZED, &
        NELEC,&
        NAST ,&
        OMEGA,&
        energies_bound,&
        energies_incident,&
        bigL,&
        mult,&
        num_points,&
        num_tran,&
        index_pointer&
        )

        integer :: NZED,NELEC,NAST,num_points,num_tran
        real*8 :: energies_bound(nast)
        integer :: bigL(nast),mult(nast)
        real*8 :: energies_incident(num_points),t1,t2
        integer :: index_pointer(num_points)
        real*8 :: omega(num_points,num_tran) 
        integer :: ii,index
        call cpu_time(t1)
        open(1,file='OMEGAZ')

        WRITE(1,*) NZED,NELEC
        WRITE(1,*)NAST,num_points,num_tran
        WRITE(1,*)(mult(ii),bigL(ii),ii=1,nast)
        WRITE(1,270)(energies_bound(ii),ii=1,nast)
        do ii =1,num_points
            
            index = index_pointer(ii)
            write(1,380) energies_incident(ii),omega(index,:)
        end do 

        close(1)
        call cpu_time(t2)
        write(*,15)t2-t1
        15 format ('Write to OMEGAZ time: ',F5.1,' seconds.')

        270   FORMAT(1p5e16.6)
        380   FORMAT(1pe14.8,6(1pe11.3)/(14x,6(e11.3)))
    end subroutine

    SUBROUTINE qsort(a, n, t)

        !I copy and pasted this from the original code. L

        !    NON-RECURSIVE STACK VERSION OF QUICKSORT FROM N.WIRTH'S PASCAL
        !    BOOK, 'ALGORITHMS + DATA STRUCTURES = PROGRAMS'.
        !    REAL*8 PRECISION, ALSO CHANGES THE ORDER OF THE ASSOCIATED ARRAY T.
                IMPLICIT NONE
                INTEGER, INTENT(IN)    :: n
                REAL*8, INTENT(INOUT)    :: a(n)
                INTEGER, INTENT(INOUT) :: t(n)
        !    Local Variables
                INTEGER    :: i, j, k, l, r, s, stackl(15), stackr(15), ww
                REAL*8    :: w, x
                s = 1
                stackl(1) = 1
                stackr(1) = n
        !    KEEP TAKING THE TOP REQUEST FROM THE STACK UNTIL S = 0.
         10     CONTINUE
                l = stackl(s)
                r = stackr(s)
                s = s - 1
        !    KEEP SPLITTING A(L), ... , A(R) UNTIL L >= R.
         20     CONTINUE
                i = l
                j = r
                k = (l+r) / 2
                x = a(k)
        !    REPEAT UNTIL I > J.
              DO
              DO
               IF (a(i).LT.x) THEN    ! Search from lower end
              i = i + 1
              CYCLE
               ELSE
               EXIT
               END IF
              END DO
              DO
              IF (x.LT.a(j)) THEN    ! Search from upper end
              j = j - 1
              CYCLE
              ELSE
              EXIT
              END IF
              END DO
              IF (i.LE.j) THEN    ! Swap positions i & j
              w = a(i)
              ww = t(i)
              a(i) = a(j)
              t(i) = t(j)
              a(j) = w
              t(j) = ww
              i = i + 1
              j = j - 1
              IF (i.GT.j) EXIT
              ELSE
              EXIT
              END IF
              END DO
              IF (j-l.GE.r-i) THEN
              IF (l.LT.j) THEN
              s = s + 1
              stackl(s) = l
              stackr(s) = j
              END IF
              l = i
              ELSE
              IF (i.LT.r) THEN
              s = s + 1
              stackl(s) = i
              stackr(s) = r
              END IF
              r = j
              END IF
              IF (l.LT.r) GO TO 20
              IF (s.NE.0) GO TO 10
              RETURN
              END SUBROUTINE qsort

end program 
!*****************************************************************************80  
!
! MAIN is the main program ( >> mpirun -np 8 ./EvoluAlgo )
!
!-------------------------------------------------------------------------------  
program main

    use precis_mod
    use mod_tools
    use mod_tree
    use mod_evolu
    use mpi     

    implicit none
        
    integer                     :: i, j, k, l, rank_meta
    integer                     :: generation, Ncurrent_archive, Nmeta_population
    integer                     :: Ninitial, ios
    logical                     :: dominated
    character(255)              :: cfile
    character(1512)             :: buffer
    real(kind=working_precis)   :: pi = acos(-1.0), rndom 
    real(kind=working_precis)   :: t0, t1, portion_archive
    integer, dimension(:), allocatable                     :: rank, rank_mating 
    integer, dimension(:), allocatable                     :: rank_init, ind, permut
    real(kind=working_precis), dimension(:,:), allocatable :: population, archive, init
    real(kind=working_precis), dimension(:,:), allocatable :: meta_population
    real(kind=working_precis), dimension(:,:), allocatable :: meta_performance
    real(kind=working_precis), dimension(:,:), allocatable :: performance, mating
    real(kind=working_precis), dimension(:,:), allocatable :: myperformance, myperformance_init
    real(kind=working_precis), dimension(:,:), allocatable :: archinit, performance_init
    real(kind=working_precis), dimension(:), allocatable   :: I1,I2,child,Itemp

!   Random numbers
    integer :: n_rnd 
    integer,allocatable,dimension(:) :: last_seed   

!   MPI related
    integer :: ierr, taskid, numtasks, lenMPI, islave
    character(MPI_MAX_PROCESSOR_NAME) hostname

!   Input parameters (from the file 'Evolution.ini')
    integer                     :: Nevolution = 10000
    integer                     :: Npopulation = 32
    integer                     :: Narchive = 80
    integer                     :: Ninit = 1000
    integer                     :: Nsteps = 50
    integer                     :: Nruns = 4
    real(kind=working_precis)   :: crossover_rate  = 0.1_working_precis
    real(kind=working_precis)   :: Pcrossover = 0.5_working_precis
    real(kind=working_precis)   :: Pmutation_evolu = 0.05_working_precis
    real(kind=working_precis)   :: STDmutation_evolu = 0.01_working_precis
    character(len=30)           :: save_file_evolu = 'ZZZresults', F_param
    integer                     :: save_rate_evolu = 1
    logical                     :: saveParetoOnly = .false.
    character(len=30)           :: BestDir = ''
    character(len=30)           :: read_file = 'YYY_arch00000000.dat'
    integer                     :: i_read = 1 
    logical                     :: FromEvoluAlgo = .false.

    namelist /param2/ Nevolution, Npopulation, Narchive, Ninit, Nsteps, &
        Nruns, crossover_rate, Pcrossover, Pmutation_evolu, Pmutation_evolu, &
        STDmutation_evolu, save_file_evolu, save_rate_evolu, &
        saveParetoOnly, BestDir, read_file, i_read, FromEvoluAlgo

!*****************************************************************************
!   Reading parameters from Evolution.ini
    open(unit=99, file='Evolution.ini', status='OLD')
    read(99, NML=param2)
    close(99)

    F_param = trim(save_file_evolu) // '_param.dat'
    open(unit=88, file=F_param, status="replace", position="rewind")
    write(88, NML=param2)

!   Allocations
    allocate(population(Npopulation,Ndof), archive(Narchive+Npopulation,Ndof+Neval))
    allocate(archinit(100*Ninit,Ndof+Neval))
    allocate(init(Ninit,Ndof),performance_init(Ninit,Neval))
    allocate(meta_population(Npopulation+Narchive,Ndof))
    allocate(meta_performance(Npopulation+Narchive,Neval))
    allocate(performance(Npopulation,Neval), mating(Npopulation,Ndof))
    allocate(myperformance(Npopulation,Neval), myperformance_init(Ninit,Neval))
    allocate(rank(Npopulation),rank_mating(Npopulation),rank_init(Ninit),ind(Ninit))
    allocate(I1(Ndof),I2(Ndof),Itemp(Ndof),child(Ndof))
    allocate(permut(Npopulation+Narchive))
    call cpu_time(t0)
    call random_seed(SIZE=n_rnd)
    allocate(last_seed(n_rnd))

!   MPI initialization
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, taskid, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)
    call MPI_GET_PROCESSOR_NAME(hostname, lenMPI, ierr)
    call init_random_seed(taskid) 
!     print *, 'taskid', taskid

!   Initial population and its evaluation (MPI)
    if (taskid == 0)  then
!       Random initial population
        if ( len(trim(BestDir)) > 0 ) then
            print *, 'Directory:', BestDir
            call system('ls -1 -p '//TRIM(BestDir)//'/*.dat > list')
            open(15,file="list")
            i = 1
            do
                read(15,'(A)',END=10) cfile
                print*, 'File: ', TRIM(cfile)
                open(16,file=cfile)
                ios = 0
                do while (ios==0)
                    read(16, '(A)', iostat=ios) buffer
                    if (ios ==0) then
                        read(buffer, *, iostat=ios) archinit(i,:)
                        i = i+1
                    end if
                end do    
                close(16)
            end do
10          continue
            print *,'Archive datas loaded'
            Ninitial = i-1
            if (Ninitial > Ninit) then
                call clean_archive2(archinit,Ninit,Ninitial,Ndof,Neval)
                init(1:Ninit,:) = archinit(1:Ninit,1:Ndof)
            else if (Ninitial < Ninit) then
                init(1:Ninitial,1:Ndof) = archinit(1:Ninitial,1:Ndof)
                call random_number( init(Ninitial+1:Ninit,:) )
            else if (Ninitial == Ninit) then
                init(1:Ninit,:) = archinit(1:Ninit,1:Ndof)
            end if
        else 
            print *, 'Random initial population'
            call random_number(init)
        end if
        close(15)
        deallocate(archinit)
    endif
!   Send population to slaves
    call MPI_BCAST (init,Ninit*Ndof,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

!   Evaluation of each individual
    myperformance = 0
    call random_seed(GET=last_seed)
    do i = taskid + 1, Ninit, numtasks
        call evaluate_individual(init(i,:),Ndof,Neval,myperformance_init(i,:),Nsteps,Nruns)
!         print *, 'i:', i, 'perf:', myperformance_init(i,:)
    end do
    call MPI_REDUCE(myperformance_init,performance_init,Ninit*Neval, &
                    MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    call random_seed(PUT=last_seed)

!   Sorting and creation of the archive
    if (taskid == 0) then
    Ncurrent_archive = 0
    do i = 1, Ninit
        rank_init(i) = -1
        do j = 1, Ninit
            dominated = .true.
            do k = 1, Neval
                if (i.ne.j .and. (performance_init(i,k) > performance_init(j,k) - 1E-12)) &
&                   dominated = .false.
            end do
            if ( dominated ) rank_init(i) = rank_init(i) + 1
        end do
        if (rank_init(i)==0) then
            Ncurrent_archive = Ncurrent_archive + 1
            archive(Ncurrent_archive,1:Ndof) = init(i,:)
            archive(Ncurrent_archive,Ndof+1:Ndof+Neval) = performance_init(i,:)
        end if
    end do
    if (saveParetoOnly) then
        call save_pareto(archive(1:Ncurrent_archive,:),Ncurrent_archive, &
&                            Ndof,Neval,save_file_evolu,0)
    else
         call save_archive(archive(1:Ncurrent_archive,:), &
&                            Ncurrent_archive, Ndof,Neval,save_file_evolu,0)
    end if
    
!   Initial population
    call qsorti(rank_init,ind,Ninit)
    do i = 1,Npopulation   
        population(i,1:Ndof) = init(ind(i),1:Ndof)
    end do

    print *, achar(7)
    call cpu_time(t1)
    write(*,'(a,I6,a,I6,a,F5.1,a,F7.1,a,F8.1,a,F8.1)')                      &
&           ' Generation: '         , 0,                                    &
&           ' / # archive: '        , Ncurrent_archive,                     & 
&           ' / From archive: '     , 100*portion_archive,                  & 
&           '% / CPUtime: '         , t1 - t0,                              &
&           ' / F1: '               , maxval(performance_init(:,1)),        &
&           ' / F2: '               , maxval(performance_init(:,2))
    end if
    
!****************************************
!   Main loop
!****************************************
    do generation = 1, Nevolution
        master_thread1: if (taskid == 0) then        
!       Creation of the mating pool 
!        print *, 'mating pool'
        call mating_pool(population,archive,Npopulation,Ncurrent_archive, &
&                        Ndof,Neval,rank,mating,rank_mating,portion_archive)

!       Selection, crossover and mutation
!        print *, 'STD mutation:', STDmutation_evolu
        do i = 1, Npopulation
            call random_number(rndom)
            if (rndom < Pcrossover) then
                call tournament(mating,rank_mating,Npopulation,Ndof,I1)
                call tournament(mating,rank_mating,Npopulation,Ndof,I2)
                call cross_over(I1,I2,Ndof,crossover_rate,Itemp)
                call mutation(Itemp,Ndof,0.5*STDmutation_evolu,Pmutation_evolu,child)
            else
                call tournament(mating,rank_mating,Npopulation,Ndof,Itemp)
                call mutation(Itemp,Ndof,STDmutation_evolu,Pmutation_evolu,child)
            end if
            population(i,:) = child
        end do

!       Enforce bounds and constraints
        do i = 1, Npopulation
            do j = 1,Ndof
                if (population(i,j)<0) population(i,j) = 0.
                if (population(i,j)>1) population(i,j) = 1.
            end do
        end do
        end if master_thread1
        call MPI_BCAST (population,Npopulation*Ndof,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

!       Evaluation of each individual
        call random_seed(GET=last_seed)
        myperformance = 0
        do i = taskid + 1, Npopulation, numtasks
            call evaluate_individual(population(i,:),Ndof,Neval,myperformance(i,:),Nsteps,Nruns)
        end do
        call MPI_REDUCE(myperformance,performance,Npopulation*Neval, &
                        MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call random_seed(PUT=last_seed)

        master_thread2: if (taskid == 0) then        
!       Creation of a meta population (population + archive)
!         print *, 'Creation of a meta population'
        meta_population (1:Npopulation,1:Ndof)  = population
        meta_performance(1:Npopulation,1:Neval) = performance
        meta_population (Npopulation+1:Npopulation+Ncurrent_archive,1:Ndof)  = &
&                archive(1:Ncurrent_archive,1:Ndof)
        meta_performance(Npopulation+1:Npopulation+Ncurrent_archive,1:Neval) = &
&                archive(1:Ncurrent_archive,Ndof+1:Ndof+Neval)
        Nmeta_population = Npopulation + Ncurrent_archive

!       Ranking of the metapopulation and update of the archive
        Ncurrent_archive = 0
        call rperm2(Nmeta_population,permut)
        do l = 1, Nmeta_population
            i = permut(l)
            rank_meta = -1
            do j = 1, Nmeta_population
                dominated = .true.
                do k = 1, Neval
                    if (i.ne.j .and. (meta_performance(i,k) > meta_performance(j,k) - 1E-12)) &
                                dominated = .false.
                end do
                if (dominated) rank_meta = rank_meta + 1
            end do
            if (i .le. Npopulation) rank(i) = rank_meta
            if ( (rank_meta==0) .and. (Ncurrent_archive < Narchive+Npopulation) ) then
                Ncurrent_archive = Ncurrent_archive + 1
                archive(Ncurrent_archive,1:Ndof) = meta_population(i,:)
                archive(Ncurrent_archive,Ndof+1:Ndof+Neval) = meta_performance(i,:)
            end if
        end do

!       Cleaning archive from doubles
        call archive_cleaning(archive,Ncurrent_archive,Ndof,Neval)

!       Remove excess individuals in the archive        
        if (Ncurrent_archive > Narchive) then
            call clean_archive(archive,Narchive,Ncurrent_archive,Ndof,Neval)
            Ncurrent_archive = Narchive
        end if 

!       Save and display
        write(*,'(a)',advance='no') ' '
        if (mod(generation,save_rate_evolu) == 0) then
            if (saveParetoOnly) then
                call save_pareto(archive(1:Ncurrent_archive,:),Ncurrent_archive, &
&                            Ndof,Neval,save_file_evolu,generation)
            else
                call save_archive(archive(1:Ncurrent_archive,:), &
&                            Ncurrent_archive, Ndof,Neval,save_file_evolu,generation)
            end if
        end if
        call cpu_time(t1)
        write(*,'( a,I6, a,I6, a,F5.1, a,F7.1, a,F8.1, a,F8.1 )')           &
&           'Generation: '          , generation,                           &
&           ' / # archive: '        , Ncurrent_archive,                     & 
&           ' / From archive: '     , 100*portion_archive,                  & 
&           '% / CPUtime: '         , t1 - t0,                              &
&           ' / F1: '           , maxval(meta_performance(1:Nmeta_population,1)), &
&           ' / F2: '           , maxval(meta_performance(1:Nmeta_population,2))
        end if master_thread2
    end do  

!   Deallocation of arrays
    deallocate(population, archive)
    deallocate(meta_population)
    deallocate(meta_performance)
    deallocate(performance, mating,myperformance)
    deallocate(init,myperformance_init,performance_init)
    deallocate(rank,rank_mating)
    deallocate(I1,I2,Itemp,child)
    deallocate(permut)
    call MPI_FINALIZE(ierr)
end program main

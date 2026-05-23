!*****************************************************************************80  
!
!! MAIN is the main program
!
program main

    use precis_mod
    use mod_tools
    use mod_tree

    implicit none
    real(kind=working_precis)   :: genome(Ndof)
    real(kind=working_precis)   :: performance(Neval)
!*****************************************************************************80  
    integer                     :: i, j, k
    integer                     :: i_run
    integer                     :: generation
    integer                     :: Ntrees
    integer                     :: Nleaves
    integer                     :: n_rnd 
    integer                     :: Nseeds
    integer                     :: Ntwigs
    integer                     :: Npruned
    logical                     :: err
    character(6)                :: number_file
    character(30)               :: save_local
    character(255)              :: cfile
    character(1512)             :: buffer
    real(kind=working_precis)   :: pi = acos(-1.0)
    real(kind=working_precis)   :: t0, t1, OMPtime
    real(kind=working_precis)   :: location(3)
    real(kind=working_precis)   :: angles_flight(2)
    real(kind=working_precis)   :: new_genes(Ndof)
    real(kind=working_precis)   :: elev(Nelev*Nazim), azim(Nelev*Nazim)
    real(kind=working_precis)   :: angle, amplitude
    real(kind=working_precis)   :: VolumeTwig, VolumePerLeaf
    real(kind=working_precis)   :: U_pruning(3)
    real(kind=working_precis)   :: toto_read(6)
    type(tree), pointer         :: t
    type(branch), pointer       :: b
    type(tree_pointer)          :: trees(2)
    logical                     :: trees_logical(2)

!   Allocable
    type(leaf), dimension(:), allocatable                   :: Leaves
    real(kind=working_precis), dimension(:), allocatable    :: rnd_generation
    real(kind=working_precis), dimension(:), allocatable    :: perf_generation
    real(kind=working_precis), dimension(:), allocatable    :: perf_generation2
    real(kind=working_precis), dimension(:,:), allocatable  :: perf_run

!   Input parameters (from the file 'Forest.ini')
    real(kind=working_precis)   :: Pmutation = 0.5_working_precis
    real(kind=working_precis)   :: STDmutation = 0.01_working_precis
    real(kind=working_precis)   :: TwigLength = 1.0_working_precis
    real(kind=working_precis)   :: TwigDiameter = 0.1_working_precis
    real(kind=working_precis)   :: SizeLeaf = 1.0_working_precis
    real(kind=working_precis)   :: LeafSurface = 0.25_working_precis
    real(kind=working_precis)   :: Cauchy = 10000.0_working_precis
    real(kind=working_precis)   :: VolumeRatioLeaf = 4.0_working_precis
    real(kind=working_precis)   :: MaintenanceH = 0.02_working_precis
    integer                     :: Nmax = 10000
    integer                     :: Ngeneration = 10000
    real(kind=working_precis)   :: SizeForest = 100.0_working_precis
    integer                     :: Ntrees_ini = 100
    integer                     :: Ntrees_max = 10000
    character(len=30)           :: save_file = 'ZZZ', F_param 
    integer                     :: save_rate = 1
    character(len=30)           :: tree_init = ''
    real(kind=working_precis)   :: F2min = 100.0_working_precis
    real(kind=working_precis)   :: F2max = 1000.0_working_precis
    logical                     :: FromForest = .true.

    namelist /param1/ Pmutation, STDmutation, TwigLength, TwigDiameter, &
        SizeLeaf, LeafSurface, Cauchy, VolumeRatioLeaf, MaintenanceH, Nmax, &
        Ngeneration, SizeForest, Ntrees_ini, Ntrees_max, save_file, save_rate, &
        tree_init, F2min, F2max, FromForest

!*****************************************************************************
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
    character(len=30)           :: save_file_evolu = 'ZZZresults'
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
!   Reading parameters from Forest.ini
    open(unit=98, file='Forest.ini', status='OLD')
    read(98, NML=param1)
    close(98)
!   Reading parameters from Evolution.ini
    open(unit=99, file='Evolution.ini', status='OLD')
    read(99, NML=param2)
    close(99)

!*****************************************************************************
    ! Reading genome from archive file READ_FILE
    if ( FromEvoluAlgo ) then
        open(17, file=read_file,action="read",position="rewind")
        do i=1,i_read
            read(unit=17, fmt=*) genome, performance
        end do 
        print *, 'performance:', performance
    else
        open(17, file=read_file,action="read",position="rewind")
        do i=1,i_read
            read(unit=17, fmt=*) toto_read, genome
        end do 
        print *, 'i_read', i_read, '(x,y)', toto_read(3:4)
    end if        

!   Allocations + Initialization
    trees_logical = .false.
    allocate(rnd_generation(Nsteps))
    allocate(perf_generation(Nsteps),perf_generation2(Nsteps),perf_run(Neval,Nruns))
    VolumeTwig = 0.25 * pi * TwigLength * TwigDiameter**2
    VolumePerLeaf = VolumeRatioLeaf * VolumeTwig
    perf_run = 0.0
    performance = 0.0

    k = 0
    do i = 1, Nelev
        do j = 1, Nazim
            k = k + 1
            elev(k) = acos( (i-0.5)/Nelev ) 
            azim(k) = pi*2.0*(j-1)/Nazim
        end do
    end do

print *, 'Start runs'
!*****************************************************************************
!  DIFFERENT RUNS
!*****************************************************************************
different_runs: do i_run = 1, Nruns
!     print *, 'run', i_run
!   Initialization
!     call random_set(i_run)
    call init_random_seed(i_run+6)
    call random_number(rnd_generation)
    location = (/ 0.0d0, 0.0d0, 0.0d0 /)
    call new_tree(trees,trees_logical,1,location,0.0d0,genome,TwigLength,TwigDiameter)
    Ntrees = 1
    t => trees(1)%p

    write(number_file,fmt='(a,I4.4)') '1r', i_run
    save_local = trim(save_file) // number_file // 'gen'
    call order_tree(t,err)
    call save_tree(t%branches,t%n_branches,save_local,0,t%Reserve/VolumeTwig)
    call make_statistics(t%n_leaves,t%branches,t%n_branches)
    call save_statistics(t%n_leaves,t%branches,t%n_branches,save_local,0,10)
    call save_area(t,save_local,0)
    open(unit=23, file='ZAllocation.dat', status='UNKNOWN')
    close(unit=23, status="DELETE")

!****************************************
!   Main loop
!****************************************
time_step: do generation = 1, Nsteps
!       Light interception
        call count_leaves(trees,trees_logical,Nleaves)    
        allocate(Leaves(Nleaves))
        call leaves_extract(trees,trees_logical,Leaves)
        do k = 1, Nelev*Nazim
            call light_interception(Leaves,Nleaves,SizeLeaf,elev(k),azim(k),k)
        end do
        call light_on_trees(Leaves,trees,trees_logical)
        deallocate(Leaves)
        t => trees(1)%p

!       Calculating stresses, request biomass, secondary growth
        call calculate_stresses(t, LeafSurface, Cauchy)
        call requested_growth(t, MaintenanceH)
        call secondary_growth(t, VolumePerLeaf)

!       Pruning + deleting
        angle = 1.0 * generation
        amplitude = 0.835 - log( rnd_generation(generation) ) / 6.0 
        U_pruning = (/ amplitude * cos(angle), amplitude * sin(angle), 0.0d0 /)
        Npruned = t%n_branches
        call pruning(t,LeafSurface,Cauchy,U_pruning)
        call order_tree(t,err)
        Npruned = Npruned - t%n_branches

!       New seeds + new leafs
        Ntwigs = ceiling(t%Reserve/VolumeTwig)
        call primary_growth(t,TwigLength,TwigDiameter,generation,Nseeds)
        Ntwigs = Ntwigs - ceiling(t%Reserve/VolumeTwig)
        t%Reserve = t%Reserve - 5.0 * VolumeTwig * Nseeds

!       Ordering
        call order_tree(t,err)
        call count_leaves(trees,trees_logical,Nleaves)    
        allocate(Leaves(Nleaves))
        call leaves_extract(trees,trees_logical,Leaves)
        do k = 1, Nelev*Nazim
            call light_interception(Leaves,Nleaves,SizeLeaf,elev(k),azim(k),k)
        end do
        call light_on_trees(Leaves,trees,trees_logical)
        deallocate(Leaves)

!       Saving
        print *,'generation', generation, 'N branches', t%n_branches, 'N leaves', t%n_leaves
        open(unit=19,file='Z_history.dat', position="append")
        write (19,*) 'generation', generation, 'N branches', t%n_branches, 'N leaves', t%n_leaves
        close(19) 
        call save_tree(t%branches,t%n_branches,save_local,generation,t%Reserve/VolumeTwig)
        call make_statistics(t%n_leaves,t%branches,t%n_branches)
        call save_statistics(t%n_leaves,t%branches,t%n_branches,save_local,generation,10)
        call save_area(t,save_local,generation)
        call save_allocation(t,amplitude,Ntwigs,Nseeds,Npruned,save_local,generation)

!       Evaluation of performance
        perf_generation(generation) = 0.0
        do i = 1, t%n_leaves
            perf_generation(generation) = perf_generation(generation) + &
&               t%leaves(i)%p%location(3) + t%leaves(i)%p%length * t%leaves(i)%p%unit_t(3)
        end do  
        perf_generation2(generation) = Nseeds
!         print *, 'perf_generations', generation, perf_generation(generation), &
!             perf_generation2(generation)
    end do time_step
!****************************************

    do i = Nsteps - 19, Nsteps
        perf_run(1,i_run) = perf_run(1,i_run) + 1.0 / 20.0 * perf_generation(i) 
    end do
    do i = 1, Nsteps
        perf_run(2,i_run) = perf_run(2,i_run) + 1.0 * perf_generation2(i) 
    end do
    ! deallocation of tree
    call delete_tree(trees,trees_logical,1)
    print *, save_local, perf_run(:,i_run)
end do different_runs
!*****************************************************************************

    do i = 1, Nruns
        performance(1) = performance(1) + perf_run(1,i) / Nruns
        performance(2) = performance(2) + perf_run(2,i) / Nruns
    end do
    deallocate(rnd_generation)

    print *, 'performance:', performance
end program main



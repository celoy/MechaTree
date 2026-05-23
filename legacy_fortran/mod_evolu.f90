module mod_evolu

    use precis_mod
    use mod_tools
    use mod_tree
     
!*****************************************************************************80  
CONTAINS 
!*****************************************************************************80  

!*****************************************************************************80  
subroutine mating_pool(population,archive,Npopulation,Narchive,Ndof,Neval,rank, &
    mating,rank_mating,portion_archive)
    implicit none    
    integer, intent(in)                     :: Npopulation,Narchive,Ndof,Neval
    integer, dimension(:), intent(in)       :: rank
    real(kind=working_precis), intent(out)  :: portion_archive
    real(kind=working_precis), dimension(:,:), intent(in)   :: population,archive
    real(kind=working_precis), dimension(:,:), intent(out)  :: mating(Npopulation,Ndof)
    integer, dimension(:), intent(out)                      :: rank_mating(Npopulation)
!*****************************************************************************80  
! MATING_POOL creates a mating pool population from the current population and
!   archive
!-------------------------------------------------------------------------------  
    logical                                     :: same_box
    integer                                     :: i, j, k, l, Nmating
    integer, dimension(Narchive)                :: permut
    integer, dimension(Npopulation)             :: permut2
    integer, dimension(Narchive,Neval)          :: Eval_discrete
    integer, dimension(Npopulation,Neval)       :: Eval_mating
    real(kind=working_precis), dimension(Neval) :: Fmin, Fmax, deltaF 

!   Calculate the size of the boxes on the Pareto front
    do i = 1, Neval
        Fmin(i) = minval(archive(:,Ndof+i))
        Fmax(i) = maxval(archive(:,Ndof+i))
        deltaF(i) = 4.0*(Fmax(i)-Fmin(i))/Npopulation**(1.0/(Neval-1.0))
    end do

    call rperm2(Narchive, permut) 

!   Place each individual of the archive in one box and fill the mating pool
    Nmating = 0
    do k = 1, Narchive
        i = permut(k)
        do j = 1, Neval
            Eval_discrete(i,j) = ceiling((archive(i,Ndof+j)-Fmin(j)) / deltaF(j))
        end do
        ! print *,'individual:', i , ' Eval discrete:', Eval_discrete(i,:)
        same_box = .false.
        do l = 1, Nmating
            same_box = .true.
            do j = 1, Neval
                if (Eval_mating(l,j) .ne. Eval_discrete(i,j)) same_box = .false.
            end do
            if (same_box) exit
        end do
        if (.not.same_box .and. Nmating < Npopulation) then
            Nmating = Nmating +1
            Eval_mating(Nmating,:) = Eval_discrete(i,:)
            mating(Nmating,:) = archive(i,1:Ndof)
            rank_mating(Nmating) = 0 
        end if
    end do
    ! print *, 'N mating pool:', Nmating
 
    portion_archive = 1.0 * Nmating / Npopulation 

!   Complete the mating pool with random individuals from the population
    if (Nmating < Npopulation) then
        call rperm2(Npopulation,permut2)
        do i = Nmating+1, Npopulation
            mating(i,:) = population(permut2(i),:)
            rank_mating(i) = rank(permut2(i))
        end do
    end if  
end subroutine mating_pool

!*****************************************************************************80  
subroutine tournament(mating,rank,Nmating,Ndof,winner)
    implicit none    
    integer, intent(in)                                         :: Nmating,Ndof
    integer, dimension(:), intent(in)                           :: rank
    real(kind=working_precis), dimension(:,:), intent(in)       :: mating
    real(kind=working_precis), dimension(Ndof), intent(out)     :: winner
!*****************************************************************************80  
! TOURNAMENT picks two random indidividuals of the mating pool and select the 
!   winner of the lowest rank
!-------------------------------------------------------------------------------  
    integer                     :: i1, i2
    real(kind=working_precis)   :: rndom(2)

    call random_number(rndom)
    i1 = ceiling(Nmating*rndom(1))
    i2 = ceiling(Nmating*rndom(2))

    if (rank(i1) < rank(i2)) then
        winner(:) = mating(i1,:)
    else 
        winner(:) = mating(i2,:)
    end if
end subroutine tournament

!*****************************************************************************80  
subroutine cross_over(I1,I2,Ndof,cross_rate,child)
    implicit none    
    integer, intent(in)                                       :: Ndof
    real(kind=working_precis), intent(in)                     :: cross_rate
    real(kind=working_precis), dimension(:), intent(in)       :: I1,I2
    real(kind=working_precis), dimension(Ndof), intent(out)   :: child
!*****************************************************************************80  
! CROSSOVER calculates the child from two parents I1 and I2
!-------------------------------------------------------------------------------  
    integer                     :: i
    real(kind=working_precis)   :: rndom(Ndof), cursor

    call random_number(cursor)
    call random_number(rndom)
    do i = 1, Ndof
        if (rndom(i) < cross_rate) then
            child(i) = I1(i) + cursor * (I2(i) - I1(i))
        else
            child(i) = I1(i)
        end if
    end do
end subroutine cross_over

!*****************************************************************************80  
subroutine mutation(parent,Ndof,STDmutation,Pmutation,child)
    implicit none    
    integer, intent(in)                                       :: Ndof
    real(kind=working_precis), intent(in)                     :: STDmutation,Pmutation
    real(kind=working_precis), dimension(:), intent(in)       :: parent
    real(kind=working_precis), dimension(Ndof), intent(out)   :: child
!*****************************************************************************80  
! MUTATION mutates the PARENT into a CHILD
!-------------------------------------------------------------------------------  
    integer                     :: i
    real(kind=working_precis)   :: rndom(Ndof)
    logical                     :: change

    change = .true. 
    call random_number(rndom)
    do i = 1, Ndof
        if (rndom(i) < Pmutation) then
            child(i) = parent(i) + STDmutation*random_normal()
            child(i) = max(0.0,min(1.0,child(i)))
            change = .false.
        else
            child(i) = parent(i)
        end if
    end do
    i = ceiling( rndom(Ndof) * Ndof )
    i = max(1, min(Ndof, i)) 
    if (change) then
        child(i) = parent(i) + STDmutation*random_normal() 
        child(i) = max(0.0,min(1.0,child(i)))
    end if
end subroutine mutation

!*****************************************************************************80  
subroutine archive_cleaning(archive,Narchive,Ndof,Neval)
    implicit none    
    integer, intent(inout)                                    :: Narchive
    integer, intent(in)                                       :: Ndof, Neval
    real(kind=working_precis), dimension(:,:), intent(inout)  :: archive
!*****************************************************************************80  
! ARCHIVE_CLEANING removes the double of individuals in the archive 
!-------------------------------------------------------------------------------
    integer                   :: i, j, N
    real(kind=working_precis) :: distance(Narchive)

    N = 2
    do while (N < Narchive)
        do i = 1, N-1
            distance(i) = norm2( archive(N,Ndof+1:Ndof+Neval) - archive(i,Ndof+1:Ndof+Neval) ) + &
                          norm2( archive(N,1:Ndof/2)          - archive(i,1:Ndof/2) ) 
        end do
        if ( minval(distance(1:N-1)) < 1E-6 ) then
            do j = N, Narchive - 1
                archive(j,:) = archive(j+1,:)
            end do
            Narchive = Narchive - 1
            ! print *, 'doublon'
        else
            N = N + 1
        end if
    end do
end subroutine

!*****************************************************************************80  
subroutine clean_archive(archive,Narchive,Ncurrent,Ndof,Neval)
    implicit none    
    integer, intent(in)                                       :: Narchive,Ndof
    integer, intent(in)                                       :: Neval,Ncurrent
    real(kind=working_precis), dimension(:,:), intent(inout)  :: archive
!*****************************************************************************80  
! CLEAN_ARCHIVE removes the excess of individuals in the archive to have exactly
!   Narchive individuals
!-------------------------------------------------------------------------------
    integer                   :: i, j, k, N, minI, maxAbs(Neval) 
    real(kind=working_precis) :: distance(Ncurrent), mindistance(Ncurrent)
    logical                   :: MaxEval

    N = Ncurrent
    do while (N > Narchive)
        do i = 1, Neval
            maxAbs(i) = maxloc(archive(1:N,Ndof+i),1)
        end do
        do i = 1, N
            do j = 1, N
                distance(j) = norm2( archive(i,Ndof+1:Ndof+Neval) - archive(j,Ndof+1:Ndof+Neval) )
            end do
            distance(i) = 1E16
            mindistance(i) = minval(distance)
            MaxEval = .false.
            do k = 1, Neval
                if (i==maxAbs(k)) MaxEval=.true.
            end do 
            if (MaxEval) mindistance(i) = 1E16
        end do
        minI = minloc(mindistance(1:N),1)
        archive(minI,:) = archive(N,:)
        N = N - 1
    end do
end subroutine

!*****************************************************************************80  
subroutine clean_archive2(archive,Narchive,Ncurrent,Ndof,Neval)
    implicit none    
    integer, intent(in)                                       :: Narchive,Ndof
    integer, intent(in)                                       :: Neval,Ncurrent
    real(kind=working_precis), dimension(:,:), intent(inout)  :: archive
!*****************************************************************************80  
! CLEAN_ARCHIVE removes the excess of individuals in the archive to have exactly
!   Narchive individuals
!-------------------------------------------------------------------------------
    integer :: i,j,N,minI
    real(kind=working_precis) :: distance(Ncurrent), mindistance(Ncurrent)

    N = Ncurrent
    do while (N > Narchive)
        do i = 1, N
            do j = 1, N
                distance(j) = norm2( archive(i,Ndof+1:Ndof+Neval) - archive(j,1:Ndof+1:Ndof+Neval) )
            end do
            distance(i) = 1E16
            mindistance(i) = minval(distance)
        end do
        minI = minloc(mindistance,1)
        archive(minI,:) = archive(N,:)
        N = N - 1
    end do
end subroutine

!*****************************************************************************80  
subroutine save_pareto(archive,Narchive,Ndof,Neval,save_file,iout)
    implicit none    
    integer, intent(in)                                       :: Narchive,Ndof
    integer, intent(in)                                       :: Neval,iout
    real(kind=working_precis), dimension(:,:), intent(in)     :: archive
    character(len=30), intent(in)                             :: save_file
!*****************************************************************************80  
! SAVE_PARETO saves the coordinates of the individuals of the archive
!-------------------------------------------------------------------------------  
    integer                     :: i, j
    character(8)                :: number_file

    write(number_file,fmt='(I8.8)') iout
    open(unit=10,file=trim(save_file) // number_file // '.dat', &
&     status="replace", position="rewind")
    
    do i = 1, Narchive
        do j = 1, Neval
            write(10,fmt='(F16.8)',advance='no') archive(i,Ndof+j)
        end do
        write(10,*)
    end do
    close( unit=10 )
end subroutine save_pareto

!*****************************************************************************80  
subroutine save_archive(archive,Narchive,Ndof,Neval,save_file,iout)
    implicit none    
    integer, intent(in)                                       :: Narchive,Ndof
    integer, intent(in)                                       :: Neval,iout
    real(kind=working_precis), dimension(:,:), intent(in)     :: archive
    character(len=30), intent(in)                             :: save_file
!*****************************************************************************80  
! SAVE_ARCHIVE saves the entire archive
!-------------------------------------------------------------------------------  
    integer                     :: i
    character(8)                :: number_file

    write(number_file,fmt='(I8.8)') iout
    open(unit=11,file=trim(save_file) // '_arch' // number_file // '.dat', &
&     status="replace", position="rewind")

    do i = 1, Narchive
            write(11,*) archive(i,1:Ndof+Neval)
    end do
    close( unit=11 )
end subroutine save_archive

!*****************************************************************************80  
subroutine save_forest(trees,Ntrees,save_file,iout)
    implicit none    
    type(tree_pointer), dimension(:), intent(in)    :: trees
    integer, intent(in)                             :: Ntrees, iout
    character(len=30), intent(in)                   :: save_file
!*****************************************************************************80  
!
!! SAVE_FOREST save the data of the current tree in the file:
!            SAVE_FILE // I // '.dat'
!
    integer      :: i,j
    character(8) :: number_file    

    write(number_file,fmt='(I8.8)') iout
    open(unit=13,file=trim(save_file) // number_file // '.dat', &
&     status="replace", position="rewind")
!     print *, save_file, trim(save_file)

    do i = 1, Ntrees
        do j = 1, trees(i)%p%n_branches
            write(13,*) &
&                   trees(i)%p%branches(j)%p%generation, &  ! 0
&                   trees(i)%p%branches(j)%p%unit_t,     &  ! 1-3
&                   trees(i)%p%branches(j)%p%unit_b,     &  ! 4-6
&                   trees(i)%p%branches(j)%p%location,   &  ! 7-9
&                   trees(i)%p%branches(j)%p%light,      &  ! 10
&                   trees(i)%p%branches(j)%p%diameter,   &  ! 11
&                   i,                                   &  ! 12
&                   trees(i)%p%branches(j)%p%Strahler,   &  ! 13
&                   trees(i)%p%branches(j)%p%botanic        ! 14
        end do
    end do
    close(13)
end subroutine save_forest

!*****************************************************************************80  
subroutine save_simplified(trees,Ntrees,save_file,iout)
    implicit none    
    type(tree_pointer), dimension(:), intent(in)    :: trees
    integer, intent(in)                             :: Ntrees, iout
    character(len=30), intent(in)                   :: save_file
!*****************************************************************************80  
!
!! SAVE_FOREST save the data of the current tree in the file:
!            SAVE_FILE // I // '.dat'
!
    integer                     :: i, j
    character(8)                :: number_file    
    real(kind=working_precis)   :: max_height, X, Y, Xmean, Ymean, Rstd
    real(kind=working_precis)   :: max_H, biomass, pi = acos(-1.0)
    integer                     :: max_branches

    write(number_file,fmt='(I8.8)') iout
    open(unit=13,file=trim(save_file) // number_file // '.dat', &
&     status="replace", position="rewind")

    max_H = 0.0
    max_branches = 0
    do i = 1, Ntrees
        Xmean = 0.0
        Ymean = 0.0
        Rstd  = 0.0
        max_height = 0.0
        biomass    = 0.0
        X = trees(i)%p%trunk%location(1)
        Y = trees(i)%p%trunk%location(2)
        do j = 1, trees(i)%p%n_leaves
            Xmean = Xmean + trees(i)%p%leaves(j)%p%location(1) / trees(i)%p%n_leaves
            Ymean = Ymean + trees(i)%p%leaves(j)%p%location(2) / trees(i)%p%n_leaves
            if ( trees(i)%p%leaves(j)%p%location(3) > max_height ) then
                max_height = trees(i)%p%leaves(j)%p%location(3)
            end if
        end do
        do j = 1, trees(i)%p%n_leaves
            Rstd =  ( trees(i)%p%leaves(j)%p%location(1) - Xmean )**2.0 + &
                    ( trees(i)%p%leaves(j)%p%location(2) - Ymean )**2.0 + Rstd
        end do
        do j = 1, trees(i)%p%n_branches
            biomass = biomass + trees(i)%p%branches(j)%p%length*0.25*pi*trees(i)%p%branches(j)%p%diameter**2.0
        end do
        Rstd = sqrt(Rstd/trees(i)%p%n_leaves)
        if (trees(i)%p%n_branches > max_branches) max_branches = trees(i)%p%n_branches
        if (max_height > max_H) max_H = max_height
        write(13,*) X, Y, Xmean, Ymean, Rstd, max_height, trees(i)%p%genes, &
                    trees(i)%p%trunk%diameter, trees(i)%p%n_leaves, trees(i)%p%n_branches, biomass
    end do
    close(13)
    print *, 'Gen:', iout, 'Ntrees:', Ntrees, 'max_H:', max_H, 'max_branches', max_branches
end subroutine save_simplified

















!*****************************************************************************80  
!
! EVALUATE_INDIVIDUAL evaluates the performances of an individual of genome
!   GENOME and put the results into the array PERFORMANCE
!
!-------------------------------------------------------------------------------  
subroutine evaluate_individual(genome,Ndof,Neval,performance,Nsteps,Nruns)
    implicit none
    integer, intent(in)                     :: Ndof, Neval
    real(kind=working_precis), intent(in)   :: genome(Ndof)
    real(kind=working_precis), intent(out)  :: performance(Neval)
    integer, intent(in)                     :: Nsteps, Nruns
!*****************************************************************************80  
    integer                     :: i, j, k
    integer                     :: i_run
    integer                     :: generation
    integer                     :: Ntrees
    integer                     :: Nleaves
    integer                     :: n_rnd 
    integer                     :: Nseeds
    logical                     :: err
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
    real(kind=working_precis)   :: perf_generation(Nsteps)
    real(kind=working_precis)   :: perf_generation2(Nsteps)
    real(kind=working_precis)   :: perf_run(Neval,Nruns)
    type(tree), pointer         :: t
    type(branch), pointer       :: b
    type(tree_pointer)          :: trees(2)
    logical                     :: trees_logical(2)

!   Allocable
    type(leaf), dimension(:), allocatable                   :: Leaves
    real(kind=working_precis), dimension(:), allocatable    :: rnd_generation

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

!     print *, 'START EVALUATION'
!*****************************************************************************
!   Reading parameters from Forest.ini
    open(unit=99, file='Forest.ini', status='OLD')
    read(99, NML=param1)
    close(99)

!   Allocations + Initialization
    trees_logical = .false.
    allocate(rnd_generation(Nsteps))
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

!*****************************************************************************
!  DIFFERENT RUNS
!*****************************************************************************
different_runs: do i_run = 1, Nruns
!     print *, 'run', i_run
! 	Initialization
    call random_set(i_run + 0)
    call random_number(rnd_generation)
    location = (/ 0.0d0, 0.0d0, 0.0d0 /)
    call new_tree(trees,trees_logical,1,location,0.0d0,genome,TwigLength,TwigDiameter)
    Ntrees = 1

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
        call pruning(t,LeafSurface,Cauchy,U_pruning)
        call order_tree(t,err)

!       New seeds + new leafs
        call primary_growth(t,TwigLength,TwigDiameter,generation,Nseeds)
        t%Reserve = t%Reserve - 5.0 * VolumeTwig * Nseeds

! 		Ordering
        call order_tree(t,err)

! 		Evaluation of performance
        perf_generation(generation) = 0.0
        do i = 1, t%n_leaves
            perf_generation(generation) = perf_generation(generation) + &
&               t%leaves(i)%p%location(3) + t%leaves(i)%p%length * t%leaves(i)%p%unit_t(3)
        end do  
        perf_generation2(generation) = Nseeds
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
end do different_runs
!*****************************************************************************

    do i = 1, Nruns
        performance(1) = performance(1) + perf_run(1,i) / Nruns
        performance(2) = performance(2) + perf_run(2,i) / Nruns
    end do
    deallocate(rnd_generation)
end subroutine evaluate_individual

end module mod_evolu
program main

    use precis_mod
    use mod_tools
    use mod_tree
    use mod_evolu
    use omp_lib       

    implicit none
        
    integer                     :: i, j, k
    integer                     :: generation
    integer                     :: Ntrees, Nbigtrees
    integer                     :: Nleaves
    integer                     :: n_rnd 
    integer                     :: Nseeds
    integer                     :: Ninitial, ios
    integer                     :: Nbmax, ibmax
    logical                     :: err
    character(255)              :: cfile
    character(1512)             :: buffer
    character(30)               :: save_local = 'Zone'
    real(kind=working_precis)   :: pi = acos(-1.0)
    real(kind=working_precis)   :: t0, t1, OMPtime
    real(kind=working_precis)   :: location(3)
    real(kind=working_precis)   :: angles_flight(2)
    real(kind=working_precis)   :: new_genes(Ndof)
    real(kind=working_precis)   :: elev(Nelev*Nazim), azim(Nelev*Nazim)
    real(kind=working_precis)   :: angle, amplitude
    real(kind=working_precis)   :: VolumeTwig, VolumePerLeaf
    real(kind=working_precis)   :: U_pruning(3)
    real(kind=working_precis)   :: biomass_tot, biomass_square_tot
    type(tree), pointer         :: t
    type(branch), pointer       :: b

!   Allocable
    type(tree_pointer), dimension(:), allocatable           :: trees
    logical, dimension(:), allocatable                      :: trees_logical
    real(kind=working_precis), dimension(:), allocatable    :: Xposition_ini
    real(kind=working_precis), dimension(:), allocatable    :: Yposition_ini
    real(kind=working_precis), dimension(:), allocatable    :: angle_ini   
    real(kind=working_precis), dimension(:), allocatable    :: rand_leaves
    real(kind=working_precis), dimension(:), allocatable    :: rnd_generation
    real(kind=working_precis), dimension(:,:), allocatable  :: genome_ini 
    real(kind=working_precis), dimension(:,:), allocatable  :: archinit    
    real(kind=working_precis), dimension(:,:), allocatable  :: bestforest    
    type(leaf), dimension(:), allocatable                   :: Leaves

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
!   Reading parameters from Forest.ini
    open(unit=99, file='Forest.ini', status='OLD')
    read(99, NML=param1)
    close(99)

    F_param = trim(Save_file) // '_param.dat'
    open(unit=88, file=F_param, status="replace", position="rewind")
    write(88, NML=param1)

!   Initialization
    call cpu_time(t0)
!     OMPtime = OMP_get_wtime()
!     call OMP_SET_NUM_THREADS(8)
    call random_seed(SIZE=n_rnd)
    call init_random_seed(0) 
!     print *, '# of OMP threads:', OMP_get_num_threads()
    VolumeTwig = 0.25 * pi * TwigLength * TwigDiameter**2
    VolumePerLeaf = VolumeRatioLeaf * VolumeTwig
    open(unit=23, file='SelfThinning.dat', status='UNKNOWN')
    close(unit=23, status="DELETE")

!   Allocations
    allocate(trees(Ntrees_max))
    allocate(trees_logical(Ntrees_max))
    allocate(rand_leaves(Nmax))
    trees_logical = .false.
    allocate(Xposition_ini(Ntrees_ini),Yposition_ini(Ntrees_ini),angle_ini(Ntrees_ini))
    allocate(genome_ini(Ntrees_ini,Ndof))
    allocate(rnd_generation(Ngeneration))
    allocate(archinit(Ntrees_ini,Ndof+Neval))
    allocate(bestforest(Ntrees_ini,6+Ndof))
    call random_number(Xposition_ini)
    call random_number(Yposition_ini)
    call random_number(angle_ini)
    call random_number(genome_ini)
    call random_number(rnd_generation)

    k = 0
    do i = 1, Nelev
        do j = 1, Nazim
            k = k + 1
            elev(k) = acos( (i-0.5)/Nelev ) 
            azim(k) = pi*2.0*(j-1)/Nazim
        end do
    end do

!   Initial population
    if ( len(trim(tree_init)) > 0 ) then
        print *, 'Directory:', tree_init
        call system('ls -1 -p '//TRIM(tree_init)//'/*.dat > list')
        open(15,file="list")
        i = 1
        do
            j = 0
            read(15,'(A)',END=10) cfile
            print*, 'File: ', TRIM(cfile)
            open(16,file=cfile)
            ios = 0
!             do while (ios==0 .and. i.le.Ntrees_ini)
            do while (ios==0 .and. i.le.Ntrees_ini .and. j < 2000)
                read(16, '(A)', iostat=ios) buffer
                if (ios ==0) then
                    if ( FromForest ) then
                        read(buffer, *, iostat=ios) bestforest(i,:)
                        i = i + 1
                        j = j + 1
                    else
                        read(buffer, *, iostat=ios) archinit(i,:)
                        if ( archinit(i,Ndof+2) < F2max .and. archinit(i,Ndof+2) > F2min ) then
                            print *,'Neval' , archinit(i,Ndof+1:Ndof+2), 'i', i
                            i = i + 1
                        end if
                    end if
                end if
            end do    
            close(16)
        end do
10      continue
        Ninitial = i-1
        print *,'Archive datas loaded', Ninitial
        if (Ninitial < Ntrees_ini) then
            if ( FromForest ) then
                do i = 1, Ntrees_ini
                    genome_ini(i,1:Ndof) = bestforest(1 + mod(i,Ninitial), 7:6+Ndof)
                end do
!                 genome_ini(1:Ninitial,1:Ndof) = bestforest(1:Ninitial,7:6+Ndof)
            else
                genome_ini(1:Ninitial,1:Ndof) = archinit(1:Ninitial,1:Ndof)
            end if
        else
            if ( FromForest ) then
                genome_ini(1:Ninitial,1:Ndof) = bestforest(1:Ntrees_ini,7:6+Ndof)
            else
                genome_ini(1:Ninitial,1:Ndof) = archinit(1:Ntrees_ini,1:Ndof)
            end if
        end if
    end if
    close(15)

!   New trees
    do i = 1, Ntrees_ini
        location = (/ SizeForest * sqrt(Xposition_ini(i)) * cos( Yposition_ini(i)*pi*2.0 ), & 
                      SizeForest * sqrt(Xposition_ini(i)) * sin( Yposition_ini(i)*pi*2.0 ), &  
                      0.0_working_precis /)
!         location = 1.77 * SizeForest / sqrt(1.0*Ntrees_ini) * &
!             (/ 1.0 * mod(i-1, int(sqrt(1.0*Ntrees_ini))),  &
!                1.0 * int((i-1)/sqrt(1.0*Ntrees_ini)), &
!                0.0/)
        call new_tree(trees,trees_logical,Ntrees_max,location, &
                      2.0*pi*angle_ini(i), genome_ini(i,:),TwigLength,TwigDiameter)
    end do
    Ntrees = Ntrees_ini
    call save_simplified(trees,Ntrees,save_file,0)
!     call save_forest(trees,Ntrees,save_file,0)

!*****************************************************************************
!   Main loop
!*****************************************************************************
    do generation = 1, Ngeneration
!       Light interception
!         print *, 'Light interception'
        call count_leaves(trees,trees_logical,Nleaves)    
        allocate(Leaves(Nleaves))
!         print *, 'Light interception 1', Nleaves
        call leaves_extract(trees,trees_logical,Leaves)
!$OMP PARALLEL DO
        do k = 1, Nelev*Nazim
            call light_interception(Leaves,Nleaves,SizeLeaf,elev(k),azim(k),k)
        end do
!$OMP END PARALLEL DO
        call light_on_trees(Leaves,trees,trees_logical)
        deallocate(Leaves)

!       Save  
!         print *, 'Saving'
        Nbmax = 0
        ibmax = 0
        do i = 1, Ntrees
            t => trees(i)%p
            if (generation > 999) call make_statistics(t%n_leaves,t%branches,t%n_branches)
            if ( (t%n_branches > Nbmax) .and. (norm2(t%trunk%location(1:2)) < 0.5 * SizeForest) ) then
                Nbmax = trees(i)%p%n_branches
                ibmax = i
            end if 
        end do
        t => trees(ibmax)%p
        print *,'generation', generation, 'N branches', t%n_branches, 'N trees', Ntrees
        open(unit=19,file='Z_history.dat', position="append")
!         print *,'generation', generation, 'location', t%trunk%location(1:2), &
!             'DoB', t%branches(2)%p%generation, 'N branches', t%n_branches
!         open(unit=19,file='Z_history.dat', position="append")        
!         write (19,*) 'generation', generation, 'location', t%trunk%location(1:2), &
!             'DoB', t%branches(2)%p%generation, 'N branches', t%n_branches
!         close(19)    

        if (mod(generation,save_rate)==0 .and. generation > 999) then
            call save_statistics(t%n_leaves,t%branches,t%n_branches,save_local,generation,10)
            call save_area(t,save_local,generation)
            call save_simplified(trees,Ntrees,save_file,generation)
!             call save_forest(trees,Ntrees,save_file,generation)
        end if   
        call cpu_time(t1)
        if (Nleaves < 1) exit

        biomass_tot = 0.0
        biomass_square_tot = 0.0
        Nbigtrees = 0
        do i = 1, Ntrees
            call biomass_calculation(trees(i)%p,biomass_tot,biomass_square_tot)
            if ( trees(i)%p%n_branches > 10 ) Nbigtrees = Nbigtrees + 1
        end do
        open(unit=23, file='SelfThinning.dat', status='UNKNOWN', position='APPEND')
        write(23,*) generation, Ntrees, Nbigtrees, &
 &                  biomass_tot**2.0/biomass_square_tot, &
 &                  biomass_square_tot/biomass_tot
        close(23)

!       Calculating stresses, request biomass, secondary growth
!         print *, 'Secondary Growth'
!$OMP PARALLEL DO
        do i = 1, Ntrees
            trees(i)%p%age = trees(i)%p%age + 1
            call calculate_stresses(trees(i)%p, LeafSurface, Cauchy)
            call requested_growth(trees(i)%p, MaintenanceH)
            call secondary_growth(trees(i)%p, VolumePerLeaf)
        end do
!$OMP END PARALLEL DO

!       Pruning + deleting
!         print *, 'Pruning'
        angle = 1.0 * generation
        amplitude = 0.835 - log( rnd_generation(generation) ) / 6.0 
        U_pruning = (/ amplitude * cos(angle), amplitude * sin(angle), 0.0d0 /)
!$OMP PARALLEL DO
        do i = 1, Ntrees
            trees(i)%p%n_max = max(trees(i)%p%n_branches, trees(i)%p%n_max)
            call pruning(trees(i)%p,LeafSurface,Cauchy,U_pruning)
            call order_tree(trees(i)%p,err)
        end do 
!$OMP END PARALLEL DO
!         print *,  'Generation: ', generation, '  / Ntrees:', Ntrees, &
!                   '   / Nleaves:', Nleaves,   '  / wind:', amplitude                         
            
!         print *, 'Deleting'
        do i = 1, Ntrees
            t => trees(i)%p
            if ( (t%n_branches < 11 .and. t%age > 5) .or. &
                 (t%age > 1000) ) then
!                 print *, 'delete', t%n_branches, t%n_max, t%age
                call delete_tree(trees,trees_logical,i)
            end if
        end do

!         print *, 'Ordering'
        k = 0
        do i = 1, Ntrees
            if (trees_logical(i)) then
                k = k + 1
                trees(k)%p => trees(i)%p
                trees_logical(k) = .true.
            end if
        end do 
        do i = k+1, Ntrees
            trees_logical(i) = .false.
        end do 
        Ntrees = k

!       New seeds + new leafs
!         print *, 'New seeds and leaves'
        do i = 1, Ntrees
            t => trees(i)%p
            call primary_growth(t,TwigLength,TwigDiameter,generation,Nseeds)
            t%Reserve = t%Reserve - 5.0 * VolumeTwig * Nseeds
!             if (Nseeds > 10) print *, 'Nseeds', Nseeds, 'Reserve', t%Reserve/VolumeTwig, t%n_branches
            Nseeds = min(Nseeds/1, Ntrees_max - Ntrees) 
            call random_number(rand_leaves(1:Nseeds))
            do j = 1, Nseeds
                k = ceiling(rand_leaves(j)*t%n_leaves)
                b => t%leaves(k)%p 
                call random_number(angles_flight)
                angles_flight = angles_flight * pi * 2.0
                location(1) = b%location(1) + b%location(3)*cos(angles_flight(1))
                location(2) = b%location(2) + b%location(3)*sin(angles_flight(1))
                location(3) = 0.0
                if ( (location(1)**2 + location(2)**2) < SizeForest**2 ) then
                    call mutation(t%genes,Ndof,STDmutation,Pmutation,new_genes)
                    call new_tree(trees,trees_logical,Ntrees_max, &
                                  location,angles_flight(2),new_genes, &
                                  TwigLength,TwigDiameter)
                    Ntrees = Ntrees + 1
                end if 
            end do
        end do

!       Ordering
!         print *, 'Ordering'
!$OMP PARALLEL DO
        do i = 1, Ntrees
            trees(i)%p%n_max = max(trees(i)%p%n_branches, trees(i)%p%n_max)
            call order_tree(trees(i)%p,err)
        end do 
!$OMP END PARALLEL DO        
    end do

!*****************************************************************************
!   Deallocations
    print *, achar(7)
    deallocate(trees,trees_logical)
    deallocate(rand_leaves)
    deallocate(Xposition_ini,Yposition_ini,angle_ini,genome_ini)
    deallocate(rnd_generation)
end program main

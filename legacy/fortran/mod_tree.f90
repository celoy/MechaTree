module mod_tree

    use precis_mod
    use mod_tools

    integer, parameter :: max_tree = 10000
    integer, parameter :: Nelev = 4, Nazim = 8
    integer, parameter :: Ndof = 40, Neval = 2

    type branch                             ! Definition of a new type: branch
        integer                                   :: generation                
        integer                                   :: nb_leaves= 1  
        integer                                   :: Strahler = 0
        integer                                   :: botanic = 0
        integer                                   :: distance_ground               
        real(kind=working_precis)                 :: length   = 1.0
        real(kind=working_precis)                 :: diameter = 0.1
        real(kind=working_precis)                 :: light    = 0.0
        real(kind=working_precis)                 :: sun      = 0.0
        real(kind=working_precis)                 :: stress   = 0.0
        real(kind=working_precis)                 :: distance_leaves
        real(kind=working_precis)                 :: max_stress      
        real(kind=working_precis)                 :: vol_growth      
        real(kind=working_precis)                 :: vol_summed
        real(kind=working_precis)                 :: maintenance_vol      
        real(kind=working_precis), dimension(1:3) :: unit_t   = (/ 0.0, 0.0, 1.0 /)
        real(kind=working_precis), dimension(1:3) :: unit_b   = (/ 1.0, 0.0, 0.0 /)
        real(kind=working_precis), dimension(1:3) :: location = (/ 0.0, 0.0, 0.0 /)
        real(kind=working_precis), dimension(1:3) :: T        = (/ 0.0, 0.0, 0.0 /)
        real(kind=working_precis), dimension(1:3) :: M        = (/ 0.0, 0.0, 0.0 /)
        logical                                   :: marker   = .false.
        type (branch), pointer                    :: parent   => null()
        type (branch), pointer                    :: left     => null()
        type (branch), pointer                    :: right    => null()
    end type branch
    
    type branch_pointer
        type (branch), pointer :: p => null()
    end type branch_pointer

    type leaf
        real, dimension(1:3)   :: location = (/ 0.0, 0.0, 0.0 /)
        integer                :: light(Nelev*Nazim)
    end type leaf

    type tree
        real(kind=working_precis)                   :: genes(Ndof) 
        real(kind=working_precis)                   :: theta1, theta2 
        real(kind=working_precis)                   :: gamma1, gamma2
        real(kind=working_precis)                   :: Reserve        
        real(kind=working_precis)                   :: NNbranch(10)
        real(kind=working_precis)                   :: NNreserve(18)
        integer                                     :: age
        integer                                     :: n_branches, n_leaves
        integer                                     :: n_max
        type (branch), pointer                      :: trunk
        type (branch_pointer), dimension(max_tree)  :: branches, leaves
    end type tree

    type tree_pointer
        type (tree), pointer :: p => null()
    end type tree_pointer

!*******************************************************************************    
    CONTAINS !******************************************************************
!*******************************************************************************

!*****************************************************************************80  
subroutine new_tree (trees,trees_logical,Ntrees_max,location,angle,genome, &
                    TwigLength,TwigDiameter)
    implicit none    
    type(tree_pointer), dimension(:), intent(inout)         :: trees
    logical, dimension(:), intent(inout)                    :: trees_logical
    integer, intent(in)                                     :: Ntrees_max
    real(kind=working_precis), dimension(3), intent(in)     :: location
    real(kind=working_precis), intent(in)                   :: angle
    real(kind=working_precis), intent(in)                   :: TwigLength
    real(kind=working_precis), intent(in)                   :: TwigDiameter
    real(kind=working_precis), dimension(Ndof), intent(in)  :: genome
!*****************************************************************************80  
!
!! NEW_TREE creates and allocates a new tree
!    
    integer                     :: i
    type(tree), pointer         :: t  
    real(kind=working_precis)   :: pi = acos(-1.0)

    i = 1
    do while (trees_logical(i))
        i = i + 1
    end do

    allocate(trees(i)%p)
    trees_logical(i) = .true.
    t => trees(i)%p 
    allocate(t%trunk)
    t%trunk%generation = 1
    t%trunk%diameter = TwigDiameter
    t%branches(1)%p => t%trunk
    t%n_branches = 1
    t%leaves(1)%p => t%trunk
    t%n_leaves = 1
    t%age = 0
    t%n_max = 0
    t%Reserve = 2.0 * ( 0.25 * pi * TwigLength * TwigDiameter**2.0 )

    t%genes  =  genome
    t%theta1 =  genome(1) * pi / 2.0
    t%theta2 = -genome(2) * pi / 2.0
    t%gamma1 =  genome(3) * 2.0*pi
    t%gamma2 =  genome(3) * 2.0*pi
    t%NNbranch = genome(4:13)
    t%NNreserve= genome(14:31)

    t%trunk%location = location
    t%trunk%unit_b = (/ cos(angle), sin(angle), 0.0d0 /)
end subroutine new_tree

!*****************************************************************************80  
subroutine delete_tree (trees,trees_logical,ith)
    implicit none    
    type(tree_pointer), dimension(:), intent(inout)     :: trees
    logical, dimension(:), intent(inout)                :: trees_logical
    integer, intent(in)                                 :: ith
!*****************************************************************************80  
!
!! DELETE_TREE deletes and deallocates the i_th tree
!    
    integer                                          :: i
    type(tree), pointer                              :: t  

    t => trees(ith)%p 

    do i = 1, t%n_branches
        if ( associated(t%branches(i)%p) ) deallocate(t%branches(i)%p)
    end do

    deallocate(trees(ith)%p)
    trees_logical(ith) = .false.
end subroutine delete_tree

!*****************************************************************************80  
subroutine count_leaves(trees,trees_logical,Nleaves)
    implicit none    
    type(tree_pointer), dimension(:), intent(in)    :: trees
    logical, dimension(:), intent(in)               :: trees_logical
    integer, intent(out)                            :: Nleaves
!*****************************************************************************80  
!
!! COUNT_LEAVES count the total number of leaves (on all trees)
!    
    integer                                         :: i

    Nleaves = 0
    i = 1
    do while (trees_logical(i))
        Nleaves = Nleaves + trees(i)%p%n_leaves
        i = i + 1
    end do
end subroutine count_leaves

!*****************************************************************************80  
subroutine biomass_calculation(t,B1,B2)
    implicit none    
    type(tree), pointer, intent(in)             :: t
    real(kind=working_precis), intent(inout)    :: B1 
    real(kind=working_precis), intent(inout)    :: B2
!*****************************************************************************80  
!
!! BIOMASS_CALCULATION calculates the biomass and its second moment of a tree
!    
    integer                   :: i
    real(kind=working_precis) :: biomass
    real(kind=working_precis) :: pi = acos(-1.0)
    type (branch), pointer    :: b

    biomass = 0
    do i = 1, t%n_branches
        b => t%branches(i)%p
        biomass = biomass + b%length*0.25*pi*b%diameter**2.0
    end do
    B1 = B1 + biomass
    B2 = B2 + biomass**2.0
end subroutine biomass_calculation

!*****************************************************************************80  
subroutine leaves_extract(trees,trees_logical,Leaves)
    implicit none    
    type(tree_pointer), dimension(:), intent(in)    :: trees
    logical, dimension(:), intent(in)               :: trees_logical
    type(leaf), dimension(:), intent(inout)         :: Leaves
!*****************************************************************************80  
!
!! LEAVES_EXTRACT build the array Leaves to store the location and light on each
!       leaf of every tree
!    
    integer                     :: i, j, k, n
    type (branch), pointer      :: b
    real(kind=working_precis)   :: X, Y, Z

    i = 1
    k = 0
    do while (trees_logical(i))
        n = trees(i)%p%n_leaves
        do j = 1, n
            k = k + 1
            b => trees(i)%p%leaves(j)%p
            X = b%location(1) + b%length * b%unit_t(1)
            Y = b%location(2) + b%length * b%unit_t(2)
            Z = b%location(3) + b%length * b%unit_t(3)
            Leaves(k)%location = (/ X, Y, Z /)
            Leaves(k)%light = 0
        end do
        i = i + 1
    end do
end subroutine leaves_extract

!*****************************************************************************80  
subroutine light_interception(Leaves,Nleaves,SizeLeaf,elev,azim,k_light)
    implicit none    
    type(leaf), dimension(:), intent(inout)         :: Leaves
    integer, intent(in)                             :: Nleaves, k_light
    real(kind=working_precis), intent(in)           :: SizeLeaf
    real(kind=working_precis), intent(in)           :: elev, azim
!*****************************************************************************80  
!
!! LIGHT_INTERCEPTION calculates the light reaching the leaves
!    
    integer                                         :: i, j, xmin, xmax, ymin, ymax, x, y
    integer, dimension(Nleaves)                     :: ind
    integer, dimension(:,:), allocatable            :: shadow
    real                                            :: X0, Y0, Z0, Xp
    real, dimension(Nleaves)                        :: Xposition, Yposition, Zposition
    
!     print *, 'light 1'
    do i = 1, Nleaves
        X0 = Leaves(i)%location(1)
        Y0 = Leaves(i)%location(2)
        Z0 = Leaves(i)%location(3)
        Xp = X0 * cos(azim) + Y0 * sin(azim)
        Xposition(i) = Xp * cos(elev) + Z0 * sin(elev)
        Yposition(i) =-X0 * sin(azim) + Y0 * cos(azim)
        Zposition(i) =-Xp * sin(elev) + Z0 * cos(elev)
    end do
    
!     print *, 'light 2'
    call qsortr(Zposition, ind, Nleaves)
    
!     print *, 'light 3'
    xmin = nint(minval(Xposition/SizeLeaf))
    xmax = nint(maxval(Xposition/SizeLeaf))
    ymin = nint(minval(Yposition/SizeLeaf))
    ymax = nint(maxval(Yposition/SizeLeaf))
    
!     print *, 'allocate', k_light
    allocate(shadow(xmax-xmin+1,ymax-ymin+1))
    shadow = 0
!     print *, 'size:', xmax-xmin+1, ymax-ymin+1
    
    do i = 1, Nleaves
        j = ind(Nleaves - i + 1)
        x = nint(Xposition(j)/SizeLeaf) - xmin + 1
        y = nint(Yposition(j)/SizeLeaf) - ymin + 1
        Leaves(j)%light(k_light) =  max(1 - shadow(x,y), 0) 
        shadow(x,y) = shadow(x,y) + 1
    end do

    deallocate(shadow)    
end subroutine light_interception


!*****************************************************************************80  
subroutine light_on_trees(Leaves,trees,trees_logical)
    implicit none    
    type(leaf), dimension(:), intent(in)               :: Leaves
    type(tree_pointer), dimension(:), intent(inout)    :: trees
    logical, dimension(:), intent(in)                  :: trees_logical
!*****************************************************************************80  
!
!! COUNT_LEAVES count the total number of leaves (on all trees)
!    
    integer                     :: i, j, k, n
    type (branch), pointer      :: b

    i = 1
    k = 0
    do while (trees_logical(i))
        n = trees(i)%p%n_leaves
        do j = 1, n
            k = k + 1
            trees(i)%p%leaves(j)%p%light = 1.0 * sum( Leaves(k)%light(:) ) / (Nazim*Nelev)
!             print *, 'i', i, 'j', j, 'LI', sum( Leaves(k)%light(:) ), 'LR', trees(i)%p%leaves(j)%p%light
        end do
        i = i + 1
    end do
end subroutine light_on_trees

!*****************************************************************************80  
subroutine order_tree (t,err)
    implicit none    
    type(tree), pointer, intent(in)         :: t
    logical, intent(out)                    :: err
!*****************************************************************************80  
!
!! ORDER_TREE returns vectors of pointers BRANCHES(:) and LEAVES (:) 
!      with branches ordered from trunk to leaves
!    
    integer                     :: i, j, size_current, n_b
    logical                     :: next_step
    type (branch), pointer      :: b
    type (branch_pointer), dimension(:), allocatable :: currents, nexts

    allocate(currents(t%n_branches), nexts(t%n_branches))
    
    currents(1)%p => t%trunk
    size_current = 1
    err = .false.

    t%n_leaves = 0
    do                                  ! Find the leaves in the architecture
        j = 0                           ! starting from the trunk and going up
        do i = 1, size_current
            if ( associated(currents(i)%p%left)  ) then
                j = j+1
                nexts(j)%p => currents(i)%p%left 
            end if
            if ( associated(currents(i)%p%right) ) then
                j = j+1
                nexts(j)%p => currents(i)%p%right 
            end if                
        end do

        if ( j == 0 ) then 
            if ( associated(currents(1)%p, t%trunk) ) then
                t%n_leaves = 1
                t%leaves(1) = currents(1)
            end if
            exit
        end if
                
        size_current = j
        do i = 1, size_current
            currents(i)%p => nexts(i)%p
            if (.not. associated(currents(i)%p%left ) .and. &
                .not. associated(currents(i)%p%right)) then
                t%n_leaves = t%n_leaves + 1
                t%leaves(t%n_leaves) = currents(i)
            end if
        end do    
    end do

    do i = 1, t%n_leaves                  ! Start with the leaves
        t%leaves(i)%p%nb_leaves = 1
        currents(i)%p => t%leaves(i)%p
        t%branches(i)%p => t%leaves(i)%p    
        currents(i)%p%marker = .true.
        ! print *, 'leaf nb:', currents(i)%p%generation
    end do
    size_current = t%n_leaves
    n_b          = t%n_leaves
        
    do                                  ! Open loop, it exits when it reach
        j = 0                           ! the trunk
        do i = 1, size_current          ! This DO loop finds the branches
            b => currents(i)%p%parent   ! below the current ones and puts
            next_step = .true.          ! them in the NEXTS pointer array
            ! print *, ' parent:', b%generation
            if ( associated(b) ) then
                if ( associated(b%left) ) then
                    if ( .not. b%left%marker )  next_step = .false.
                end if
                if ( associated(b%right) ) then
                    if ( .not. b%right%marker )  next_step = .false.
                end if
                if ( b%marker ) next_step = .false.
            else
                next_step = .false.
            end if
            if ( next_step ) then
                j = j + 1
                nexts(j)%p => b
                b%marker = .true.
            end if
        end do
            
        if ( j == 0 ) then              ! If there is nothing in the NEXTS
            exit                        ! array, it exits the DO loop
        end if

        size_current = j
        currents(1:size_current) = nexts(1:size_current)

        do i = 1, size_current
            b => currents(i)%p
            b%nb_leaves = 0
            if ( associated(b%left) ) then
                b%nb_leaves = b%nb_leaves + b%left%nb_leaves     
            end if
            if ( associated(b%right) ) then
                b%nb_leaves = b%nb_leaves + b%right%nb_leaves      
            end if
        end do

        t%branches(n_b+1:n_b+size_current) = nexts(1:size_current)
        n_b = n_b+size_current
    end do
    
    currents(1:n_b) = t%branches(1:n_b)
    t%branches(1:n_b) = currents(n_b:1:-1)
    
    if (n_b .ne. t%n_branches) then
        err = .true.
        print *, 'order_tree: n_b .ne. n_branches'
    end if
    if (.not. associated(t%branches(1)%p,t%trunk)) err = .true.
    deallocate(currents, nexts)
    
    do i=1, n_b
        t%branches(i)%p%marker = .false.
    end do
    
end subroutine order_tree

!*****************************************************************************80  
subroutine primary_growth(t,TwigLength,TwigDiameter,generation,Nseeds)
    implicit none
    type(tree), intent(inout)               :: t 
    real(kind=working_precis), intent(in)   :: TwigLength, TwigDiameter
    integer, intent(in)                     :: generation
    integer, intent(out)                    :: Nseeds
!*****************************************************************************80  
!
!! PRIMARY_GROWTH creates new branches from available ressources
!    
    integer                   :: i, j, n, Nnewleaves
    integer                   :: k, perm(t%n_leaves)
    integer                   :: ind(t%n_leaves)
    real(kind=working_precis) :: pi = acos(-1.0)
    real(kind=working_precis) :: VolumeTwig
    real(kind=working_precis) :: light(t%n_leaves)
    real(kind=working_precis) :: Pnewseeds, Pnewleaves, Phototropism
    real(kind=working_precis) :: t1, t2, g1, g2
    logical                   :: Primary, Ground
    type (branch), pointer    :: b

    VolumeTwig = 0.25 * pi * TwigLength * TwigDiameter**2
    call neural_reserve(t%n_leaves,t%Reserve/t%n_leaves/VolumeTwig,t%NNreserve, &
                        Pnewseeds,Pnewleaves,Phototropism)
    Nseeds = floor(Pnewseeds*t%Reserve / (5.0*VolumeTwig))    
    Nnewleaves = floor(Pnewleaves*t%Reserve / (2.0*VolumeTwig))   

    do i = 1, t%n_leaves
        light(i) = t%leaves(i)%p%light
    end do 
    call qsortd(light, ind, t%n_leaves)
    n = min(t%n_leaves, Nnewleaves)
    k = n + nint( (1.0 - Phototropism) * (t%n_leaves - n) )

    call rperm2(k,perm(1:k))
!     print *, 'n',n,'k', k, 'perm', perm(1:k), 'ind', ind

    do i = 1, n
        j = ind(t%n_leaves - perm(i) + 1)
!         print *, 'j', j
        b => t%leaves(j)%p
        t%Reserve = t%Reserve - 2.0 * VolumeTwig
        Ground = ( b%location(3) + b%length * b%unit_t(3) > 0.99 * b%length )
        if ( Ground ) then
            t1 = t%theta1 + 0.174532925 * random_normal()
            t2 = t%theta2 + 0.174532925 * random_normal()
            g1 = t%gamma1 + 0.174532925 * random_normal()
            g2 = t%gamma2 + 0.174532925 * random_normal() 
!             print *, 't1, t2, g1, g2' , t1, t2, g1, g2          
            call new_branches(b, t%branches, t%n_branches, t1, t2, &
&                             g1, g2, TwigLength,generation)
        end if        
    end do  
!     print *, 'Reserve', t%Reserve/VolumeTwig, 'n_branches', t%n_branches, 'n', n
end subroutine primary_growth

!*****************************************************************************80  
subroutine new_branches(mother,branches,n_branches,theta1,theta2, &
&                        gamma1,gamma2,length,generation)
    implicit none    
    type (branch), pointer, intent(in)                 :: mother
    integer, intent(inout)                             :: n_branches
    type (branch_pointer), dimension(:), intent(inout) :: branches
    real(kind=working_precis), intent(in)              :: theta1, theta2
    real(kind=working_precis), intent(in)              :: gamma1, gamma2
    real(kind=working_precis), intent(in)              :: length
    integer, intent(in)                                :: generation
!*****************************************************************************80  
!
!! NEW_BRANCHES creates two new daughter branches from a mother branch and
!   update the arrays LEAVES and BRANCHES
!
    real(kind=working_precis), dimension(1:3)    :: location, t1, b1, t2, b2
    type (branch), pointer  :: left => null(), right => null()
    
    location = mother%location + mother%length * mother%unit_t
    call daughter_unit_vectors(t1,b1,mother,theta1,gamma1)
    call daughter_unit_vectors(t2,b2,mother,theta2,gamma2)

    allocate(left)
    left%parent => mother
    left%generation = generation
    left%length = length
    left%unit_t = t1
    left%unit_b = b1
    left%location = location

    allocate(right)
    right%parent => mother
    right%generation = generation
    right%length = length
    right%unit_t = t2
    right%unit_b = b2
    right%location = location

    mother%left  => left
    mother%right => right
    mother%light = 0.0

    branches(n_branches+1)%p => left
    branches(n_branches+2)%p => right     
    n_branches = n_branches + 2
end subroutine new_branches

!*****************************************************************************80  
subroutine daughter_unit_vectors (t,n,mother,theta,gamma)
    implicit none    
    type (branch), pointer, intent(in)                      :: mother
    real(kind=working_precis), intent(in)                   :: theta, gamma
    real(kind=working_precis), dimension(1:3), intent(out)  :: t, n
!*****************************************************************************80  
!
!! DAUGHTER_UNIT_VECTORS computes the unit vectors of daughter branches from
!   the mother branch
!
    real(kind=working_precis), dimension(1:3) :: b1, b2
    
    call cross_product(mother%unit_t, mother%unit_b, b1) 
    t = cos(theta) * mother%unit_t + sin(theta) * b1
    call cross_product(t, mother%unit_b, b2) 
    n = cos(gamma) * mother%unit_b + sin(gamma) * b2

end subroutine daughter_unit_vectors

!*****************************************************************************80  
subroutine cut_branch (start,n_branches,lefttrue,branches_cut)   
    implicit none    
    type (branch), pointer, intent(in)  :: start
    integer, intent(in)                 :: n_branches
    logical, intent(in)                 :: lefttrue 
    integer, intent(out)                :: branches_cut
!*****************************************************************************80  
!
!! CUT_BRANCH does:
!  - deallocate all branches after START
!  - nullify START%left or START%right (depending on LEFTTRUE value)
!  - update the value of N_BRANCHES (but not the vector of pointers BRANCHES(:))
!
    integer :: i, j, size_current
    type (branch_pointer), dimension(:), allocatable :: currents, nexts
    
    allocate(currents(n_branches), nexts(n_branches))
    branches_cut = 0
    
    if (lefttrue .and. associated(start%left)) then
        currents(1)%p => start%left
        start%left => null()
    elseif (.not. lefttrue .and. associated(start%right)) then
        currents(1)%p => start%right
        start%right => null()
    else
        write (*,*) ' no branch to cut'
        return
    end if
    
    size_current = 1    
    branches_cut = 0
    do 
        j = 0
        do i = 1, size_current
            if ( associated  (currents(i)%p%left)  ) then
                j = j+1
                nexts(j)%p => currents(i)%p%left 
            end if
            if ( associated  (currents(i)%p%right) ) then
                j = j+1
                nexts(j)%p => currents(i)%p%right 
            end if
            ! print *, 'branch cut geneneration:', currents(i)%p%generation
            ! call print_branch(currents(i)%p)
            if(associated(currents(i)%p)) deallocate(currents(i)%p)
            branches_cut = branches_cut + 1                
        end do
        
        if ( j == 0 ) then 
            exit
        end if
                
        size_current = j
        do i = 1, size_current
            currents(i)%p => nexts(i)%p
        end do    
    end do

    deallocate(currents,nexts)
end subroutine cut_branch

!*****************************************************************************80  
subroutine wind_force (br,V,force,moment)
    implicit none
    real(kind=working_precis), dimension(3), intent(in)     :: V
    real(kind=working_precis), dimension(3), intent(out)    :: force, moment
    type (branch), pointer, intent(in)                      :: br
!*****************************************************************************80  
!
!! WIND_FORCE calculates the force and moment exerted by the wind of velocity
!      V on a given branch
!    
    real(kind=working_precis)                :: L, d, costheta
    real(kind=working_precis), dimension(3)  :: x, u, t, n, b, Nn
    
    x = br%location
    L = br%length
    t = br%unit_t
    d = br%Diameter
    
    u = V/norm2(V)
    call cross_product(t,u,Nn)
    costheta = norm2(Nn)
    n = Nn/costheta
    call cross_product(n,t,b)
    
    force = norm2(V)**2.0 * d * L * costheta**2.0 * b
    call cross_product(0.5*L*t,force,moment) 
end subroutine wind_force

!*****************************************************************************80  
subroutine calculate_stresses (t,S0,Cy)
    implicit none    
    type(tree), pointer, intent(in)         :: t
    real(kind=working_precis), intent(in)   :: S0, Cy
!*****************************************************************************80  
!
!! CALCULATE_STRESSES calculates the maximum stresses in each branch when
!      the wind is incremented of 45 degrees
!    
    integer                                 :: i, angle
    real(kind=working_precis)               :: pi = acos(-1.0)
    real(kind=working_precis), dimension(3) :: U = (/ 0, 0, 0 /), force, moment 
    real(kind=working_precis), dimension(3) :: bend_moment, torqueU
    real(kind=working_precis), dimension(3) :: T1 = (/ 0, 0, 0 /), T2 = (/ 0, 0, 0 /)
    real(kind=working_precis), dimension(3) :: M1 = (/ 0, 0, 0 /), M2 = (/ 0, 0, 0 /)
    type (branch), pointer                  :: b
    
    do i = 1, t%n_branches
        t%branches(i)%p%max_stress = 0
    end do
    
    do angle = 1, 4                           ! There are 4 angles considered
        U(1:2) = (/ cos(pi*angle/4), sin(pi*angle/4) /)
        do i = 1, t%n_leaves                  ! Compute the force on each final
            b => t%leaves(i)%p                ! branch (leaf+branch contribute)
            call wind_force(b,U,force,moment)
            b%T = force + S0*U*norm2(U)
            call cross_product(b%length*b%unit_t,S0*U,torqueU)
            b%M = moment + torqueU
            call cross_product(b%unit_t,b%M,bend_moment)
            b%stress = 16.0/pi * Cy * norm2(bend_moment) / (b%diameter**3.0)
            b%max_stress = max(b%stress,b%max_stress)
!             print *, 'L Stress', b%stress, b%diameter, Cy, norm2(bend_moment)
        end do
        
        do  i = 1, t%n_branches - t%n_leaves    ! Compute the force on each branch
            b => t%branches(t%n_branches - t%n_leaves - i + 1)%p
            if ( associated(b%left) ) then
                T1 = b%left%T
                M1 = b%left%M
            end if
            if ( associated(b%right) ) then
                T2 = b%right%T
                M2 = b%right%M
            end if
            call wind_force(b,U,force,moment)
            b%T = force + T1 + T2
            call cross_product(b%length * b%unit_t, T1 + T2, torqueU)
            b%M = moment + M1 + M2 + torqueU
            call cross_product(b%unit_t, b%M, bend_moment)
            b%stress = 16.0/pi * Cy * norm2(bend_moment) / (b%diameter**3.0)
            b%max_stress = max(b%stress,b%max_stress)
!             print *, 'B Stress', b%stress, b%diameter, Cy, norm2(bend_moment)
            T1 = (/ 0, 0, 0 /)
            T2 = (/ 0, 0, 0 /)
            M1 = (/ 0, 0, 0 /)
            M2 = (/ 0, 0, 0 /)
        end do
    end do
end subroutine calculate_stresses

!*****************************************************************************80  
subroutine requested_growth (t,MaintenanceH)
    implicit none    
    type(tree), pointer, intent(in)         :: t
    real(kind=working_precis), intent(in)   :: MaintenanceH  
!*****************************************************************************80
!
!! REQUESTED_GROWTH calculates the growth volume asked by each branch based
!      on the value of MAX_STRESS
!    
    integer                   :: i
    real(kind=working_precis) :: pi = acos(-1.0)
    real(kind=working_precis) :: vol_actual, vol_wished, Safety
    type (branch), pointer    :: b

    do i = 1, t%n_branches
        b => t%branches(i)%p
        vol_actual = 0.25 * pi * b%diameter**2 * b%length
        b%maintenance_vol = pi * b%diameter * MaintenanceH
        call neural_branch(b%nb_leaves,b%max_stress,t%NNbranch,Safety)
        vol_wished = Safety * vol_actual * b%max_stress**(2.0/3.0)
        b%vol_growth = max(0.,vol_wished-vol_actual) + b%maintenance_vol
        if ( associated (b%parent) ) then
            b%vol_summed = b%parent%vol_summed*b%nb_leaves/b%parent%nb_leaves+ &
&                          b%vol_growth
        else 
            b%vol_summed = b%vol_growth
        end if
    end do
end subroutine requested_growth

!*****************************************************************************80  
subroutine neural_branch(nb_leaves,max_stress,NNbranch,Safety)
    implicit none    
    integer, intent(in)                                 :: nb_leaves
    real(kind=working_precis), intent(in)               :: max_stress
    real(kind=working_precis), dimension(:), intent(in) :: NNbranch
    real(kind=working_precis), intent(out)              :: Safety
!*****************************************************************************80
!
!! NEURAL_BRANCH calculates the volume requested by a branch given certain inputs
!    
    integer, parameter          :: Nhidden=3, Nin=2, Nout=1
    integer                     :: i
    real(kind=working_precis)   :: pi=acos(-1.0)
    real(kind=working_precis)   :: M1(Nhidden,Nin),M2(Nout,Nhidden+1)
    real(kind=working_precis)   :: X(Nin), Z(Nhidden), Zp(Nhidden+1), F(Nout)
    real(kind=working_precis)   :: toto(Nhidden*(Nin+Nout)+Nout)

    do i = 1, Nhidden*(Nin+Nout) + Nout
        toto(i) = tan( (NNbranch(i) - 0.5) * pi * 0.99 )
    end do
    M1 = reshape( toto(            1:Nhidden*Nin       ), (/ Nhidden, Nin  /)) 
    M2 = reshape( toto(Nhidden*Nin+1:Nhidden*(Nin+Nout)+Nout), (/ Nout, Nhidden + 1/))
    M1(1,1) = 0.0
    M1(Nhidden,Nin) = 0.0
    X(1) = 0.01*nb_leaves
    X(2) = max_stress
    Z = matmul(M1,X)
    do i = 1,Nhidden
        Zp(i) = tanh( 5.0 * Z(i) ) / Nhidden
    end do
    Zp(Nhidden + 1) = 1.0 / Nhidden
    F = matmul(M2,Zp)
    Safety   = max(0.0,F(1)+1.0)
end subroutine neural_branch

!*****************************************************************************80  
subroutine neural_reserve(nb_leaves,vol_relative,NNreserve,Pnewseeds,Pnewleaves,Phototropism)
    implicit none    
    integer, intent(in)                                 :: nb_leaves
    real(kind=working_precis), intent(in)               :: vol_relative
    real(kind=working_precis), dimension(:), intent(in) :: NNreserve
    real(kind=working_precis), intent(out)              :: Pnewseeds
    real(kind=working_precis), intent(out)              :: Pnewleaves
    real(kind=working_precis), intent(out)              :: Phototropism
!*****************************************************************************80
!
!! NEURAL_BRANCH calculates the volume requested by a branch given certain inputs
!    
    integer, parameter          :: Nhidden=3, Nin=2, Nout=3
    integer                     :: i
    real(kind=working_precis)   :: pi=acos(-1.0)
    real(kind=working_precis)   :: sum_p
    real(kind=working_precis)   :: M1(Nhidden,Nin),M2(Nout,Nhidden+1)
    real(kind=working_precis)   :: X(Nin), Z(Nhidden), Zp(Nhidden+1), F(Nout)
    real(kind=working_precis)   :: toto(Nhidden*(Nin+Nout)+Nout)

    do i = 1, Nhidden*(Nin+Nout) + Nout
        toto(i) = tan( (NNreserve(i) - 0.5) * pi * 0.99 )
    end do
    M1 = reshape( toto(            1:Nhidden*Nin       ), (/ Nhidden, Nin  /)) 
    M2 = reshape( toto(Nhidden*Nin+1:Nhidden*(Nin+Nout)+Nout), (/ Nout, Nhidden + 1/))
    M1(1,1) = 0.0
    M1(Nhidden,Nin) = 0.0
    X(1) = 0.01*nb_leaves
    X(2) = vol_relative
    Z = matmul(M1,X)
    do i = 1,Nhidden
        Zp(i) = tanh( 5.0 * Z(i) ) / Nhidden
    end do
    Zp(Nhidden + 1) = 1.0 / Nhidden
    F = matmul(M2,Zp)
    Pnewleaves   = min( max(0.0,F(1)+2.0), 4.0 ) / 4.0
    Pnewseeds    = min( max(0.0,F(2)+2.0), 4.0 ) / 4.0
    Phototropism = min( max(0.0,F(3)+2.0), 4.0 ) / 4.0
    sum_p = Pnewseeds + Pnewleaves
    if (sum_p > 1) then
        Pnewleaves = Pnewleaves / sum_p
        Pnewseeds  = Pnewseeds  / sum_p
    end if 
end subroutine neural_reserve

!*****************************************************************************80  
subroutine secondary_growth (t,VolumePerLeaf)
    implicit none    
    type(tree), pointer, intent(in)         :: t
    real(kind=working_precis), intent(in)   :: VolumePerLeaf
!*****************************************************************************80  
!
!! SECONDARY_GROWTH implement the radial growth of branches
!   
    integer                   :: i
    real(kind=working_precis) :: pi = acos(-1.0)
    real(kind=working_precis) :: fraction, vol_growth_branches, growth
    real(kind=working_precis) :: Total_reserv
    type (branch), pointer    :: b

    ! allocation from leaves
    do i = 1, t%n_leaves
        b => t%leaves(i)%p
        vol_growth_branches = min(b%light*VolumePerLeaf, b%vol_summed)
        t%Reserve = t%Reserve + (b%light*VolumePerLeaf - vol_growth_branches)
!         if (t%Reserve>10) print *, '***', t%Reserve, b%light, VolumePerLeaf, vol_growth_branches,b%vol_summed
        fraction = b%vol_growth/b%nb_leaves/b%vol_summed
        growth = vol_growth_branches * fraction - b%maintenance_vol/b%nb_leaves
        b%diameter = sqrt( max(1E-3, growth/b%length*4/pi + b%diameter**2) )
        do while ( associated(b%parent) )
            b => b%parent
            fraction = b%vol_growth/b%nb_leaves/t%leaves(i)%p%vol_summed
            growth = vol_growth_branches * fraction - b%maintenance_vol/b%nb_leaves
            b%diameter = sqrt( max(0., growth/b%length*4/pi + b%diameter**2) )
        end do
    end do
end subroutine secondary_growth

!*****************************************************************************80  
subroutine pruning (t,S0,Cy,U)
    implicit none    
    type(tree), pointer, intent(in)                     :: t
    real(kind=working_precis), intent(in)               :: S0, Cy
    real(kind=working_precis), dimension(3), intent(in) :: U
!*****************************************************************************80  
!
!! PRUNING calculates the stresses in each branch and determines which branches
!   are pruned
!    
    integer                                 :: i
    integer                                 :: branches_cut, ncut
    real(kind=working_precis)               :: pi = acos(-1.0), m = 10.0
    real(kind=working_precis)               :: rndm(t%n_branches)
    real(kind=working_precis)               :: vol_relative
    real(kind=working_precis), dimension(3) :: force, moment, bend_moment,torqueU
    real(kind=working_precis), dimension(3) :: T1 = (/ 0, 0, 0 /), T2 = (/ 0, 0, 0 /)
    real(kind=working_precis), dimension(3) :: M1 = (/ 0, 0, 0 /), M2 = (/ 0, 0, 0 /)
    logical                                 :: cut_logical(t%n_branches)
    type (branch), pointer                  :: b
        
    cut_logical = .false.
    call random_number(rndm)

    do i = 1, t%n_leaves                  ! Compute the force on each final
        b => t%leaves(i)%p                ! branch (leaf+branch contribute)
        call wind_force(b,U,force,moment)
        b%T = force + S0*U*norm2(U)
        call cross_product(b%length*b%unit_t,S0*U,torqueU)
        b%M = moment + torqueU
        call cross_product(b%unit_t,b%M,bend_moment)
        b%stress = 16.0/pi * Cy * norm2(bend_moment) / (b%diameter**3.0)
        vol_relative = (b%diameter/0.1)**2.0
        if ( 1.0 - exp(-vol_relative*b%stress**m) > rndm(i) ) then
            cut_logical(t%n_branches - i + 1) = .true.
!             print *, 'CUT TWIG, diameter', &
!               t%branches(t%n_branches - i + 1)%p%diameter, b%diameter
        end if  
    end do

!     print *, 'branches to cut'
    do  i = 1, t%n_branches - t%n_leaves - 1    ! Compute the force on each branch
        b => t%branches(t%n_branches - t%n_leaves - i + 1)%p
        if ( associated(b%left) ) then
            T1 = b%left%T
            M1 = b%left%M
        end if
        if ( associated(b%right) ) then
            T2 = b%right%T
            M2 = b%right%M
        end if
        call wind_force(b,U,force,moment)
        b%T = force + T1 + T2
        call cross_product(b%length * b%unit_t, T1 + T2, torqueU)
        b%M = moment + M1 + M2 + torqueU
        call cross_product(b%unit_t, b%M, bend_moment)
        b%stress = 16.0 / pi * Cy * norm2(bend_moment) / (b%diameter**3.0)
        T1 = (/ 0, 0, 0 /)
        T2 = (/ 0, 0, 0 /)
        M1 = (/ 0, 0, 0 /)
        M2 = (/ 0, 0, 0 /)
        vol_relative = (b%diameter/0.1)**2.0
        if ( 1.0 - exp(-vol_relative*b%stress**m) > rndm(i+t%n_leaves) ) then
            cut_logical(t%n_branches - t%n_leaves - i + 1) = .true.
        end if        
    end do

!     print *,'cutting branches'
    ncut = 0    
    do  i = 1, t%n_branches - 1    ! Compute the force on each branch
        if ( cut_logical(t%n_branches - i + 1) ) then
            if ( associated(t%branches(t%n_branches - i + 1)%p) ) then
                b => t%branches(t%n_branches - i + 1)%p
                if ( associated(b%parent%left,b) ) then
                    call cut_branch(b%parent,t%n_branches,.true.,branches_cut)
                    ncut = ncut + branches_cut
                elseif ( associated(b%parent%right,b) ) then
                    call cut_branch(b%parent,t%n_branches,.false.,branches_cut)
                    ncut = ncut + branches_cut
                end if
            end if
        end if
    end do
    t%n_branches = t%n_branches - ncut    
!     if (ncut>0) print *,'done cutting branches', ncut
end subroutine pruning

!*****************************************************************************80  
subroutine save_tree (branches,n_branches,save_file,iout,reserve)
    implicit none    
    integer, intent(in)                             :: n_branches, iout
    type (branch_pointer), dimension(:), intent(in) :: branches
    character(len=30), intent(in)                   :: save_file
    real(kind=working_precis), intent(in)           :: reserve    
!*****************************************************************************80  
!
!! SAVE_TREE save the data of the current tree in the file:
!            SAVE_FILE // I // '.dat'
!
    integer      :: i
    character(8) :: number_file
    
    write(number_file,fmt='(I8.8)') iout
    open(unit=13,file=trim(save_file) // number_file // '.dat', &
&     status="replace", position="rewind")

    do i = 1, n_branches
        write(13,*) &
&                   branches(i)%p%generation, &
&                   branches(i)%p%diameter,   &
&                   branches(i)%p%length,     &
&                   branches(i)%p%unit_t,     &
&                   branches(i)%p%unit_b,     &
&                   branches(i)%p%location,   &
&                   branches(i)%p%light,      &
&                   reserve
    end do
    close(13)
end subroutine save_tree

!*****************************************************************************80  
subroutine make_statistics(n_leaves,branches,n_branches)
    implicit none    
    integer, intent(in)                             :: n_branches, n_leaves
    type (branch_pointer), dimension(:), intent(in) :: branches
!*****************************************************************************80  
!
!! MAKE_STATISTICS computes Strahler and botanic ranks of the tree branches
!
    integer                     :: i, j, k
    type (branch), pointer      :: b

    do i = 1,n_branches
        branches(i)%p%Strahler = 0
        branches(i)%p%botanic = 0
    end do

    do i = 1,n_leaves
        branches(n_branches-i+1)%p%Strahler = 1
    end do 
    branches(1)%p%botanic = 1

    do i = 1,n_branches
        b => branches(i)%p 
        if (associated(b%left) .and. .not.associated(b%right) .and. b%botanic>0) &
&               b%left%botanic =  b%botanic            
        if (associated(b%right) .and. .not.associated(b%left) .and. b%botanic>0) &
&               b%right%botanic =  b%botanic  
        if (associated(b%right) .and. associated(b%left) .and. b%botanic>0) then
            if (b%left%diameter**2.0 .ge. 1.0*b%right%diameter**2.0) then
                b%left%botanic  = b%botanic
                b%right%botanic = b%botanic + 1
            else if (b%right%diameter**2.0 > 1.0*b%left%diameter**2.0) then
                b%left%botanic  = b%botanic + 1
                b%right%botanic = b%botanic
            else
                b%left%botanic  = b%botanic + 1
                b%right%botanic = b%botanic + 1
            end if
        end if
    end do   

    do while (branches(1)%p%Strahler == 0)
        do i = n_branches-n_leaves,1,-1
            b => branches(i)%p 
            if (associated(b%left) .and. .not.associated(b%right) .and. b%left%Strahler>0) &
&               b%Strahler =  b%left%Strahler            
            if (associated(b%right) .and. .not.associated(b%left) .and. b%right%Strahler>0) &
&               b%Strahler =  b%right%Strahler  
            if (associated(b%right) .and. associated(b%left)  &
&               .and. b%left%Strahler>0 .and. b%right%Strahler>0) then
                if (b%left%Strahler == b%right%Strahler) then
                    b%Strahler = b%left%Strahler + 1
                else 
                    b%Strahler = max(b%left%Strahler,b%right%Strahler)
                end if
            end if
        end do     
    end do
    
end subroutine make_statistics

!*****************************************************************************80  
subroutine save_statistics(n_leaves,branches,n_branches,save_file,iout,Nmax)
    implicit none    
    integer, intent(in)                             :: n_branches, n_leaves
    integer, intent(in)                             :: iout, Nmax
    type (branch_pointer), dimension(:), intent(in) :: branches
    character(len=30), intent(in)                   :: save_file
!*****************************************************************************80  
!
!! SAVE_STATISTICS save in a file statistics of the tree
!
    integer                     :: i, j, k
    integer                     :: NsegmentsS(Nmax), NsegmentsB(Nmax)
    integer                     :: Strahlermax, Sorder
    real(kind=working_precis)   :: Tokunaga(Nmax,Nmax)
    real(kind=working_precis)   :: length_Strahler(Nmax)
    real(kind=working_precis)   :: length_botanic(Nmax)
    character(8)                :: number_file
    type (branch), pointer      :: b

    NsegmentsS = 0
    NsegmentsB = 0
    length_Strahler = 0.0
    length_botanic = 0.0
    Tokunaga = 0
    NsegmentsS(1) = n_leaves
    NsegmentsB(1) = 1

    do i = 1,n_branches-n_leaves
        b => branches(i)%p 
        if ( associated(b%right) .and. associated(b%left) ) then
            NsegmentsB(b%botanic + 1) = NsegmentsB(b%botanic + 1) + 1
            if (b%left%botanic  == b%right%botanic ) then
                NsegmentsB(b%botanic + 1) = NsegmentsB(b%botanic + 1) + 1
            end if
            if (b%left%Strahler == b%right%Strahler) then
                Sorder = b%left%Strahler + 1
                NsegmentsS(Sorder) = NsegmentsS(Sorder) + 1
            end if
        end if
    end do  

    Strahlermax = branches(1)%p%Strahler
    if (iout==1) then
        open(unit=11,file='Z_Nsegments_' // trim(save_file) // '.dat', status="replace")
    else
        open(unit=11,file='Z_Nsegments_' // trim(save_file) // '.dat', position="append")
    end if
    write (11,*) iout, NsegmentsS
    close(11)

    if (iout==1) then
        open(unit=11,file='Z_Nsegments2_' // trim(save_file) // '.dat', status="replace")
    else
        open(unit=11,file='Z_Nsegments2_' // trim(save_file) // '.dat', position="append")
    end if
    write (11,*) iout, NsegmentsB
    close(11)

    do i = 1, Strahlermax
        j = 0 
        do k = 1,n_branches
            if (branches(k)%p%Strahler == i) j = j + 1
        end do
        length_Strahler(i) = 1.0 * j / NsegmentsS(i)
    end do

    do i = 1, Nmax
        j = 0 
        do k = 1,n_branches
            if (branches(k)%p%botanic == i) j = j + 1
        end do
        length_botanic(i) = 1.0 * j / NsegmentsB(i)
    end do

    if (iout==1) then
        open(unit=11,file='Z_length_' // trim(save_file) // '.dat', status="replace")
    else
        open(unit=11,file='Z_length_' // trim(save_file) // '.dat', position="append")
    end if
    write (11,*) iout, length_Strahler
    close(11)

    if (iout==1) then
        open(unit=11,file='Z_length2_' // trim(save_file) // '.dat', status="replace")
    else
        open(unit=11,file='Z_length2_' // trim(save_file) // '.dat', position="append")
    end if
    write (11,*) iout, length_botanic
    close(11)

    do i = 2, n_branches
        j = branches(i)%p%Strahler
        k = branches(i)%p%parent%Strahler
        if ( j < k ) then
            Tokunaga(k,j) = Tokunaga(k,j) + 1.0/NsegmentsS(k)
        end if
    end do
    do i = 2, Strahlermax
        Tokunaga(i,i-1) = Tokunaga(i,i-1) - 2.0
    end do
    if (iout==1) then
        open(unit=11,file='Z_Tokunaga_' // trim(save_file) // '.dat', status="replace")
    else
        open(unit=11,file='Z_Tokunaga_' // trim(save_file) // '.dat', position="append")
    end if
    write (11,*) iout, Tokunaga
    close(11)

    write(number_file,fmt='(I8.8)') iout
    open(unit=13,file='ZStat_' // trim(save_file) // number_file // '.dat', &
&     status="replace", position="rewind")

    do i = 1, n_branches
        write(13,*) &
&                   branches(i)%p%Strahler,   &
&                   branches(i)%p%diameter,   &
&                   branches(i)%p%length,     &
&                   branches(i)%p%unit_t,     &
&                   branches(i)%p%unit_b,     &
&                   branches(i)%p%location,   &
&                   branches(i)%p%light,      &
&                   branches(i)%p%botanic
    end do
    close(13)
end subroutine save_statistics

!*****************************************************************************80  
subroutine save_area(t,save_file,iout)
    implicit none    
    type(tree), intent(inout)               :: t 
    integer, intent(in)                     :: iout
    character(len=30), intent(in)           :: save_file
!*****************************************************************************80  
!
!! SAVE_AREA saves area-preserving statistics
!
    integer                     :: i
    character(8)                :: number_file
    type (branch), pointer      :: b
    real(kind=working_precis)   :: pi = acos(-1.0)

    do i = 1,t%n_branches
        t%branches(i)%p%distance_leaves = -1.0
    end do

    do i = 1,t%n_leaves
        b => t%leaves(i)%p
        b%distance_leaves = 0.5
        if ( associated(b%parent) ) then
            b%distance_ground = b%parent%distance_ground + 1
        else
            b%distance_ground = 1
!             print *, 'TRUNK', i
        end if
    end do

    do while (t%branches(1)%p%distance_leaves < 0)
        do i = t%n_branches - t%n_leaves, 1, -1
            b => t%branches(i)%p 
            if (associated(b%left) .and. .not.associated(b%right) .and. b%left%distance_leaves>0) &
&               b%distance_leaves =  b%left%distance_leaves + 1.0           
            if (associated(b%right) .and. .not.associated(b%left) .and. b%right%distance_leaves>0) &
&               b%distance_leaves =  b%right%distance_leaves + 1.0 
            if (associated(b%right) .and. associated(b%left)  &
&               .and. b%left%distance_leaves>0 .and. b%right%distance_leaves>0) then
                b%distance_leaves =  ( (b%left%distance_leaves  + 1.0)*b%left%nb_leaves + &
                                       (b%right%distance_leaves + 1.0)*b%right%nb_leaves ) / &
                                     ( b%left%nb_leaves + b%right%nb_leaves)
            end if
        end do     
    end do

    write(number_file,fmt='(I8.8)') iout
    open(unit=13,file='ZStat2_' // trim(save_file) // number_file // '.dat', &
&     status="replace", position="rewind")
    open(unit=11,file='ZAreaRatio_' // trim(save_file) // number_file // '.dat', &
&     status="replace", position="rewind")

    do i = 1, t%n_branches
        b => t%branches(i)%p
        write(13,*) &
&                   b%distance_leaves,      &
&                   b%diameter,             &
&                   b%length,               &
&                   b%unit_t,               &
&                   b%unit_b,               &
&                   b%location,             &
&                   b%distance_ground 

        if ( associated(b%right) .and. associated(b%left) ) then
            write(11,*) &
&                   b%distance_leaves,                                  &
&                   b%distance_ground,                                  &
&                   0.25 * pi *b%length       * b%diameter**2.0,        &
&                   0.25 * pi *b%left%length  * b%left%diameter**2.0,   &
&                   0.25 * pi *b%right%length * b%right%diameter**2.0
        end if 
    end do
    close(13)
    close(11)
end subroutine save_area

!*****************************************************************************80  
subroutine save_allocation(t,wind,Ntwigs,Nseeds,Npruned,save_file,iout)
    implicit none    
    type(tree), intent(inout)               :: t 
    real(kind=working_precis), intent(in)   :: wind
    integer, intent(in)                     :: Ntwigs,Nseeds,Npruned,iout
    character(len=30), intent(in)           :: save_file
!*****************************************************************************80  
!
!! SAVE_ALLOCATION saves time evolution of allocations
!
    integer                     :: i
    type (branch), pointer      :: b
    real(kind=working_precis)   :: pi = acos(-1.0)

    open(unit=16,file='ZAllocation.dat', status="unknown",     position="append") 
    write(16,*) iout, t%n_branches, t%n_leaves, wind, Ntwigs,Nseeds,Npruned, t%Reserve
    close(16)
end subroutine save_allocation

end module mod_tree

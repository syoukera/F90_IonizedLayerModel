module variables_module
    implicit none

    ! parameters for grid
    integer, parameter :: nr = 101
    integer, parameter :: nz = 101
    double precision, parameter :: dr        = 1.0d-4 ! distance between grid points [m]
    double precision, parameter :: dz        = 1.0d-4 ! distance between grid points [m]
    double precision, parameter :: length_r  = (nr-1)*dr ! length of calclation domain [m]
    double precision, parameter :: length_z  = (nz-1)*dz ! length of calclation domain [m]
    
    ! constants
    double precision, parameter :: epsilon_0    = 8.854187817d-12  ! vacuum permittivity [C/V m]
    double precision, parameter :: q_e          = 1.602176634d-19  ! elementary charge [C]
    double precision, parameter :: pi           = 3.141592653589d0
    double precision, parameter :: k_B          = 1.380649d-23 ! Boltzmann constant [J/K] 

    ! parameters for Gaussian profile
    double precision, parameter :: height_flame = length_z/2.0 ! height of flame [m]
    double precision, parameter :: a_thickness  = 1.0d-3 ! thickness parameter [m]

    ! parameters for transport and reactions
    double precision, parameter :: k_i = 1.0d20 ! rate coeficient of ionization ions/m3s
    double precision, parameter :: k_r = 2.4d-13 ! rate coeficient of recombination m3/ions s
    double precision, parameter :: K_pos = 2.9d-4 ! mobility of positive ions [m2/s V]
    double precision, parameter :: K_neg = 2.9d-4 ! mobility of negative ions [m2/s V]
    double precision, parameter :: K_ele = 0.4  ! mobility of electrons [m2/s V]
    double precision, parameter :: alpha = 0.0  ! ratio of electrons among the negatively charges species [0-1]
    double precision, parameter :: T = 298d0  ! temperature [K]

    ! diffusion coefficients is drived from Einstein Eq.
    double precision, parameter :: D_pos = K_pos*k_B*T/q_e ! diffusion coefficients of positive ions [m2/s]
    double precision, parameter :: D_neg = K_neg*k_B*T/q_e ! diffusion coefficients of negative ions [m2/s]
    double precision, parameter :: D_ele = K_ele*k_B*T/q_e ! diffusion coefficients of electrons [m2/s]
    
    ! parameters for boundary conditions
    double precision, parameter :: V_start      = 300D0 ! valtage for initial point [V]
    double precision, parameter :: V_end        = 0.0d0 ! voltage for end point [V]

    ! parameters for computation
    integer, parameter :: k_start = 1
    integer, parameter :: k_end   = 1000000
    integer, parameter :: k_step  = 100000
    double precision, parameter :: tolerance = 2.0d-6

    ! double precision, parameter :: omega_V   = 0.5d0 ! relaxation coefficient (1 < omega < 2)
    ! double precision, parameter :: omega_pos = 0.005d0 ! relaxation coefficient (1 < omega < 2)
    ! double precision, parameter :: omega_neg = 0.005d0 ! relaxation coefficient (1 < omega < 2)
    ! double precision, parameter :: omega_ele = 0.005d0 ! relaxation coefficient (1 < omega < 2)

    double precision, parameter :: omega_V   = 0.5d0 ! relaxation coefficient (1 < omega < 2)
    double precision, parameter :: omega_pos = 0.005d0 ! relaxation coefficient (1 < omega < 2)
    double precision, parameter :: omega_neg = 0.005d0 ! relaxation coefficient (1 < omega < 2)
    double precision, parameter :: omega_ele = 0.005d0 ! relaxation coefficient (1 < omega < 2)

    double precision :: error

    ! set variables arrays
    double precision :: distance_r(nr, nz) ! position of eac grid point [m]
    double precision :: distance_z(nr, nz) ! position of eac grid point [m]
    double precision :: V(nr, nz) ! electric potential [V]
    double precision :: E_r(nr, nz) ! electric field [V/m]
    double precision :: E_z(nr, nz) ! electric field [V/m]
    double precision :: n_pos(nr, nz) ! number density of positive ions [m-3]
    double precision :: n_neg(nr, nz) ! number density of negative ions [m-3]
    double precision :: n_ele(nr, nz) ! number density of electrons [m-3]
    double precision :: rho(nr, nz) ! density of electric charge [C/m3]
    
    double precision :: V_old(nr, nz)
    double precision :: n_pos_old (nr, nz)
    double precision, dimension(nr, nz) :: n_neg_old 
    double precision, dimension(nr, nz) :: n_ele_old 

    ! input data for flame height
    double precision :: normalized_flame_height(nr)
    double precision :: normalized_intensity(nr)

    ! output variables
    double precision :: current_density(nr, nz) ! current density [A/m3]
    double precision :: body_force(nr, nz) ! electric body force [N]

    contains

    subroutine initialize_variables()
        integer :: i, j
        
        do i = 1, nr
            do j = 1, nz

                distance_r(i, j) = (i-1) * dr 
                distance_z(i, j) = (j-1) * dz

                ! n_pos(i, j) = 0.0d0
                n_pos(i, j) = 1.0d13*max(exp(- pi*(distance_z(i, j) - height_flame)**2/a_thickness**2), 0.0)
                ! n_neg(i, j) = 0.0d0
                n_neg(i, j) = 1.0d13*max(exp(- pi*(distance_z(i, j) - height_flame)**2/a_thickness**2), 0.0)
                n_ele(i, j) = 0.0d0

                rho(i, j) = (n_pos(i, j) - n_neg(i, j) - n_ele(i, j))*q_e
            end do
        end do

        ! set boundary on V
        do i = 1, nr
            V(i, 1) = V_start
            V(i, nz) = V_end
        end do

        ! set initial conditions of V
        do i = 1, nr
            do j = 2, nz-1
                V(i, j) = 0.0d0
                ! V(i, j) = V_start + (V_end - V_start)*((j-1.0)/(nz-1.0))
            end do
        end do

        call update_electric_field()

    end subroutine initialize_variables

    subroutine load_flame_height()

        integer :: i
        character(len=100) :: filename
    
        filename = "input/flame_height_ch_10kV_posi_nr101"  ! 読み込むファイル名
    
        open(unit=10, file=filename, status="old", action="read")
    
        ! データを1行ずつ読み込み
        do i = 1, nr
            read(10, *) normalized_flame_height(i), normalized_intensity(i)
        end do
    
        close(10)
    
        ! ! 結果を確認
        ! do i = 1, nr
        !     print *, i, normalized_flame_height(i), normalized_intensity(i)
        ! end do

    end subroutine load_flame_height

    subroutine update_electric_field()
        integer :: i, j

        ! calclate electric field (E = -dV/dz)
        do i = 2, nr-1
            do j = 2, nz-1
                E_r(i, j) = -(V(i+1, j) - V(i-1, j)) / (2.0d0 * dr)
                E_z(i, j) = -(V(i, j+1) - V(i, j-1)) / (2.0d0 * dz)
            end do
        end do

        ! boundary conditions on z axis
        do i = 1, nr

            ! lower boundary
            E_r(i, 1) = 0.0d0
            E_z(i, 1) = -(V(i, 2) - V(i, 1)) / dz
            
            ! upper boundary
            E_r(i, nz) = 0.0d0
            E_z(i, nz) = -(V(i, nz) - V(i, nz-1)) / dz

        end do

        ! boundary conditions on r axis
        do i = 2, nz-1

            ! lower boundary
            E_r(1, j) = -(V(2, j) - V(1, j)) / dr
            E_z(1, j) = -(V(1, j+1) - V(1, j-1)) / (2.0d0 * dz)
            
            ! upper boundary
            E_r(nr, j) = -(V(nr, j) - V(nr-1, j)) / dr
            E_z(nr, j) = -(V(nr, j+1) - V(nr, j-1)) / (2.0d0 * dz)

        end do

    end subroutine update_electric_field

    subroutine export_variables(k)
        implicit none
        integer, intent(in) :: k
        integer :: i, j
        character(len=60) :: filename
        
        ! ! calculate current density
        ! do i = 2, nz-1

        !     current_density(i) = (D_pos*((n_pos(i+1)-n_pos(i-1))/(2.0*dz)) - K_pos*n_pos(i)*E(i))*(+q_e) &
        !                        + (D_pos*((n_neg(i+1)-n_pos(i-1))/(2.0*dz)) + K_neg*n_neg(i)*E(i))*(-q_e) &
        !                        + (D_ele*((n_ele(i+1)-n_pos(i-1))/(2.0*dz)) + K_ele*n_ele(i)*E(i))*(-q_e)

        ! end do

        
        ! create a unique filename using the integer i
        write(filename, '("potential_1d_", I0, ".dat")') k

        print *, "Output to file: ", filename

        ! output
        open(unit=1, file=filename, status='replace')
        do j = nz, 1, -1
            ! write(1, *) (V(i, j), i = 1, nr)
            ! write(1, *) (n_pos(i, j), i = 1, nr)
            write(1, *) (n_neg(i, j), i = 1, nr)
        end do
        close(1)

    end subroutine export_variables

    ! subroutine import_variables(k)
    !     implicit none
    !     integer, intent(in) :: k
    !     integer :: i
    !     character(len=60) :: filename
        
    !     ! ! calculate current density
    !     ! do i = 2, nz-1

    !     !     current_density(i) = (D_pos*((n_pos(i+1)-n_pos(i-1))/(2.0*dz)) - K_pos*n_pos(i)*E(i))*(+q_e) &
    !     !                        + (D_pos*((n_neg(i+1)-n_pos(i-1))/(2.0*dz)) + K_neg*n_neg(i)*E(i))*(-q_e) &
    !     !                        + (D_ele*((n_ele(i+1)-n_pos(i-1))/(2.0*dz)) + K_ele*n_ele(i)*E(i))*(-q_e)

    !     ! end do

        
    !     ! ! create a unique filename using the integer i
    !     ! write(filename, '("potential_1d_", I0, ".dat")') k
    !     filename = 'output/1kV_omega_V1.0_omega_ion0.05/potential_1d_100000.dat'
        
    !     print *, "Import from file: ", filename

    !     ! open file
    !     open(unit=1, file=filename, status='old')

    !     ! skip header
    !     read(1, '(A)')        
        
    !     do i = 1, nz
    !         read(1, '(7E24.16)') X(i), rho(i), V(i), E(i), n_pos(i), n_neg(i), n_ele(i)
    !     end do

    !     close(1)

    ! end subroutine import_variables

end module variables_module
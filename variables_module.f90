module variables_module
    implicit none

    ! parameters for grid
    integer, parameter :: nx = 101
    double precision, parameter :: dx        = 1.0d-4 ! distance between grid points [m]
    double precision, parameter :: length_x  = (nx-1)*dx ! length of calclation domain [m]
    
    ! constants
    double precision, parameter :: epsilon_0    = 8.854187817d-12  ! vacuum permittivity [C/V m]
    double precision, parameter :: q_e          = 1.602176634d-19  ! elementary charge [C]
    double precision, parameter :: pi           = 3.141592653589d0
    double precision, parameter :: k_B          = 1.380649d-23 ! Boltzmann constant [J/K] 

    ! parameters for Gaussian profile
    double precision, parameter :: height_flame = length_x/2.0 ! height of flame [m]
    double precision, parameter :: a_thickness  = 1.0d-3 ! thickness parameter [m]

    ! parameters for transport and reactions
    double precision, parameter :: k_i = 1.0d20 ! rate coeficient of ionization ions/m3s
    double precision, parameter :: k_r = 2.4d-13 ! rate coeficient of recombination m3/ions s
    double precision, parameter :: K_pos = 2.9d-4 ! mobility of positive ions [m2/s V]
    double precision, parameter :: K_neg = 2.9d-4 ! mobility of negative ions [m2/s V]
    double precision, parameter :: K_ele = 0.4d0  ! mobility of electrons [m2/s V]
    double precision, parameter :: alpha = 0.0  ! ratio of electrons among the negatively charges species [0-1]
    double precision, parameter :: T = 298d0  ! temperature [K]

    ! diffusion coefficients is drived from Einstein Eq.
    double precision, parameter :: D_pos = K_pos*k_B*T/q_e ! diffusion coefficients of positive ions [m2/s]
    double precision, parameter :: D_neg = K_neg*k_B*T/q_e ! diffusion coefficients of negative ions [m2/s]
    double precision, parameter :: D_ele = K_ele*k_B*T/q_e ! diffusion coefficients of electrons [m2/s]
    
    ! parameters for boundary conditions
    double precision, parameter :: V_start      = 1.0d3 ! valtage for initial point [V]
    double precision, parameter :: V_end        = 0.0d0 ! voltage for end point [V]

    ! parameters for computation
    integer, parameter :: k_start = 1
    integer, parameter :: k_end   = 100000
    integer, parameter :: k_step  = 10000
    double precision, parameter :: tolerance = 2.0d-6
    double precision, parameter :: omega_V   = 0.1d0 ! relaxation coefficient (1 < omega < 2)
    double precision, parameter :: omega_pos = 0.05d0 ! relaxation coefficient (1 < omega < 2)
    double precision, parameter :: omega_neg = 0.05d0 ! relaxation coefficient (1 < omega < 2)
    double precision, parameter :: omega_ele = 1.0d0 ! relaxation coefficient (1 < omega < 2)
    double precision :: error

    ! set variables arrays
    double precision :: X(nx) ! position of eac grid point [m]
    double precision :: V(nx) ! electric potential [V]
    double precision :: E(nx) ! electric field [V/m]
    double precision :: n_pos(nx) ! number density of positive ions [m-3]
    double precision :: n_neg(nx) ! number density of negative ions [m-3]
    double precision :: n_ele(nx) ! number density of electrons [m-3]
    double precision :: rho(nx) ! density of electric charge [C/m3]

    ! output variables
    double precision :: current_density(nx) ! current density [A/m3]
    double precision :: body_force(nx) ! electric body force [N]

    contains

    subroutine initialize_variables()
        integer :: i
        
        do i = 1, nx
            X(i) = (i-1) * dx
            ! n_pos(i) = 0.0d0
            ! n_neg(i) = 0.0d0
            n_pos(i) = 1.0d13*max(exp(- pi*(X(i) - height_flame)**2/a_thickness**2), 0.0)
            n_neg(i) = 1.0d13*max(exp(- pi*(X(i) - height_flame)**2/a_thickness**2), 0.0)
            n_ele(i) = 0.0d0
            rho(i) = (n_pos(i) - n_neg(i) - n_ele(i))*q_e
        end do

        ! set boundary on V
        V(1) = V_start
        V(nx) = V_end

        ! set initial conditions of V
        do i = 2, nx-1
            V(i) = V_start + (V_end - V_start)*((i-1.0)/(nx-1.0))
        end do

        call update_electric_field()

    end subroutine initialize_variables

    subroutine update_electric_field()
        integer :: i

        ! calclate electric field (E = -dV/dx)
        do i = 2, nx-1
        E(i) = -(V(i+1) - V(i-1)) / (2.0d0 * dx)
        end do

        ! boundary conditions
        E(1) = -(V(2) - V(1)) / dx
        E(nx) = -(V(nx) - V(nx-1)) / dx
    end subroutine update_electric_field

    subroutine export_variables(k)
        implicit none
        integer, intent(in) :: k
        integer :: i
        character(len=60) :: filename
        
        ! ! calculate current density
        ! do i = 2, nx-1

        !     current_density(i) = (D_pos*((n_pos(i+1)-n_pos(i-1))/(2.0*dx)) - K_pos*n_pos(i)*E(i))*(+q_e) &
        !                        + (D_pos*((n_neg(i+1)-n_pos(i-1))/(2.0*dx)) + K_neg*n_neg(i)*E(i))*(-q_e) &
        !                        + (D_ele*((n_ele(i+1)-n_pos(i-1))/(2.0*dx)) + K_ele*n_ele(i)*E(i))*(-q_e)

        ! end do

        
        ! create a unique filename using the integer i
        write(filename, '("potential_1d_", I0, ".dat")') k

        print *, "Output to file: ", filename

        ! output
        open(unit=1, file=filename, status='replace')
        write(1,*) "X[m] rho[C/m3] V[V] E[V/m] n_pos[ions/m3] n_neg[ions/m3] n_ele[ions/m3]"
        do i = 1, nx
            write(1, '(7E24.16)') X(i), rho(i), V(i), E(i), n_pos(i), n_neg(i), n_ele(i)
        end do
        close(1)

    end subroutine export_variables

    subroutine import_variables(k)
        implicit none
        integer, intent(in) :: k
        integer :: i
        character(len=60) :: filename
        
        ! ! calculate current density
        ! do i = 2, nx-1

        !     current_density(i) = (D_pos*((n_pos(i+1)-n_pos(i-1))/(2.0*dx)) - K_pos*n_pos(i)*E(i))*(+q_e) &
        !                        + (D_pos*((n_neg(i+1)-n_pos(i-1))/(2.0*dx)) + K_neg*n_neg(i)*E(i))*(-q_e) &
        !                        + (D_ele*((n_ele(i+1)-n_pos(i-1))/(2.0*dx)) + K_ele*n_ele(i)*E(i))*(-q_e)

        ! end do

        
        ! ! create a unique filename using the integer i
        ! write(filename, '("potential_1d_", I0, ".dat")') k
        filename = 'output/1kV_omega_V1.0_omega_ion0.05/potential_1d_100000.dat'
        
        print *, "Import from file: ", filename

        ! open file
        open(unit=1, file=filename, status='old')

        ! skip header
        read(1, '(A)')        
        
        do i = 1, nx
            read(1, '(7E24.16)') X(i), rho(i), V(i), E(i), n_pos(i), n_neg(i), n_ele(i)
        end do

        close(1)

    end subroutine import_variables

end module variables_module
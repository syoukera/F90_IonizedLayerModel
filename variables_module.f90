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
    double precision, parameter :: T = 298d0  ! temperature [K]

    ! diffusion coefficients is drived from Einstein Eq.
    double precision, parameter :: D_pos = K_pos*k_B*T/q_e ! diffusion coefficients of positive ions [m2/s]
    double precision, parameter :: D_neg = K_neg*k_B*T/q_e ! diffusion coefficients of negative ions [m2/s]
    double precision, parameter :: D_ele = K_ele*k_B*T/q_e ! diffusion coefficients of electrons [m2/s]
    
    ! parameters for boundary conditions
    double precision, parameter :: V_start      = 0.0d0 ! valtage for initial point [V]
    double precision, parameter :: V_end        = 1.0d3 ! voltage for end point [V]

    ! parameters for computation
    integer, parameter :: max_iter = 10000
    double precision, parameter :: tolerance = 1.0d-6
    double precision, parameter :: omega_V   = 1.8d0 ! relaxation coefficient (1 < omega < 2)
    double precision, parameter :: omega_pos = 1.8d0 ! relaxation coefficient (1 < omega < 2)
    double precision, parameter :: omega_neg = 1.8d0 ! relaxation coefficient (1 < omega < 2)
    double precision, parameter :: omega_ele = 1.8d0 ! relaxation coefficient (1 < omega < 2)
    double precision :: error

    ! set variables arrays
    double precision :: X(nx) ! position of eac grid point [m]
    double precision :: V(nx) ! electric potential [V]
    double precision :: n_pos(nx) ! number density of positive ions [m-3]
    double precision :: n_neg(nx) ! number density of negative ions [m-3]
    double precision :: n_ele(nx) ! number density of electrons [m-3]
    double precision :: rho(nx) ! density of electric charge [C/m3]

    contains

    subroutine initialize_variables()
        integer :: i
        
        do i = 1, nx
            X(i) = (i-1) * dx
            ! rho(i) = 0.0d0
            ! n_pos(i) = 0.0d0
            n_pos(i) = 1.0d16*max(exp(- pi*(X(i) - height_flame)**2/a_thickness**2), 0.0)
            n_neg(i) = 0.0d0
            n_ele(i) = 0.0d0
        end do

        ! set boundary on V
        V(1) = V_start
        V(nx) = V_end

        ! set initial conditions of V
        do i = 2, nx-1
            V(i) = 0.0d0
        end do
        
    end subroutine

end module variables_module
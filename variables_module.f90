module variables_module
    implicit none

    ! parameters for grid
    integer, parameter :: nx = 101
    double precision, parameter :: dx        = 1.0d-2 ! distance between grid points [m]
    double precision, parameter :: length_x  = (nx-1)*dx ! length of calclation domain [m]
    
    ! parameters for Gaussian profile
    double precision, parameter :: height_flame = length_x/2.0 ! height of flame [m]
    double precision, parameter :: a_thickness  = 1.0d-3 ! thickness parameter [m]
    double precision, parameter :: epsilon_0    = 8.854187817d-12  ! vacuum permittivity [C/V m]
    double precision, parameter :: q_e          = 1.602176634d-19  ! elementary charge [C]
    double precision, parameter :: pi           = 3.141592653589d0
    
    ! parameters for computation
    double precision, parameter :: tolerance = 1.0d-6
    double precision, parameter :: omega     = 1.8d0 ! relaxation coefficient (1 < omega < 2)

    ! set variables arrays
    double precision :: X(nx) ! position of eac grid point [m]
    double precision :: V(nx) ! electric potential [V]
    double precision :: n_p(nx) ! number density of positive ions
    double precision :: n_n(nx) ! number density of negative ions
    double precision :: n_e(nx) ! number density of electrons
    double precision :: rho(nx) ! density of electric charge

    contains

    subroutine initialize_variables()
        integer :: i
        
        do i = 1, nx
            X(i) = (i-1) * dx
            ! rho(i) = 0.0d0
            rho(i) = 1.0d-6*max(exp(- pi*(X(i) - height_flame)**2/a_thickness**2), 0.0)
        end do

        ! set boundary on V
        V(1) = 0.0d0
        V(nx) = 1.0d0

        ! set initial conditions of V
        do i = 1, nx
            V(i) = 0.0d0
        end do
        
    end subroutine

end module variables_module
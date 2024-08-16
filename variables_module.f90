module variables_module
    implicit none
    integer, parameter :: nx = 101
    double precision :: X(nx) ! position of eac grid point [m]
    double precision :: V(nx) ! electric potential [V]
    double precision :: n_p(nx) ! number density of positive ions
    double precision :: n_n(nx) ! number density of negative ions
    double precision :: n_e(nx) ! number density of electrons
    double precision :: rho(nx) ! density of electric charge

    ! パラメータの設定
    double precision, parameter :: dx        = 1.0d-2 ! x方向の格子間隔 [m]
    double precision, parameter :: tolerance = 1.0d-6
    double precision, parameter :: omega     = 1.8d0 ! 過緩和係数 (1 < omega < 2)
    double precision, parameter :: epsilon_0 = 8.854187817d-12  ! vacuum permittivity [C/V m]
    double precision, parameter :: q_e       = 1.602176634d-19  ! elementary charge [C]

    contains

    subroutine initialize_variables()
        integer :: i
        
        do i = 1, nx
            X(i) = (i-1) * dx
            rho(i) = 0.0d0
            ! rho(i) = max(exp(- pi*(x - H)/alpha**2), 0.0)
        end do

        ! 例として、中央に正負の電荷を配置
        rho(nx/2) = 1.0d-6
        rho(nx/2+1) = -1.0d-6
    end subroutine

end module variables_module
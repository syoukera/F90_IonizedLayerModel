module variables_module
    implicit none
    integer, parameter :: nx = 101
    double precision :: V(nx), rho(nx)

    ! パラメータの設定
    double precision, parameter :: dx = 1.0d-2   ! x方向の格子間隔
    double precision, parameter :: tolerance = 1.0d-6
    double precision, parameter :: omega = 1.8d0    ! 過緩和係数 (1 < omega < 2)
    double precision, parameter :: epsilon_0 = 8.854187817d-12  ! 真空の誘電率

end module variables_module
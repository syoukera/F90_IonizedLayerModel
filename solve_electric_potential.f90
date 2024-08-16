program solve_electric_potential
  implicit none
  integer, parameter :: nx = 101
  integer :: i, k, max_iter
  real(8) :: V(nx), rho(nx), V_old(nx), dx, tolerance, omega, error, epsilon_0

  ! パラメータの設定
  dx = 1.0d-2   ! x方向の格子間隔
  max_iter = 10000
  tolerance = 1.0d-6
  omega = 1.8d0    ! 過緩和係数 (1 < omega < 2)
  epsilon_0 = 8.854187817d-12  ! 真空の誘電率

  ! 初期条件と境界条件の設定
  V = 0.0d0
  rho = 0.0d0

  ! 例として、中央に正負の電荷を配置
  rho(nx/2) = 1.0d-6
  rho(nx/2+1) = -1.0d-6

  ! SOR法でポアソン方程式を解く
  do k = 1, max_iter
    V_old = V
    do i = 2, nx-1
      V(i) = (1.0d0 - omega) * V(i) + omega * 0.5d0 * &
             (V(i+1) + V(i-1) + dx*dx*rho(i)/epsilon_0)
    end do

    ! 収束チェック
    error = maxval(abs(V - V_old))
    if (error < tolerance) then
      print *, 'Converged after ', k, ' iterations.'
      exit
    end if
  end do

  if (k == max_iter) then
    print *, 'Did not converge after ', max_iter, ' iterations.'
  end if

  ! 結果を出力 (簡易出力)
  open(unit=1, file='potential_1d.dat', status='replace')
  do i = 1, nx
    write(1,*) i*dx, V(i)
  end do
  close(1)

end program solve_electric_potential

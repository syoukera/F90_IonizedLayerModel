program solve_electric_potential
    use variables_module
    implicit none
    ! integer, parameter :: nx = 101
    integer :: i, k
    integer, parameter :: max_iter = 10000
    ! double precision :: V(nx), rho(nx), V_old(nx), dx, tolerance, omega, error, epsilon_0
    double precision :: V_old(nx), error

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

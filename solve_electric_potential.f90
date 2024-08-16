program solve_electric_potential
    use variables_module
    implicit none
    ! integer, parameter :: nx = 101
    integer :: i, k
    integer, parameter :: max_iter = 10000
    ! double precision :: V(nx), rho(nx), V_old(nx), dx, tolerance, omega, error, epsilon_0
    double precision :: V_old(nx), error
    double precision :: rho(nx) ! density of electric charge [C/m3]

    call initialize_variables()

    ! calclate density of electric charge
    do i = 1, nx
        rho(i) = (n_pos(i) - n_neg(i) - n_neg(i))*q_e
    end do

    ! solve poison equation by SOR method
    do k = 1, max_iter
        V_old = V
        do i = 2, nx-1
            V(i) = (1.0d0 - omega) * V(i) + omega * 0.5d0 * &
                    (V(i+1) + V(i-1) + dx*dx*rho(i)/epsilon_0)
        end do

        ! check convergence
        error = maxval(abs(V - V_old))
        if (error < tolerance) then
            print *, 'Converged after ', k, ' iterations.'
            exit
        end if
    end do

    if (k == max_iter) then
        print *, 'Did not converge after ', max_iter, ' iterations.'
    end if

    ! output
    open(unit=1, file='potential_1d.dat', status='replace')
    write(1,*) "X[m] rho[C/m3] V[V]"
    do i = 1, nx
        write(1,*) i*dx, rho(i), V(i)
    end do
    close(1)

end program solve_electric_potential

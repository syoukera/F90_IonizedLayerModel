program main
    use variables_module
    implicit none
    integer :: i, k

    call initialize_variables()

    ! iteration by SOR method
    do k = 1, max_iter

        error = 0.0

        ! calclate density of electric charge
        do i = 1, nx
            rho(i) = (n_pos(i) - n_neg(i) - n_neg(i))*q_e
        end do

        call solve_poisson_equation()

        call solve_ion_pos_conservation()

        print *, 'step ', i, ' error: ', error

        ! check convergence
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
    write(1,*) "X[m] rho[C/m3] V[V] n_pos[ions/m3] n_neg[ions/m3] n_ele[ions/m3]"
    do i = 1, nx
        write(1,*) i*dx, rho(i), V(i), n_pos(i), n_neg(i), n_ele(i)
    end do
    close(1)

end program main

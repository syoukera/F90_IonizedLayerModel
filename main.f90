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
            rho(i) = (n_pos(i) - n_neg(i) - n_ele(i))*q_e
        end do

        call solve_poisson_equation()

        call update_electric_field()

        call solve_ion_pos_conservation()

        call solve_ion_neg_conservation()

        print *, 'step ', k, ' error: ', error

        ! check convergence
        if (error < tolerance) then
            print *, 'Converged after ', k, ' iterations.'
            exit
        end if

    end do

    if (k == max_iter) then
        print *, 'Did not converge after ', max_iter, ' iterations.'
    end if

    ! calculate current density
    do i = 2, nx-1

        current_density(i) = (D_pos*((n_pos(i+1)-n_pos(i-1))/(2.0*dx)) - K_pos*n_pos(i)*E(i))*(+q_e) &
                           + (D_pos*((n_neg(i+1)-n_pos(i-1))/(2.0*dx)) + K_neg*n_neg(i)*E(i))*(-q_e) &
                           + (D_ele*((n_ele(i+1)-n_pos(i-1))/(2.0*dx)) + K_ele*n_ele(i)*E(i))*(-q_e)

    end do

    ! output
    open(unit=1, file='potential_1d.dat', status='replace')
    write(1,*) "X[m] rho[C/m3] V[V] E[V/m] n_pos[ions/m3] n_neg[ions/m3] n_ele[ions/m3] j[A/m2]"
    do i = 1, nx
        write(1,*) i*dx, rho(i), V(i), E(i), n_pos(i), n_neg(i), n_ele(i), current_density(i)
    end do
    close(1)

end program main

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

        ! ! output_variables for fixed duration
        ! if (mod(k, 1000000) == 0) then
        !     call output_variables(k)
        ! end if

    end do

    if (k == max_iter) then
        print *, 'Did not converge after ', max_iter, ' iterations.'
    end if

    call output_variables(k)

end program main
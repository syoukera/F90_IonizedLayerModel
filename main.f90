program main
    use variables_module
    implicit none
    integer :: i, k

    call initialize_variables()
    ! call import_variables('output/m300V/omega_V0.5_omega_pos0.05/potential_1d_100000.dat')
    
    ! export initial conditions
    call export_variables(k_start)

    ! iteration by SOR method
    do k = k_start, k_end

        error = 0.0

        ! calclate density of electric charge
        do i = 1, nx
            rho(i) = (n_pos(i) - n_neg(i) - n_ele(i))*q_e
        end do

        call solve_poisson_equation()

        call update_electric_field()

        call solve_ion_pos_conservation()

        call solve_ion_neg_conservation()

        call solve_electron_conservation()

        print *, 'step ', k, ' error: ', error

        ! check convergence
        if (error < tolerance) then
            print *, 'Converged after ', k, ' iterations.'
            exit
        end if

        ! export_variables for fixed duration
        if (mod(k, k_step) == 0) then
            call export_variables(k)
        end if

    end do

    if (k == k_end) then
        print *, 'Did not converge after ', k_end, ' iterations.'
    end if

    call export_variables(k)

end program main
program main
    use variables_module
    implicit none
    integer :: k

    call initialize_variables()
    ! call import_variables('output/m300V/omega_V0.5_omega_pos0.05/potential_1d_100000.dat')
    
    ! call load_flame_height()

    call calculate_reaction_profile()

    ! ! export initial conditions
    ! call export_variables()

    ! iteration by SOR method
    do k = k_start, k_end

        error = 0.0

        call update_charge_density()

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

        ! ! export_variables for fixed duration
        ! if (mod(k, k_step) == 0) then
        !     call export_variables()
        ! end if

    end do

    if (k == k_end) then
        print *, 'Did not converge after ', k_end, ' iterations.'
    end if

    call calculate_output_variables()

    call export_variables()

end program main
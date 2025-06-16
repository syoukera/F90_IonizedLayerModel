program main
    use variables_module
    use,intrinsic :: iso_fortran_env
    implicit none
    integer :: k
    double precision :: time_begin_s,time_end_s

    call initialize_variables()
    ! call import_variables('output/m300V/omega_V0.5_omega_pos0.05/potential_1d_100000.dat')
    
    ! call load_flame_height()

    call calculate_reaction_profile()

    ! ! export initial conditions
    ! call export_variables()
    
    call cpu_time(time_begin_s)

    ! iteration by SOR method
    do k = k_start, k_end

        error = 0.0

        call update_charge_density()

        call solve_poisson_equation()

        call update_electric_field()

        call update_source_term()

        call solve_ion_conservation(n_pos, n_pos_old, K_pos, Z_pos, D_pos, Sp_pos, Su_pos, omega_pos)

        call solve_ion_conservation(n_neg, n_neg_old, K_neg, Z_neg, D_neg, Sp_neg, Su_neg, omega_neg)
        
        call solve_ion_conservation(n_ele, n_ele_old, K_ele, Z_ele, D_ele, Sp_ele, Su_ele, omega_ele)

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
    
    call cpu_time(time_end_s)
    print *,"Calculation time: ", time_end_s - time_begin_s,"sec"

    call calculate_output_variables()

    call export_variables()

end program main
program main
    use variables_module
    use,intrinsic :: iso_fortran_env
    implicit none
    ! integer :: step_ion
    double precision :: time_begin_s,time_end_s
    
    double precision :: debye_length_min, debye_length
    integer :: i, j

    call initialize_variables()

    call calculate_reaction_profile()

    ! restart from solution
    ! call import_variables('output')

    call cpu_time(time_begin_s)

    ! iteration by SOR method
    do step = step_start, step_end

        error = 0.0

        call update_charge_density()

        call solve_poisson_equation()
        call update_electric_field()

        call update_source_term()
        call solve_ion_conservation(n_pos, n_pos_old, K_pos, Z_pos, D_pos, Sp_pos, Su_pos, omega_pos)
        
        call update_source_term()
        call solve_ion_conservation(n_neg, n_neg_old, K_neg, Z_neg, D_neg, Sp_neg, Su_neg, omega_neg)
        
        call update_source_term()
        call solve_ion_conservation(n_ele, n_ele_old, K_ele, Z_ele, D_ele, Sp_ele, Su_ele, omega_ele)

        ! Calculate and print the minimum Debye length
        debye_length_min = huge(1.0d0)
        do i = 1, nr
            do j = 1, nz
                ! Debye length: sqrt(epsilon_0 * k_B * T / (n_e * e^2))
                if (n_ele(i, j) > 0.0d0) then
                    debye_length = sqrt(epsilon_0 * k_B * T(i, j) / (n_ele(i, j) * q_e**2))
                    ! debye_length = sqrt(epsilon_0 * k_B * T(i, j) / (n_ele(i, j) * q_e**2))
                    if (debye_length < debye_length_min) debye_length_min = debye_length
                end if
            end do
        end do
        ! print *, 'Minimum Debye length: ', debye_length_min, ' m'

        print *, 'step ', step, ' error: ', error, 'Minimum Debye length: ', debye_length_min, ' m, ', 'dr: ', dr

        ! check convergence
        if (error < tolerance) then
            print *, 'Converged after ', step, ' iterations.'
            exit
        end if

        ! ! export_variables for fixed duration
        ! if (mod(step, step_step) == 0) then
        !     call export_variables()
        ! end if

    end do

    if (step == step_end) then
        print *, 'Did not converge after ', step_end, ' iterations.'
    end if
    
    call cpu_time(time_end_s)
    print *,"Calculation time: ", time_end_s - time_begin_s,"sec"

    call calculate_output_variables()

    call export_variables()

end program main
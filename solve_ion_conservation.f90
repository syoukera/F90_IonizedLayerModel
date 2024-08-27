subroutine solve_ion_pos_conservation
    use variables_module
    implicit none

    ! integer, intent(in) :: nx
    ! real(8), intent(in) :: dx, dt
    ! real(8), dimension(nx), intent(inout) :: n_plus, n_minus, V
    integer :: i
    real(8), dimension(nx) :: n_pos_old 
    real(8) :: g ! spacial profile of ionization
    real(8) :: a, b, c, d ! coefficients of discretised eq.

    ! store old value
    n_pos_old = n_pos

    ! solve coservation equation
    do i = 2, nx-1

        ! calclate spacial profile of ionization
        g = exp(- (pi*(X(i) - height_flame)**2)/a_thickness**2)

        ! upwind difference
        if (E(i) .ge. 0.0) then
            ! calclate coefficients of discretised eq.
            a = 2.0*D_pos/dx**2  + (K_pos/dx)*E(i) + k_r*(n_ele(i) + n_neg(i))
            b = D_pos/dx**2
            c = D_pos/dx**2 + (K_pos/dx)*E(i-1)
        else
            ! calclate coefficients of discretised eq.
            a = 2.0*D_pos/dx**2  - (K_pos/dx)*E(i) + k_r*(n_ele(i) + n_neg(i))
            b = D_pos/dx**2 - (K_pos/dx)*E(i-1)
            c = D_pos/dx**2
        endif

        d = k_i*g

        ! calclate next n_pos(i) using SOR-method
        n_pos(i) = (1.0d0 - omega_pos)*n_pos(i) &
                 + omega_pos*(1.0/a)*(b*n_pos(i+1) + c*n_pos(i-1) + d)

    end do

    error = error + maxval(abs(n_pos - n_pos_old))

end subroutine solve_ion_pos_conservation

subroutine solve_ion_neg_conservation
    use variables_module
    implicit none

    ! integer, intent(in) :: nx
    ! real(8), intent(in) :: dx, dt
    ! real(8), dimension(nx), intent(inout) :: n_plus, n_minus, V
    integer :: i
    real(8), dimension(nx) :: n_neg_old 
    real(8) :: g ! spacial profile of ionization
    real(8) :: a, b, c, d ! coefficients of discretised eq.

    ! store old value
    n_neg_old = n_neg

    ! solve coservation equation
    do i = 2, nx-1

        ! calclate spacial profile of ionization
        g = exp(- (pi*(X(i) - height_flame)**2)/a_thickness**2)

        ! upwind difference
        if (E(i) .ne. 0.0) then
            ! calclate coefficients of discretised eq.
            a = 2.0*D_neg/dx**2  - (K_neg/dx)*E(i) + k_r*n_pos(i)
            b = D_neg/dx**2
            c = D_neg/dx**2 - (K_neg/dx)*E(i-1)
        else
            ! calclate coefficients of discretised eq.
            a = 2.0*D_neg/dx**2  + (K_neg/dx)*E(i) + k_r*n_pos(i)
            b = D_neg/dx**2 + (K_neg/dx)*E(i-1)
            c = D_neg/dx**2
        endif

        d = (1 - alpha)*k_i*g

        ! calclate next n_neg(i) using SOR-method
        n_neg(i) = (1.0d0 - omega_neg)*n_neg(i) &
                 + omega_neg*(1.0/a)*(b*n_neg(i+1) + c*n_neg(i-1) + d)

    end do

    error = error + maxval(abs(n_neg - n_neg_old))

end subroutine solve_ion_neg_conservation
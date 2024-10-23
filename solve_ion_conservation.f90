subroutine solve_ion_pos_conservation
    use variables_module
    implicit none

    ! integer, intent(in) :: nz
    ! double precision, intent(in) :: dz, dt
    ! double precision, dimension(nz), intent(inout) :: n_plus, n_minus, V
    integer :: i, j
    ! double precision :: n_pos_old (nr, nz)
    double precision :: g ! spacial profile of ionization
    double precision :: a, b, c, d ! coefficients of discretised eq.

    ! store old value
    n_pos_old = n_pos

    ! solve coservation equation
    do i = 2, nr-1
        do j = 2, nz-1

            ! calclate spacial profile of ionization
            g = exp(- (pi*(distance_z(i, j) - height_flame)**2)/a_thickness**2)

            ! upwind difference
            if (E_z(i, j) .ge. 0.0) then
                ! calclate coefficients of discretised eq.
                a = 2.0*D_pos/dz**2  + (K_pos/dz)*E_z(i, j) + k_r*(n_ele(i, j) + n_neg(i, j))
                b = D_pos/dz**2
                c = D_pos/dz**2 + (K_pos/dz)*E_z(i, j-1)
            else
                ! calclate coefficients of discretised eq.
                a = 2.0*D_pos/dz**2  - (K_pos/dz)*E_z(i, j) + k_r*(n_ele(i, j) + n_neg(i, j))
                b = D_pos/dz**2 - (K_pos/dz)*E_z(i, j+1)
                c = D_pos/dz**2
            endif

            d = k_i*g

            ! calclate next n_pos(i) using SOR-method
            n_pos(i, j) = (1.0d0 - omega_pos)*n_pos(i, j) &
                    + omega_pos*(1.0/a)*(b*n_pos(i, j+1) + c*n_pos(i, j-1) + d)

        end do
    end do

    error = error + maxval(abs(n_pos - n_pos_old))

end subroutine solve_ion_pos_conservation

subroutine solve_ion_neg_conservation
    use variables_module
    implicit none

    ! integer, intent(in) :: nz
    ! double precision, intent(in) :: dz, dt
    ! double precision, dimension(nz), intent(inout) :: n_plus, n_minus, V
    integer :: i, j
    ! double precision, dimension(nr, nz) :: n_neg_old 
    double precision :: g ! spacial profile of ionization
    double precision :: a, b, c, d ! coefficients of discretised eq.

    ! store old value
    n_neg_old = n_neg

    ! solve coservation equation
    do i = 2, nr-1
        do j = 2, nz-1

            ! calclate spacial profile of ionization
            g = exp(- (pi*(distance_z(i, j) - height_flame)**2)/a_thickness**2)

            ! upwind difference
            if (E_z(i, j) .le. 0.0) then
                ! calclate coefficients of discretised eq.
                a = 2.0*D_neg/dz**2  - (K_neg/dz)*E_z(i, j) + k_r*n_pos(i, j)
                b = D_neg/dz**2
                c = D_neg/dz**2 - (K_neg/dz)*E_z(i, j-1)
            else
                ! calclate coefficients of discretised eq.
                a = 2.0*D_neg/dz**2  + (K_neg/dz)*E_z(i, j) + k_r*n_pos(i, j)
                b = D_neg/dz**2 + (K_neg/dz)*E_z(i, j+1)
                c = D_neg/dz**2
            endif

            d = (1 - alpha)*k_i*g

            ! calclate next n_neg(i) using SOR-method
            n_neg(i, j) = (1.0d0 - omega_neg)*n_neg(i, j) &
                    + omega_neg*(1.0/a)*(b*n_neg(i, j+1) + c*n_neg(i, j-1) + d)

        end do 
    end do

    error = error + maxval(abs(n_neg - n_neg_old))

end subroutine solve_ion_neg_conservation

subroutine solve_electron_conservation
    use variables_module
    implicit none

    ! integer, intent(in) :: nz
    ! double precision, intent(in) :: dz, dt
    ! double precision, dimension(nz), intent(inout) :: n_plus, n_minus, V
    integer :: i, j
    ! double precision, dimension(nr, nz) :: n_ele_old 
    double precision :: g ! spacial profile of ionization
    double precision :: a, b, c, d ! coefficients of discretised eq.

    ! store old value
    n_ele_old = n_ele

    ! solve coservation equation
    do i = 2, nr-1
        do j = 2, nz-1

            ! calclate spacial profile of ionization
            g = exp(- (pi*(distance_z(i, j) - height_flame)**2)/a_thickness**2)

            ! upwind difference
            if (E_z(i, j) .le. 0.0) then
                ! calclate coefficients of discretised eq.
                a = 2.0*D_ele/dz**2  - (K_ele/dz)*E_z(i, j) + k_r*n_pos(i, j)
                b = D_ele/dz**2
                c = D_ele/dz**2 - (K_ele/dz)*E_z(i, j-1)
            else
                ! calclate coefficients of discretised eq.
                a = 2.0*D_ele/dz**2  + (K_ele/dz)*E_z(i, j) + k_r*n_pos(i, j)
                b = D_ele/dz**2 + (K_ele/dz)*E_z(i, j+1)
                c = D_ele/dz**2
            endif

            d = alpha*k_i*g

            ! calclate next n_ele(i) using SOR-method
            n_ele(i, j) = (1.0d0 - omega_ele)*n_ele(i, j) &
                    + omega_ele*(1.0/a)*(b*n_ele(i, j+1) + c*n_ele(i, j-1) + d)

        end do
    end do

    error = error + maxval(abs(n_ele - n_ele_old))

end subroutine solve_electron_conservation
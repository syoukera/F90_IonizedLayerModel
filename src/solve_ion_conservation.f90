subroutine solve_ion_pos_conservation
    use variables_module
    implicit none

    ! integer, intent(in) :: nz
    ! double precision, intent(in) :: dz, dt
    ! double precision, dimension(nz), intent(inout) :: n_plus, n_minus, V
    integer :: i, j
    ! double precision :: n_pos_old (nr, nz)
    double precision :: a, bi, bj, ci, cj, d ! coefficients of discretised eq.

    ! store old value
    n_pos_old = n_pos

    ! solve coservation equation
    do i = 2, nr-1
        do j = 2, nz-1

            ! central difference for diffusion and source term
            a = 4.0*D_pos(i, j)/dz**2  + k_r*n_ele(i, j)
            bi = D_pos(i, j)/dz**2*(1.0d0 + dz/2.0d0/distance_r(i, j))
            ci = D_pos(i, j)/dz**2*(1.0d0 - dz/2.0d0/distance_r(i, j))
            bj = D_pos(i, j)/dz**2
            cj = D_pos(i, j)/dz**2 
            
            ! upwind difference for convection term r direction
            if (E_r(i, j) .ge. 0.0) then
                a = a + (K_pos/dr)*E_r(i, j)
                ci = ci + (K_pos/dr)*(distance_r(i-1, j)/distance_r(i, j))*E_r(i-1, j)
            else
                a = a - (K_pos/dr)*E_r(i, j)
                bi = bi - (K_pos/dr)*(distance_r(i+1, j)/distance_r(i, j))*E_r(i+1, j)
            endif

            ! upwind difference for convection term z direction
            if (E_z(i, j) .ge. 0.0) then
                ! calclate coefficients of discretised eq.
                a = a + (K_pos/dz)*E_z(i, j)
                cj = cj + (K_pos/dz)*E_z(i, j-1)
            else
                ! calclate coefficients of discretised eq.
                a = a - (K_pos/dz)*E_z(i, j)
                bj = bj - (K_pos/dz)*E_z(i, j+1)
            endif

            d = k_i*g_i(i, j)

            ! calclate next n_pos(i) using SOR-method
            n_pos(i, j) = (1.0d0 - omega_pos)*n_pos(i, j) &
                    + omega_pos*(1.0/a)*(bi*n_pos(i+1, j) + ci*n_pos(i-1, j) &
                                       + bj*n_pos(i, j+1) + cj*n_pos(i, j-1) + d)
                        
        end do
    end do
    
    ! boundary condition for z = 0 bottom
    n_pos(:, 1) = n_pos(:, 2) ! (noiman boundary)
    ! n_pos(:, 1) = 0.0d0

    ! boundary condition for z = nz top
    n_pos(:, nz) = n_pos(:, nz-1) ! (noiman boundary)
    ! zero flux on boundary
    ! n_pos(:, nz) = n_pos(:, nz-1)*(1 + K_pos*dz*E_z(:, nz-1)/D_pos(:, nz-1))

    ! boundary condition for r = 0 center axis
    n_pos(1, :) = n_pos(2, :) ! (noiman boundary)

    ! boundary condition for r = nr outside
    n_pos(nr, :) = n_pos(nr-1, :) ! (noiman boundary)

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
    double precision :: a, bi, bj, ci, cj, d ! coefficients of discretised eq.

    ! store old value
    n_neg_old = n_neg

    ! solve coservation equation
    do i = 2, nr-1
        do j = 2, nz-1
            
            ! central difference for diffusion and source term
            ! a = 4.0*D_neg(i, j)/dz**2  + k_r*n_pos(i, j) + k_a*g_a(i, j)
            a = 4.0*D_neg(i, j)/dz**2  + k_r*n_pos(i, j)
            bi = D_neg(i, j)/dz**2*(1.0d0 + dz/2.0d0/distance_r(i, j))
            ci = D_neg(i, j)/dz**2*(1.0d0 - dz/2.0d0/distance_r(i, j))
            bj = D_neg(i, j)/dz**2
            cj = D_neg(i, j)/dz**2 

            ! upwind difference for convection term r direction
            if (E_r(i, j) .le. 0.0) then
                a = a - (K_neg/dr)*E_r(i, j)
                ci = ci - (K_neg/dr)*(distance_r(i-1, j)/distance_r(i, j))*E_r(i-1, j)
            else
                a = a + (K_neg/dr)*E_r(i, j)
                bi = bi + (K_neg/dr)*(distance_r(i+1, j)/distance_r(i, j))*E_r(i+1, j)
            endif

            ! upwind difference for convection term z direction
            if (E_z(i, j) .le. 0.0) then
                ! calclate coefficients of discretised eq.
                a = a  - (K_neg/dz)*E_z(i, j)
                cj = cj - (K_neg/dz)*E_z(i, j-1)
            else
                ! calclate coefficients of discretised eq.
                a = a  + (K_neg/dz)*E_z(i, j)
                bj = bj + (K_neg/dz)*E_z(i, j+1)
            endif
            
            d = k_a*g_a(i, j)*n_ele(i, j)

            ! calclate next n_neg(i) using SOR-method
            n_neg(i, j) = (1.0d0 - omega_neg)*n_neg(i, j) &
                    + omega_neg*(1.0/a)*(bi*n_neg(i+1, j) + ci*n_neg(i-1, j) &
                                       + bj*n_neg(i, j+1) + cj*n_neg(i, j-1) + d)

        end do 
    end do
    
    ! boundary condition for z = 0 bottom
    ! n_neg(:, 1) = n_neg(:, 2) ! (noiman boundary)
    n_neg(:, 1) = 0.0d0

    ! boundary condition for z = nz top
    n_neg(:, nz) = n_neg(:, nz-1) ! (noiman boundary)
    ! zero flux on boundary
    ! n_neg(:, nz) = n_neg(:, nz-1)*(1 - K_neg*dz*E_z(:, nz-1)/D_neg(:, nz-1))

    ! boundary condition for r = 0 center axis
    n_neg(1, :) = n_neg(2, :) ! (noiman boundary)

    ! boundary condition for r = nr outside
    n_neg(nr, :) = n_neg(nr-1, :) ! (noiman boundary)


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
    double precision :: a, bi, bj, ci, cj, d ! coefficients of discretised eq.

    ! store old value
    n_ele_old = n_ele

    ! solve coservation equation
    do i = 2, nr-1
        do j = 2, nz-1
            
            ! central difference for diffusion and source term
            a = 4.0*D_ele(i, j)/dz**2 + k_r*n_pos(i, j) + k_a*g_a(i, j)
            bi = D_ele(i, j)/dz**2*(1.0d0 + dz/2.0d0/distance_r(i, j))
            ci = D_ele(i, j)/dz**2*(1.0d0 - dz/2.0d0/distance_r(i, j))
            bj = D_ele(i, j)/dz**2
            cj = D_ele(i, j)/dz**2 

            ! upwind difference for convection term r direction
            if (E_r(i, j) .le. 0.0) then
                a = a - (K_ele/dr)*E_r(i, j)
                ci = ci - (K_ele/dr)*(distance_r(i-1, j)/distance_r(i, j))*E_r(i-1, j)
            else
                a = a + (K_ele/dr)*E_r(i, j)
                bi = bi + (K_ele/dr)*(distance_r(i+1, j)/distance_r(i, j))*E_r(i+1, j)
            endif
            
            ! upwind difference for convection term z direction
            if (E_z(i, j) .le. 0.0) then
                ! calclate coefficients of discretised eq.
                a = a  - (K_ele/dz)*E_z(i, j)
                cj = cj - (K_ele/dz)*E_z(i, j-1)
            else
                ! calclate coefficients of discretised eq.
                a = a  + (K_ele/dz)*E_z(i, j)
                bj = bj + (K_ele/dz)*E_z(i, j+1)
            endif

            d = k_i*g_i(i, j)

            ! calclate next n_ele(i) using SOR-method
            n_ele(i, j) = (1.0d0 - omega_ele)*n_ele(i, j) &
                    + omega_ele*(1.0/a)*(bi*n_ele(i+1, j) + ci*n_ele(i-1, j) &
                                       + bj*n_ele(i, j+1) + cj*n_ele(i, j-1) + d)

        end do
    end do
    
    ! boundary condition for z = 0 bottom
    ! n_ele(:, 1) = n_ele(:, 2) ! (noiman boundary)
    n_ele(:, 1) = 0.0d0

    ! boundary condition for z = nz top
    n_ele(:, nz) = n_ele(:, nz-1) ! (noiman boundary)
    ! zero flux on boundary
    ! n_ele(:, nz) = n_ele(:, nz-1)*(1 - K_ele*dz*E_z(:, nz-1)/D_ele(:, nz-1))

    ! boundary condition for r = 0 center axis
    n_ele(1, :) = n_ele(2, :) ! (noiman boundary)

    ! boundary condition for r = nr outside
    n_ele(nr, :) = n_ele(nr-1, :) ! (noiman boundary)

    error = error + maxval(abs(n_ele - n_ele_old))

end subroutine solve_electron_conservation
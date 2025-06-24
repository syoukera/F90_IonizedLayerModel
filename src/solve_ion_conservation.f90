subroutine solve_ion_conservation(n_ion, n_ion_old, K_ion, Z_ion, D_ion, Sp_ion, Su_ion, omega_ion)
    use variables_module, only: nr, nz, dz, dr, distance_r, V, error 
                                
    implicit none

    ! integer, intent(in) :: nz
    ! double precision, intent(in) :: dz, dt
    ! double precision, dimension(nz), intent(inout) :: n_plus, n_minus, V
    integer :: i, j
    ! double precision :: n_ion_old (nr, nz)
    double precision :: a, bi, bj, ci, cj, d ! coefficients of discretised eq.
    double precision :: ddVdr, ddVdz ! 2nd derivetive of voltage
    double precision :: r_p, r_e, r_w ! r, z distance on center and mean value for North, South, West, East
    double precision :: E_r_e, E_r_w, E_z_n, E_z_s, E_r_p ! mean value of E

    double precision, intent(inout) :: n_ion(nr, nz)
    double precision, intent(inout) :: n_ion_old(nr, nz)
    double precision, intent(in) :: D_ion(nr, nz)
    double precision, intent(in) :: Sp_ion(nr, nz)
    double precision, intent(in) :: Su_ion(nr, nz)
    double precision, intent(in) :: K_ion
    double precision, intent(in) :: Z_ion
    double precision, intent(in) :: omega_ion

    ! store old value
    n_ion_old = n_ion

    ! solve coservation equation
    do i = 2, nr-1
        do j = 2, nz-1

            ! prepare man value of r
            r_p = distance_r(i, j)
            r_e = (distance_r(i+1, j)+ distance_r(i, j))/2.0
            r_w = (distance_r(i-1, j)+ distance_r(i, j))/2.0

            ! prepare mean value of E
            E_r_e = - (V(i+1, j) - V(i  , j))/dr
            E_r_w = - (V(i  , j) - V(i-1, j))/dr
            E_r_p = (E_r_e + E_r_w)/2.0
            E_z_n = - (V(i, j+1) - V(i, j  ))/dz
            E_z_s = - (V(i, j  ) - V(i, j-1))/dz

            ddVdr = (V(i+1, j) - 2.0*V(i, j) + V(i-1, j))/(dr**2)
            ddVdz = (V(i, j+1) - 2.0*V(i, j) + V(i, j-1))/(dz**2)

            ! central difference for diffusion and source term
            a = 2.0*D_ion(i, j)/(dz**2) + D_ion(i, j)/(r_p*dr**2)*(r_e + r_w) - Sp_ion(i, j)
            bi = D_ion(i, j)/(r_p*dr**2)*r_e
            ci = D_ion(i, j)/(r_p*dr**2)*r_w
            bj = D_ion(i, j)/(dz**2)
            cj = D_ion(i, j)/(dz**2)

            ! upwind difference for convection term r direction
            a = a + (K_ion/(r_p*dr)) * (r_e*max(Z_ion*E_r_e, 0.0) - r_w*min(Z_ion*E_r_w, 0.0)) &
              - K_ion*ddVdr + (K_ion/r_p)*Z_ion*E_r_p
            !   - K_ion*ddVdr + (K_ion/r_p)*Z_ion*E_r_p
            bi = bi - (K_ion/(r_p*dr)) * r_e * min(Z_ion*E_r_e, 0.0)
            ci = ci + (K_ion/(r_p*dr)) * r_w * max(Z_ion*E_r_w, 0.0)

            ! upwind difference for convection term z direction
            a = a  + (K_ion/dz) * (max(Z_ion*E_z_n, 0.0) - min(Z_ion*E_z_s, 0.0)) &
              - K_ion*ddVdz
            bj = bj - (K_ion/dz) * min(Z_ion*E_z_n, 0.0)
            cj = cj + (K_ion/dz) * max(Z_ion*E_z_s, 0.0)

            d = Su_ion(i, j)

            ! calclate next n_ion(i) using SOR-method
            n_ion(i, j) = (1.0d0 - omega_ion)*n_ion(i, j) &
                    + omega_ion*(1.0/a)*(bi*n_ion(i+1, j) + ci*n_ion(i-1, j) &
                                       + bj*n_ion(i, j+1) + cj*n_ion(i, j-1) + d)
                        
        end do
    end do
    

    ! boundary condition for z = 0 bottom
    ! n_ion(:, 1) = n_ion(:, 2) ! (noiman boundary)
    ! n_ion(:, 1) = 0.0d0
    ! Table 1 of Yihua Ren
    do i = 2, nr-1

        E_z_n = - (V(i, 2) - V(i, 1))/dz

        if (Z_ion*E_z_n > 0.0) then
            ! inflow flux equals zero
            n_ion(i, 1) = n_ion(i, 2)*(1.0/(1.0 + K_ion*Z_ion*E_z_n*dz/D_ion(i, 1)))
        else
            ! inflow flux from electric field
            n_ion(i, 1) = n_ion(i, 2) - Su_ion(i, 1)*dz/(K_ion*Z_ion*E_z_n)
        end if
    end do
    
    ! boundary condition for z = nz top
    ! n_ion(:, nz) = n_ion(:, nz-1) ! (noiman boundary)
    ! zero flux on boundary
    ! n_ion(:, nz) = n_ion(:, nz-1)*(1 + K_ion*dz*E_z(:, nz-1)/D_ion(:, nz-1))
    ! Table 1 of Yihua Ren
    do i = 2, nr-1
        
        E_z_s = - (V(i, nr) - V(i, nr-1))/dz
    
        if (Z_ion*E_z_s > 0.0) then
            ! inflow flux from electric field
            n_ion(i, nz) = n_ion(i, nz-1) + Su_ion(i, nz)*dz/(K_ion*Z_ion*E_z_s)
        else
            ! inflow flux equals zero
            n_ion(i, nz) = n_ion(i, nz-1)*(1.0/(1.0 - K_ion*Z_ion*E_z_s*dz/D_ion(i, nz)))
        end if
    end do

    ! boundary condition for r = 0 center axis
    ! n_ion(1, :) = n_ion(2, :) ! (noiman boundary)
    ! Table 1 of Yihua Ren
    do j = 2, nz-1
        ! 2nd derivetive of voltage
        ddVdr = (2.0*V(1, j) - 5.0*V(2, j) + 4.0*V(3, j) - V(4, j))/(dr**2)

        ! update n_ion
        n_ion(1, j) = (5.0*n_ion(2, j) - 4.0*n_ion(3, j) + n_ion(4, j))/2.0 &
                    + (-ddVdr)*K_ion*Z_ion*(dr)**2/(2.0*D_ion(1, j))
    end do

    ! boundary condition for r = nr outside
    ! n_ion(nr, :) = n_ion(nr-1, :) ! (noiman boundary)
    ! Table 1 of Yihua Ren
    do j = 2, nz-1
        ! 2nd derivetive of voltage
        ddVdr = (2.0*V(nr, j)- 5.0*V(nr-1, j) + 4.0*V(nr-2, j)  -V(nr-3, j))/(dr**2)

        ! update n_ion
        n_ion(nr, j) = (5.0*n_ion(nr-1, j) - 4.0*n_ion(nr-2, j) + n_ion(nr-3, j))/2.0  &
                    + (-ddVdr)*K_ion*Z_ion*(dr**2)/(2.0*D_ion(nr, j))
    end do

    error = error + maxval(abs(n_ion - n_ion_old))

end subroutine solve_ion_conservation

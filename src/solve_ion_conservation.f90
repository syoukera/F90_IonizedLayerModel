subroutine solve_ion_conservation(n_ion, n_ion_old, K_ion, Z_ion, D_ion, Sp_ion, Su_ion, omega_ion)
    use variables_module, only: nr, nz, dz, dr, distance_r, V, error, ddVdr_prev, ddVdz_prev, k_B, q_e, T, case_convection
    implicit none

    ! integer, intent(in) :: nz
    ! double precision, intent(in) :: dz, dt
    ! double precision, dimension(nz), intent(inout) :: n_plus, n_minus, V
    integer :: i, j
    ! double precision :: n_ion_old (nr, nz)
    double precision :: a, bi, bj, ci, cj, d ! coefficients of discretised eq.
    double precision :: ddVdr, dVdr, ddVdz, dVdz ! 2nd derivetive of voltage
    double precision :: r_p, r_e, r_w ! r, z distance on center and mean value for North, South, West, East
    double precision :: E_r_e, E_r_w, E_z_n, E_z_s, E_r_p, E_z_p ! mean value of E
    
    ! SG method
    double precision :: Vt
    double precision :: delta_r_e, delta_r_w, delta_z_n, delta_z_s
    double precision :: coeff_r, coeff_z

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
    !$omp parallel do collapse(2) &
    !$omp private(i, j, r_p, r_e, r_w, E_r_e, E_r_w, E_r_p, E_z_n, E_z_s, E_z_p) &
    !$omp private(ddVdr, dVdz, a, bi, bj, ci, cj, d, Vt) &
    !$omp private(delta_r_e, delta_r_w, delta_z_n, delta_z_s, coeff_r, coeff_z)
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
            E_z_p = (E_z_n + E_z_s)/2.0

            ddVdr = ddVdr_prev(i, j)
            ddVdz = ddVdz_prev(i, j)
            ! ddVdr = (V(i+1, j) - 2.0*V(i, j) + V(i-1, j))/(dr**2)
            ! ddVdz = (V(i, j+1) - 2.0*V(i, j) + V(i, j-1))/(dz**2)

            ! central difference for diffusion and source term
            a = 2.0*D_ion(i, j)/(dz**2) + D_ion(i, j)/(r_p*dr**2)*(r_e + r_w) - Sp_ion(i, j)
            bi = D_ion(i, j)/(r_p*dr**2)*r_e
            ci = D_ion(i, j)/(r_p*dr**2)*r_w
            bj = D_ion(i, j)/(dz**2)
            cj = D_ion(i, j)/(dz**2)

            ! ! central difference for diffusion and source term (takuma)
            ! a = 4.0*D_ion(i, j)/(dz**2) - Sp_ion(i, j)
            ! bi = D_ion(i, j)*(1.0/(dr**2) + 1.0/(2.0*r_p*dr))
            ! ci = D_ion(i, j)*(1.0/(dr**2) - 1.0/(2.0*r_p*dr))
            ! bj = D_ion(i, j)/(dr**2)
            ! cj = D_ion(i, j)/(dr**2)

            ! a = 0.0d0
            ! bi = 0.0d0
            ! ci = 0.0d0
            ! bj = 0.0d0
            ! cj = 0.0d0

            select case (case_convection)
                case (1) ! upwin discretize

                    ! upwind difference    for convection term r direction
                    a = a + (K_ion/(r_p*dr)) * (r_e*max(Z_ion*E_r_e, 0.0) - r_w*min(Z_ion*E_r_w, 0.0)) &
                    !   + (K_ion/r_p)*Z_ion*E_r_p
                      - K_ion*ddVdr + (K_ion/r_p)*Z_ion*E_r_p
                    bi = bi - (K_ion/(r_p*dr)) * r_e * min(Z_ion*E_r_e, 0.0)
                    ci = ci + (K_ion/(r_p*dr)) * r_w * max(Z_ion*E_r_w, 0.0)

                    ! upwind difference for convection term z direction
                    a = a  + (K_ion/dz) * (max(Z_ion*E_z_n, 0.0) - min(Z_ion*E_z_s, 0.0)) &
                      - K_ion*ddVdz
                    bj = bj - (K_ion/dz) * min(Z_ion*E_z_n, 0.0)
                    cj = cj + (K_ion/dz) * max(Z_ion*E_z_s, 0.0)

                case (2) ! SG method
                    
                    ! Vt: 熱電圧 (kBT/e)。拡散係数と移動度の比 D/K から計算
                    Vt = k_B * T(i, j) / q_e  ! または Vt = D_ion / K_ion
                    
                    ! 各面での無次元電位差 Delta = (Z*E*dL)/Vt を計算
                    ! ※ Z_ion*E はポテンシャルの勾配に対応
                    delta_r_e = (Z_ion * E_r_e * dr) / Vt
                    delta_r_w = (Z_ion * E_r_w * dr) / Vt
                    delta_z_n = (Z_ion * E_z_n * dz) / Vt
                    delta_z_s = (Z_ion * E_z_s * dz) / Vt
                
                    ! --- r方向のフラックス係数 (軸対称) ---
                    ! a: 中心(i,j), bi: 東(i+1,j), ci: 西(i-1,j)
                    ! 拡散項もこの中に含まれるため、元のddVdrなどの項は不要になります。
                    
                    coeff_r = (K_ion * Vt) / (r_p * dr**2)
                    
                    ! 東側(e)面からの寄与
                    a  = a  + coeff_r * r_e * bernoulli(-delta_r_e)
                    bi = bi - coeff_r * r_e * bernoulli(delta_r_e)
                    
                    ! 西側(w)面からの寄与
                    a  = a  + coeff_r * r_w * bernoulli(delta_r_w)
                    ci = ci - coeff_r * r_w * bernoulli(-delta_r_w)

                    ! --- z方向のフラックス係数 ---
                    ! a: 中心(i,j), bj: 北(i,j+1), cj: 南(i,j-1)
                    
                    coeff_z = (K_ion * Vt) / (dz**2)
                    
                    ! 北側(n)面からの寄与
                    a  = a  + coeff_z * bernoulli(-delta_z_n)
                    bj = bj - coeff_z * bernoulli(delta_z_n)
                    
                    ! 南側(s)面からの寄与
                    a  = a  + coeff_z * bernoulli(delta_z_s)
                    cj = cj - coeff_z * bernoulli(-delta_z_s)

            end select

            ! ソース項
            d = Su_ion(i, j)

            ! SOR法で次のn_ionを計算
            n_ion(i, j) = (1.0d0 - omega_ion)*n_ion(i, j) &
                    + omega_ion*(1.0/a)*(bi*n_ion(i+1, j) + ci*n_ion(i-1, j) &
                                       + bj*n_ion(i, j+1) + cj*n_ion(i, j-1) + d)
                        
        end do
    end do
    !$omp end parallel do
    

    ! boundary condition for z = 0 bottom
    ! n_ion(:, 1) = n_ion(:, 2) ! (noiman boundary)
    ! n_ion(:, 1) = 0.0d0
    ! Table 1 of Yihua Ren
    do i = 2, nr-1

        dVdz = (-3.0*V(i, 1) + 4.0*V(i, 2) - V(i, 3))/(2.0*dz)

        if (- Z_ion*dVdz < 0.0) then
            ! inflow flux equals zero
            n_ion(i, 1) = (4.0*n_ion(i, 2) - n_ion(i, 3))/(3.0 - (2.0*K_ion*Z_ion*dz/D_ion(i, 1))*dVdz)
        else
            ! inflow flux from electric field
            n_ion(i, 1) = n_ion(i, 2)
        end if

    end do

    ! boundary condition for z = nz top
    ! n_ion(:, nz) = n_ion(:, nz-1) ! (noiman boundary)
    ! zero flux on boundary
    ! n_ion(:, nz) = n_ion(:, nz-1)*(1 + K_ion*dz*E_z(:, nz-1)/D_ion(:, nz-1))
    ! Table 1 of Yihua Ren
    do i = 2, nr-1

        dVdz = (V(i, nz-2) - 4.0*V(i, nz-1) + 3.0*V(i, nz))/(2.0*dz)

        if (- Z_ion*dVdz < 0.0) then
            ! inflow flux from electric field
            n_ion(i, nz) = n_ion(i, nz-1)
        else
            ! inflow flux equals zero
            n_ion(i, nz) = (4.0*n_ion(i, nz-1) - n_ion(i, nz-2))/(3.0 + (2.0*K_ion*Z_ion*dz/D_ion(i, nz))*dVdz)
        end if

    end do

    ! boundary condition for r = 0 center axis
    n_ion(1, :) = n_ion(2, :) ! (noiman boundary)
    ! ! Table 1 of Yihua Ren
    ! do j = 2, nz-1

    !     dVdr = (-3.0*V(1, j) + 4.0*V(2, j) - V(3, j))/(2.0*dr)
    !     ddVdr = (2.0*V(1, j) - 5.0*V(2, j) + 4.0*V(3, j) - V(4, j))/(dr**2)

    !     n_ion(1, j) = D_ion(1, j)*(-5.0*n_ion(2, j) + 4.0*n_ion(3, j) - n_ion(4, j))/(dr**2) &
    !                 + K_ion*Z_ion*(4.0*n_ion(2, j) - n_ion(3, j))/(2.0*dr)

    !     n_ion(1, j) = n_ion(1, j)/(-2.0*D_ion(1, j)/(dr**2) - K_ion*Z_ion*(- 3.0/(2.0*dr)*dVdr + ddVdr))

    ! end do 

    ! boundary condition for r = nr outside
    n_ion(nr, :) = n_ion(nr-1, :) ! (noiman boundary)
    ! Table 1 of Yihua Ren
    ! do j = 2, nz-1

    !     dVdr = (V(nr-2, j) - 4.0*V(nr-1, j) + 3.0*V(nr, j))/(2.0*dr)
    !     ddVdr = (-V(nr-3, j) + 4.0*V(nr-2, j) - 5.0*V(nr-1, j) + 2.0*V(nr, j))/(dr**2)

    !     n_ion(nr, j) = D_ion(nr, j)*(-n_ion(nr-3, j) + 4.0*n_ion(nr-2, j) - 5.0*n_ion(nr-1, j))/(dr**2) &
    !                  + K_ion*Z_ion*(n_ion(nr-2, j) - 4.0*n_ion(nr-1, j))/(2.0*dr)

    !     n_ion(nr, j) = n_ion(nr, j)/(- 2.0*D_ion(nr, j)/(dr**2) - K_ion*Z_ion*(3.0/(2.0*dr)*dVdr + ddVdr))

    ! end do

    error = error + maxval(abs(n_ion - n_ion_old))

contains

    ! ベルヌーイ関数の定義 (内部関数として定義)
    pure function bernoulli(x) result(res)
        real(8), intent(in) :: x
        real(8) :: res
        if (abs(x) < 1.0d-4) then
            res = 1.0d0 - 0.5d0*x + (1.0d0/12.0d0)*x**2
        else
            res = x / (exp(x) - 1.0d0)
        end if
    end function bernoulli

end subroutine solve_ion_conservation

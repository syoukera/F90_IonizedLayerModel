subroutine solve_poisson_equation
    use variables_module
    implicit none
    integer :: i, j
    ! double precision :: V_old(nr, nz)
    double precision :: r
    double precision :: Vt_ele, Vt_pos, Vt_neg, alpha, rho_current, V_sum, V_source, V_correction, denom, V_target

    ! save old value
    V_old = V

    ! ! solve poisson_equation
    ! do i = 2, nz-1
    !     V(i) = (1.0d0 - omega_V)*V(i) + omega_V*0.5d0* &
    !             (V(i+1) + V(i-1) + dz*dz*rho(i)/epsilon_0)
    ! end do
    
    ! ! solve poisson_equation
    ! !$omp parallel do collapse(2) private(i, j, r)
    ! do i = 2, nr-1
    !     do j = 2, nz-1

    !         r = distance_r(i, j)

    !         ! cylindrical grid
    !         V(i, j) = (1 - omega_V) * V(i, j) + omega_V * 0.25D0 * &
    !         (V(i+1, j)*(1 + 0.5*dr/r) + V(i-1, j)*(1.0 - 0.5*dr/r) + V(i, j+1) + V(i, j-1) + dr*dr*rho(i, j)/epsilon_0)

    !     end do
    ! end do
    ! !$omp end parallel do

    ! solve poisson_equation by linearize
    do i = 2, nr-1
        do j = 2, nz-1
            r = distance_r(i, j)

            select case (case_poison)
                case (1) ! direct

                    ! cylindrical grid
                    V(i, j) = (1 - omega_V) * V(i, j) + omega_V * 0.25D0 * &
                    (V(i+1, j)*(1 + 0.5*dr/r) + V(i-1, j)*(1.0 - 0.5*dr/r) + V(i, j+1) + V(i, j-1) + dr*dr*rho(i, j)/epsilon_0)


                case (2) ! linearlize

                    ! --- 事前に計算しておく定数 ---
                    Vt_ele = k_B * T(i, j) / q_e !  (電子の熱電圧)
                    Vt_pos = k_B * T(i, j) / q_e !  (陽イオンの熱電圧)
                    Vt_neg = k_B * T(i, j) / q_e !  (陰イオンの熱電圧)

                    ! 1. シールド係数 (alpha) の計算
                    !    これはデバイ長の逆数の2乗 (1/lambda_D^2) に相当します。
                    !    電荷密度が高いほど alpha が大きくなり、計算が安定します。
                    alpha = (q_e / epsilon_0) * ( &
                            n_ele(i,j) / Vt_ele + &
                            n_pos(i,j) / Vt_pos + &
                            n_neg(i,j) / Vt_neg )

                    ! 2. 現在の電荷密度 rho = e(n_pos - n_neg - n_ele)
                    rho_current = q_e * (n_pos(i,j) - n_neg(i,j) - n_ele(i,j))

                    ! 3. 線形化されたポアソン方程式の適用
                    !    通常の項に加え、左辺に alpha * V_new、右辺に alpha * V_old を考慮した形に整理します。
                    
                    ! 分子の計算 (周辺格子の和 + ソース項 + 線形化補正項)
                    V_sum = V(i+1, j)*(1.0d0 + 0.5d0*dr/r) + &
                            V(i-1, j)*(1.0d0 - 0.5d0*dr/r) + &
                            V(i, j+1) + V(i, j-1)
                    
                    V_source = (dr*dr / epsilon_0) * rho_current
                    V_correction = dr*dr * alpha * V(i, j) ! 前のイテレーションのVを使用

                    ! 分母の計算 (通常は 4.0 だが、線形化により対角成分が強化される)
                    denom = 4.0d0 + dr*dr * alpha

                    ! SOR更新
                    V_target = (V_sum + V_source + V_correction) / denom
                    V(i, j) = (1.0d0 - omega_V) * V(i, j) + omega_V * V_target

            end select

        end do
    end do
    
    ! set boundary of z on V
    do i = 1, nr
        V(i, 1) = V_start
        V(i, nz) = V_end
    end do

    ! center axis (i = 1) 
    V(1, :) = V(2, :) ! (noiman boundary)
    ! ! Takuma (3.20)
    ! do j = 2, nz-1
    !     V(1, j) = (1/6.0D0)*(V(1, j+1) + 4.0D0*V(2, j) + V(1, j-1))
    ! end do

    ! outlet (i = nr)
    V(nr, :) = V(nr-1, :) ! (noiman boundary)
    ! ! Takuma (3.20)
    ! do j = 2, nz-1
    !     V(nr, j) = (1/4.0D0)*(V(nr, j+1) + 2.0D0*V(nr-1, j) + V(nr, j-1))
    !     ! V(nr, j) = V_end/10
    ! end do

    ! calclate error for check convergence
    error = error + maxval(abs(V - V_old))

end subroutine solve_poisson_equation

subroutine solve_poisson_equation
    use variables_module
    implicit none
    integer :: i, j
    double precision :: V_old(nr, nz)
    double precision :: r

    ! save old value
    V_old = V

    ! ! solve poisson_equation
    ! do i = 2, nz-1
    !     V(i) = (1.0d0 - omega_V)*V(i) + omega_V*0.5d0* &
    !             (V(i+1) + V(i-1) + dz*dz*rho(i)/epsilon_0)
    ! end do
    
    ! solve poisson_equation
    do i = 2, nr-1
        do j = 2, nz-1

            r = (i - 1) * dr

            ! cylindrical grid
            V(i, j) = (1 - omega_V) * V(i, j) + omega_V * 0.25D0 * &
            (V(i+1, j)*(1 + 0.5*dr/r) + V(i-1, j)*(1.0 - 0.5*dr/r) + V(i, j+1) + V(i, j-1))

        end do
    end do

    ! calclate error for check convergence
    error = error + maxval(abs(V - V_old))

end subroutine solve_poisson_equation

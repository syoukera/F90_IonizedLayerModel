subroutine solve_poisson_equation
    use variables_module
    implicit none
    integer :: i, j
    ! double precision :: V_old(nr, nz)
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

            r = distance_r(i, j)

            ! cylindrical grid
            V(i, j) = (1 - omega_V) * V(i, j) + omega_V * 0.25D0 * &
            (V(i+1, j)*(1 + 0.5*dr/r) + V(i-1, j)*(1.0 - 0.5*dr/r) + V(i, j+1) + V(i, j-1) + dr*dr*rho(i, j)/epsilon_0)

        end do
    end do

    ! center axis (i = 1) 
    ! dV/dr = 0.0
    do j = 2, nz-1
        V(1, j) = (1/6.0D0)*(V(1, j+1) + 4.0D0*V(2, j) + V(1, j-1))
    end do

    ! outlet (i = nr)
    ! dV/dr = 0.0
    do j = 2, nz-1
        V(nr, j) = (1/4.0D0)*(V(nr, j+1) + 2.0D0*V(nr-1, j) + V(nr, j-1))
        ! V(nr, j) = V_end/10
    end do

    ! calclate error for check convergence
    error = error + maxval(abs(V - V_old))

end subroutine solve_poisson_equation

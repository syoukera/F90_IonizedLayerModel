subroutine solve_poisson_equation
    use variables_module
    implicit none
    integer :: i
    double precision :: V_old(nx)

    ! save old value
    V_old = V

    ! solve poisson_equation
    do i = 2, nx-1
        V(i) = (1.0d0 - omega_V)*V(i) + omega_V*0.5d0* &
                (V(i+1) + V(i-1) + dx*dx*rho(i)/epsilon_0)
    end do

    ! calclate error for check convergence
    error = error + maxval(abs(V - V_old))

end subroutine solve_poisson_equation

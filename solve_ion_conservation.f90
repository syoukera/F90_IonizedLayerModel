subroutine solve_ion_pos_conservation
    use variables_module
    implicit none

    ! integer, intent(in) :: nx
    ! real(8), intent(in) :: dx, dt
    ! real(8), dimension(nx), intent(inout) :: n_plus, n_minus, V
    integer :: i
    real(8), dimension(nx) :: n_pos_old, n_pos_old

    n_pos_old = n_pos
    n_neg_old = n_neg

    ! イオン保存式の解法（陽解法の例）
    do i = 2, nx-1
        n_plus(i) = n_plus_old(i) - dt/dx * (n_plus_old(i) * (V(i+1) - V(i-1)) / (2.0d0*dx))
        n_minus(i) = n_minus_old(i) + dt/dx * (n_minus_old(i) * (V(i+1) - V(i-1)) / (2.0d0*dx))
    end do

end subroutine solve_ion_pos_conservation
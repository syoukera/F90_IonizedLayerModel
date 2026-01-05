module variables_module
    implicit none

    ! parameters for grid
    integer, parameter :: nr = 401
    integer, parameter :: nz = nr

    ! Note: dr = dz must be preserved in current imprementation    
    double precision, parameter :: length_r = 25d-3 ! length of calclation domain [m] 
    double precision, parameter :: length_z = 25d-3 ! length of calclation domain [m] 
    double precision, parameter :: dr = length_r/(nr - 1.0) ! distance between grid points [m]
    double precision, parameter :: dz = length_z/(nz - 1.0) ! distance between grid points [m]

    ! constants
    double precision, parameter :: epsilon_0    = 8.854187817d-12  ! vacuum permittivity [C/V m]
    double precision, parameter :: q_e          = 1.602176634d-19  ! elementary charge [C]
    double precision, parameter :: pi           = 3.141592653589d0
    double precision, parameter :: k_B          = 1.380649d-23 ! Boltzmann constant [J/K] 

    ! parameters for Gaussian profile
    double precision, parameter :: height_flame = length_z/2.0 ! height of flame [m]
    double precision, parameter :: a_thickness  = 1.0d-3 ! thickness parameter [m]

    ! parameters for transport and reactions
    double precision, parameter :: k_i = 1.76780381d+21 ! rate coeficient of ionization ions/m3/s
    double precision, parameter :: k_r = 1.89301454d-13 ! rate coeficient of recombination m3/ions s
    double precision, parameter :: k_a = 4.73873934d+07 ! rate coeficient of attachment  1/s
    double precision, parameter :: K_pos = 2.9d-4 ! mobility of positive ions [m2/s V]
    double precision, parameter :: K_neg = 2.9d-4 ! mobility of negative ions [m2/s V]
    double precision, parameter :: K_ele = 0.4    ! mobility of electrons [m2/s V]
    double precision, parameter :: Z_pos = 1.0    ! polarity of positive ion's charge [-]
    double precision, parameter :: Z_neg = - 1.0  ! polarity of negative ion's charge [-]
    double precision, parameter :: Z_ele = - 1.0  ! polarity of electron's charge [-]
    
    ! parameters for boundary conditions
    double precision, parameter :: V_start      = 0.0d0 ! valtage for initial point [V]
    double precision, parameter :: V_end        = -1.2d3 ! voltage for end point [V]

    ! parameters for computation
    integer :: step
    integer, parameter :: step_start = 1
    integer, parameter :: step_end   = 1000
    integer, parameter :: step_step  = 100000
    double precision, parameter :: tolerance = 1.0d-10

    double precision, parameter :: omega_V   = 0.5 ! relaxation coefficient (1 < omega < 2)
    double precision, parameter :: omega_pos = 0.05 ! relaxation coefficient (1 < omega < 2)
    double precision, parameter :: omega_neg = 0.05 ! relaxation coefficient (1 < omega < 2)
    double precision, parameter :: omega_ele = 0.05 ! relaxation coefficient (1 < omega < 2)

    integer, parameter :: case_convection = 2 ! 1 for upwind, 2 for SG method
    integer, parameter :: case_poison = 1 ! 1 for direct, 2 for linearized way

    double precision :: error

    ! set variable arrays
    double precision :: distance_r(nr, nz) ! position of eac grid point [m]
    double precision :: distance_z(nr, nz) ! position of eac grid point [m]
    double precision :: V(nr, nz) ! electric potential [V]
    double precision :: E_r(nr, nz) ! electric field [V/m]
    double precision :: E_z(nr, nz) ! electric field [V/m]
    double precision :: n_pos(nr, nz) ! number density of positive ions [m-3]
    double precision :: n_neg(nr, nz) ! number density of negative ions [m-3]
    double precision :: n_ele(nr, nz) ! number density of electrons [m-3]
    double precision :: rho(nr, nz) ! density of electric charge [C/m3]
    
    double precision :: V_old(nr, nz)
    double precision :: n_pos_old (nr, nz)
    double precision :: n_neg_old (nr, nz)
    double precision :: n_ele_old (nr, nz)

    double precision :: g_i(nr, nz) ! spacial profile of ionization
    double precision :: g_a(nr, nz) ! spacial profile of attachement
    double precision :: T(nr, nz)  ! temperature [K]

    ! diffusion coefficients is drived from Einstein Eq.
    double precision :: D_pos(nr, nz) ! diffusion coefficients of positive ions [m2/s]
    double precision :: D_neg(nr, nz) ! diffusion coefficients of negative ions [m2/s]
    double precision :: D_ele(nr, nz) ! diffusion coefficients of electrons [m2/s]
    
    ! array for source term
    double precision :: Sp_pos(nr, nz) ! source term of positive ion dependent of n
    double precision :: Su_pos(nr, nz) ! source term of positive ion independent of n
    double precision :: Sp_neg(nr, nz) ! source term of negative ion dependent of n
    double precision :: Su_neg(nr, nz) ! source term of negative ion independent of n
    double precision :: Sp_ele(nr, nz) ! source term of electron dependent of n
    double precision :: Su_ele(nr, nz) ! source term of electron independent of n
    
    ! output variable arrays
    double precision :: J_r(nr, nz) ! current density [C/m2]
    double precision :: J_z(nr, nz) ! current density [C/m2]
    double precision :: F_r(nr, nz) ! electric body force [N/m3]
    double precision :: F_z(nr, nz) ! electric body force [N/m3]

    ! input data for flame height
    double precision :: normalized_flame_height(nr)
    double precision :: normalized_intensity(nr)

    ! モジュール変数またはサブルーチン内での定義
    double precision, save :: ddVdz_prev(nr, nz) = 0.0d0
    double precision, save :: ddVdr_prev(nr, nz) = 0.0d0
    double precision       :: alpha_E = 0.05d0  ! 電界緩和係数 (非常に小さく設定)                            

contains

    subroutine initialize_variables()
        integer :: i, j
        
        do i = 1, nr
            do j = 1, nz

                distance_r(i, j) = (i-1) * dr 
                distance_z(i, j) = (j-1) * dz

                n_pos(i, j) = 0.0d0
                ! n_pos(i, j) = 1.0d13*max(exp(- pi*(distance_z(i, j) - height_flame)**2/a_thickness**2), 0.0)
                n_neg(i, j) = 0.0d0
                ! n_neg(i, j) = 1.0d13*max(exp(- pi*(distance_z(i, j) - height_flame)**2/a_thickness**2), 0.0)
                n_ele(i, j) = 0.0d0
                ! n_ele(i, j) = 1.0d13*max(exp(- pi*(distance_z(i, j) - height_flame)**2/a_thickness**2), 0.0)

            end do
        end do

        ! set boundary on V
        do i = 1, nr
            V(i, 1) = V_start
            V(i, nz) = V_end
        end do

        ! set initial conditions of V
        do i = 1, nr
            do j = 2, nz-1
                ! V(i, j) = 0.0d0
                V(i, j) = V_start + (V_end - V_start)*((j-1.0)/(nz-1.0))
            end do
        end do

        call update_charge_density()

        call update_electric_field()
        

    end subroutine initialize_variables

    subroutine calculate_reaction_profile()
        
        integer :: i, j, k
        double precision :: flame_height, distance_flame, distance_nearest, r_norm
        integer :: i_edge

        ! initial value of flame edge
        i_edge = 1

        ! solve coservation equation
        do i = 1, nr

            do j = 1, nz

                ! initialize distance as larger length
                distance_nearest = length_z

                ! find nearest distance to flame
                do k = 1, nr
                
                    ! get normalized r distance of referenced point
                    r_norm = distance_r(k, 1)/length_r

                    ! calculate flame height from fitting eqations in Logistic function
                    ! flame_height = length_z*(7.374e-01/(1 + exp(7.634e+00*(r_norm-8.732e-01))) - 1.542e-01)
                    flame_height = length_z*(-4.195 * exp(-(r_norm*25.0)**2/(2*4.054**2)) + 18.812)/25.0

                    ! calculated distance to flame height
                    distance_flame = sqrt((distance_z(i, j) - flame_height)**2.0 + (distance_r(i, j) - distance_r(k ,1))**2.0)

                    ! save narest distance
                    distance_nearest = min(distance_flame, distance_nearest)

                end do
                    
                ! get normalized r distance of current point
                r_norm = distance_r(i, 1)/length_r

                ! calculate flame height from fitting eqations in Logistic function
                ! flame_height = length_z*(7.374e-01/(1 + exp(7.634e+00*(r_norm-8.732e-01))) - 1.542e-01)
                flame_height = length_z*(-4.195 * exp(-(r_norm*25.0)**2/(2*4.054**2)) + 18.812)/25.0
                ! rounding to nearest point on z
                flame_height = floor(flame_height/dz)*dz
                ! calculate distance_flame on z direaction
                distance_flame = abs(distance_z(i, j) - flame_height)

    
                ! calclate spacial profile of ionization
                ! g_i(i, j) = exp(- (pi*distance_flame**2)/a_thickness**2)
                ! fitting to 1D PREMIX of Yuhia Ren
                g_i(i, j) = exp(- (distance_flame)**2/6.554209032d-09)

                ! calclate spacial profile of attachment
                ! fitting to 1D PREMIX of Yuhia Ren
                ! g_a(i, j) = (erf((flame_height - distance_z(i, j))/2.74084987d-04) + 1.0)/2.0
                g_a(i, j) = (erf((flame_height - distance_z(i, j))/1.82371556e-04) + 1.0)/2.0

                ! calclate temperature
                T(i, j) = 9.20603021e+02 * erf((distance_z(i, j) - flame_height) / 2.73035935e-04) + 1.21788749e+03
                ! T(i, j) = 2000.0 ! to make debye length larger

                ! calculate diffusion coefficients
                D_pos(i, j) = K_pos*k_B*T(i, j)/q_e ! diffusion coefficients of positive ions [m2/s]
                D_neg(i, j) = K_neg*k_B*T(i, j)/q_e ! diffusion coefficients of negative ions [m2/s]
                D_ele(i, j) = K_ele*k_B*T(i, j)/q_e ! diffusion coefficients of electrons [m2/s]

                ! ! calculate diffusion coefficients
                ! D_pos(i, j) = 5.0d-5 ! diffusion coefficients of positive ions [m2/s]
                ! D_neg(i, j) = 5.0d-5  ! diffusion coefficients of negative ions [m2/s]
                ! D_ele(i, j) = 6.89d-2 ! diffusion coefficients of electrons [m2/s]

                ! n_pos(i, j) = 1.0d10*g_i(i, j)
                ! n_neg(i, j) = 1.0d10*g_i(i, j)
                ! n_ele(i, j) = 1.0d10*g_i(i, j)

            end do
        end do

    end subroutine calculate_reaction_profile

    subroutine calculate_output_variables()

        integer :: i, j

        do i = 1, nr
            do j = 1, nz

                ! current density for r direction
                J_r(i, j) = n_pos(i, j)*q_e*K_pos*E_r(i, j) &
                          - n_neg(i, j)*q_e*K_neg*E_r(i, j) &
                          - n_ele(i, j)*q_e*K_ele*E_r(i, j)

                ! current density for z direction
                J_z(i, j) = n_pos(i, j)*q_e*K_pos*E_z(i, j) &
                          - n_neg(i, j)*q_e*K_neg*E_z(i, j) &
                          - n_ele(i, j)*q_e*K_ele*E_z(i, j)

                ! electric body force for r direction
                F_r = rho(i, j)*E_r(i, j)
                
                ! electric body force for z direction
                F_z = rho(i, j)*E_z(i, j)

            end do
        end do

    end subroutine calculate_output_variables

    subroutine update_charge_density()

        integer :: i, j

        ! calclate density of electric charge
        do i = 1, nr
            do j = 1, nz
                rho(i, j) = (n_pos(i, j) - n_neg(i, j) - n_ele(i, j))*q_e
            end do
        end do

    end subroutine update_charge_density

    subroutine update_electric_field()
        integer :: i, j

        ! calclate electric field (E = -dV/dz)
        do i = 2, nr-1
            do j = 2, nz-1
                E_r(i, j) = -(V(i+1, j) - V(i-1, j)) / (2.0d0 * dr)
                E_z(i, j) = -(V(i, j+1) - V(i, j-1)) / (2.0d0 * dz)
            end do
        end do

        ! boundary conditions on z axis
        do i = 1, nr

            ! lower boundary
            E_r(i, 1) = 0.0d0
            E_z(i, 1) = -(V(i, 2) - V(i, 1)) / dz
            
            ! upper boundary
            E_r(i, nz) = 0.0d0
            E_z(i, nz) = -(V(i, nz) - V(i, nz-1)) / dz

        end do

        ! boundary conditions on r axis
        do j = 2, nz-1

            ! lower boundary
            E_r(1, j) = -(V(2, j) - V(1, j)) / dr
            E_z(1, j) = -(V(1, j+1) - V(1, j-1)) / (2.0d0 * dz)
            
            ! upper boundary
            E_r(nr, j) = -(V(nr, j) - V(nr-1, j)) / dr
            E_z(nr, j) = -(V(nr, j+1) - V(nr, j-1)) / (2.0d0 * dz)

        end do

    end subroutine update_electric_field

    subroutine update_ddV()
    
        integer :: i, j
        double precision :: ddVdr_raw, ddVdz_raw, ddVdr_eff, ddVdz_eff

        do i = 2, nr-1
            do j = 2, nz-1

                ! 現在の電位 V から計算された生の勾配を ddVdz_raw とすると
                ddVdr_raw = (V(i+1, j) - 2.0*V(i, j) + V(i-1, j))/(dr**2)
                ddVdz_raw = (V(i, j+1) - 2.0*V(i, j) + V(i, j-1))/(dz**2)

                ! 下方緩和を適用して、密度計算に使う "実効的な" 勾配を求める
                ddVdr_eff = (1.0d0 - alpha_E) * ddVdr_prev(i, j) + alpha_E * ddVdr_raw
                ddVdz_eff = (1.0d0 - alpha_E) * ddVdz_prev(i, j) + alpha_E * ddVdz_raw

                ! 後のために prev を更新しておく
                ddVdr_prev(i, j) = ddVdr_eff
                ddVdz_prev(i, j) = ddVdz_eff

            end do
        end do

    end subroutine update_ddV

    subroutine update_source_term()

        integer :: i, j

        do i = 1, nr
            do j = 1, nz
                
                ! calcualte source term for positive ion
                Sp_pos(i, j) = - k_r*n_ele(i, j)
                Su_pos(i, j) = k_i*g_i(i, j)
                
                ! calcualte source term for negative ion
                Sp_neg(i, j) = 0.0
                Su_neg(i, j) = k_a*g_a(i, j)*n_ele(i, j)
                
                ! calcualte source term for electron
                Sp_ele(i, j) = - k_r*n_pos(i, j) - k_a*g_a(i, j)
                Su_ele(i, j) = k_i*g_i(i, j)

            end do
        end do

    end subroutine update_source_term

    subroutine export_variables()
        implicit none
        integer :: i, j

        ! Save 2D arrays as binary files for later Fortran reading (all in one file)
        open(unit=10, file='output/variables_all.bin', form='unformatted', access='stream', status='replace')
        write(10) V
        write(10) n_pos
        write(10) n_neg
        write(10) n_ele
        close(10)

        open(unit=1, file='output/grid_r.dat', status='replace')
        do j = nz, 1, -1
            write(1, *) (distance_r(i, j), i = 1, nr)
        end do
        close(1)

        open(unit=1, file='output/grid_z.dat', status='replace')
        do j = nz, 1, -1
            write(1, *) (distance_z(i, j), i = 1, nr)
        end do
        close(1)

        open(unit=1, file='output/charge_density.dat', status='replace')
        do j = nz, 1, -1
            write(1, *) (rho(i, j), i = 1, nr)
        end do
        close(1)

        open(unit=1, file='output/potential.dat', status='replace')
        do j = nz, 1, -1
            write(1, *) (V(i, j), i = 1, nr)
        end do
        close(1)
        
        open(unit=1, file='output/electric_field_r.dat', status='replace')
        do j = nz, 1, -1
            write(1, *) (E_r(i, j), i = 1, nr)
        end do
        close(1)
        
        open(unit=1, file='output/electric_field_z.dat', status='replace')
        do j = nz, 1, -1
            write(1, *) (E_z(i, j), i = 1, nr)
        end do
        close(1)

        open(unit=1, file='output/positive_ions.dat', status='replace')
        do j = nz, 1, -1
            write(1, *) (n_pos(i, j), i = 1, nr)
        end do
        close(1)
        
        open(unit=1, file='output/negative_ions.dat', status='replace')
        do j = nz, 1, -1
            write(1, *) (n_neg(i, j), i = 1, nr)
        end do
        close(1)

        open(unit=1, file='output/electrons.dat', status='replace')
        do j = nz, 1, -1
            write(1, *) (n_ele(i, j), i = 1, nr)
        end do
        close(1)
        
        open(unit=1, file='output/profile_ionization.dat', status='replace')
        do j = nz, 1, -1
            write(1, *) (g_i(i, j), i = 1, nr)
        end do
        close(1)
        
        open(unit=1, file='output/profile_attachment.dat', status='replace')
        do j = nz, 1, -1
            write(1, *) (g_a(i, j), i = 1, nr)
        end do
        close(1)

        open(unit=1, file='output/temperature.dat', status='replace')
        do j = nz, 1, -1
            write(1, *) (T(i, j), i = 1, nr)
        end do
        close(1)

        open(unit=1, file='output/current_density_r.dat', status='replace')
        do j = nz, 1, -1
            write(1, *) (J_r(i, j), i = 1, nr)
        end do
        close(1)
        
        open(unit=1, file='output/current_density_z.dat', status='replace')
        do j = nz, 1, -1
            write(1, *) (J_z(i, j), i = 1, nr)
        end do
        close(1)

        open(unit=1, file='output/electric_body_force_r.dat', status='replace')
        do j = nz, 1, -1
            write(1, *) (F_r(i, j), i = 1, nr)
        end do
        close(1)
        
        open(unit=1, file='output/electric_body_force_z.dat', status='replace')
        do j = nz, 1, -1
            write(1, *) (F_z(i, j), i = 1, nr)
        end do
        close(1)

    end subroutine export_variables

    subroutine import_variables(foldername)
        implicit none
        character(len=*), intent(in) :: foldername
        character(len=256) :: filepath

        ! Compose the file path using the provided folder name
        write(filepath, '(A,"/variables_all.bin")') trim(foldername)

        ! Read 2D arrays from binary file
        open(unit=10, file=filepath, form='unformatted', access='stream', status='old')
        read(10) V
        read(10) n_pos
        read(10) n_neg
        read(10) n_ele
        close(10)

    end subroutine import_variables

end module variables_module
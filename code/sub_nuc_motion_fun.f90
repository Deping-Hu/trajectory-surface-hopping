      subroutine nuc_motion_fun(current_state,n_state,x_cur,p_cur, &
                 n_mode,k_x_tmp,k_p_tmp,label_diabatic) 

      implicit none
      include "param.def"

      integer :: n_mode, i_mode, current_state, n_state, label_diabatic
      double precision :: kinetic
      double precision, dimension(n_state) :: exci_e
      double precision, dimension(n_mode) :: x_cur,p_cur, &
                                             k_x_tmp,k_p_tmp, &
                                             freq
      double precision, dimension(n_state,n_state) :: dia_hamiton, &
                                                      a_hamiton, &
                                                      lambda
      double precision, dimension(n_state,n_mode) :: k_tun
      double precision, dimension(n_mode,n_state,n_state) :: delta_v_adiabatic, &
                                                             delta_v_diabatic
      complex (kind=8), dimension(n_state,n_state) :: dtoa_u

      call gradient(n_state,n_mode,x_cur,p_cur,freq,k_tun,exci_e, &
                    lambda,kinetic,dia_hamiton,a_hamiton,dtoa_u, &
                    delta_v_adiabatic,delta_v_diabatic) 

      do i_mode = 1,n_mode

         k_x_tmp(i_mode) = p_cur(i_mode)*freq(i_mode)

         if (label_diabatic .eq. 1) then
            k_p_tmp(i_mode) = &
            -delta_v_diabatic(i_mode,current_state,current_state)

         else if (label_diabatic .eq. 2) then
            k_p_tmp(i_mode) = & 
            -delta_v_adiabatic(i_mode,current_state,current_state)

         else if (label_diabatic .eq. 3) then
            k_p_tmp(i_mode) = & 
            -delta_v_adiabatic(i_mode,current_state,current_state)

         endif

      enddo

      end subroutine nuc_motion_fun

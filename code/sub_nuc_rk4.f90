      subroutine nuc_rk4(nuc_dt,current_state,n_state, &
                 x_cur,p_cur,n_mode,label_diabatic) 

      implicit none
      include "param.def"

      integer :: n_mode, current_state, n_state, label_diabatic
      double precision :: nuc_dt
      double precision, dimension(n_mode) :: x_cur,p_cur, &
                                             x_cur_tmp,p_cur_tmp, &
                                             k_x_tmp,k_p_tmp, &
                                             k1_x, k1_p, &
                                             k2_x, k2_p, &
                                             k3_x, k3_p, &
                                             k4_x, k4_p

      x_cur_tmp = x_cur
      p_cur_tmp = p_cur
 

      call nuc_motion_fun(current_state,n_state,x_cur,p_cur, &
                          n_mode,k_x_tmp,k_p_tmp,label_diabatic)

      k1_x = k_x_tmp
      k1_p = k_p_tmp

       
      x_cur = x_cur_tmp + 0.5d0*nuc_dt*k1_x
      p_cur = p_cur_tmp + 0.5d0*nuc_dt*k1_p


      call nuc_motion_fun(current_state,n_state,x_cur,p_cur, &
                          n_mode,k_x_tmp,k_p_tmp,label_diabatic)
      k2_x = k_x_tmp
      k2_p = k_p_tmp

      x_cur = x_cur_tmp + 0.5d0*nuc_dt*k2_x
      p_cur = p_cur_tmp + 0.5d0*nuc_dt*k2_p
 

      call nuc_motion_fun(current_state,n_state,x_cur,p_cur, &
                          n_mode,k_x_tmp,k_p_tmp,label_diabatic)
      k3_x = k_x_tmp
      k3_p = k_p_tmp

      x_cur = x_cur_tmp + nuc_dt*k3_x
      p_cur = p_cur_tmp + nuc_dt*k3_p


      call nuc_motion_fun(current_state,n_state,x_cur,p_cur, &
                          n_mode,k_x_tmp,k_p_tmp,label_diabatic)
      k4_x = k_x_tmp
      k4_p = k_p_tmp

      x_cur = x_cur_tmp + &
              (1.0d0/6.0d0)*nuc_dt*(k1_x + &
              2.0d0*k2_x+2.0d0*k3_x+k4_x)
      p_cur = p_cur_tmp + &
              (1.0d0/6.0d0)*nuc_dt*(k1_p+ &
              2.0d0*k2_p+2.0d0*k3_p+k4_p)



      end subroutine nuc_rk4

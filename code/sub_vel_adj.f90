      subroutine vel_adj(n_mode,n_state,current_state, &
                         current_state_pre,x_cur,p_cur,label_diabatic)

      implicit none
      include "param.def"

      integer :: i,j,k,n_mode,n_state,current_state, &
                 current_state_pre,label_diabatic
      double precision :: adj_r,nac_vec_total,kinetic
      double precision, dimension(n_mode) :: x_cur,p_cur,freq
      double precision, dimension(n_state) :: exci_e
      double precision, dimension(n_state,n_mode) :: k_tun
      double precision, dimension(n_state,n_state) :: adj_a,adj_b, &
                                                      delta,lambda
      double precision, dimension(n_state,n_state) :: dia_hamiton, &
                                                      a_hamiton, &
                                                      nac
      double precision, dimension(n_mode,n_state,n_state) :: delta_v_adiabatic, &
                                                             delta_v_diabatic, & 
                                                             nac_vec, &
                                                             adj_w
      complex (kind=8), dimension(n_state,n_state) :: dtoa_u


      call coupling(n_state,n_mode,x_cur,p_cur,freq,k_tun, &
                    exci_e,lambda,kinetic,dia_hamiton,a_hamiton,dtoa_u, &
                    delta_v_adiabatic,delta_v_diabatic,nac,nac_vec) 
     
      i = current_state_pre
      j = current_state

      adj_a = 0.d0
      adj_b = 0.d0
      delta = 0.d0
      adj_w = sqrt(1.d0 / n_mode)

      if (label_diabatic .eq. 2) then
         nac_vec_total = 0.d0
         do k=1,n_mode
            nac_vec_total = nac_vec_total + nac_vec(k,i,j)**2 
         enddo
         do k=1,n_mode
            adj_w(k,i,j) = nac_vec(k,i,j) / sqrt(nac_vec_total)
         enddo
      endif


      do k=1,n_mode
         adj_a(i,j) = adj_a(i,j) + 0.5d0 * &
                      adj_w(k,i,j)**2*freq(k)
         adj_b(i,j) = adj_b(i,j) + p_cur(k)*adj_w(k,i,j)*freq(k)
      enddo
      

      if (label_diabatic .eq. 1) then
         delta(i,j) = adj_b(i,j)**2 + 4*adj_a(i,j)* &
             (dia_hamiton(i,i)-dia_hamiton(j,j))
      else if (label_diabatic .eq. 2) then
         delta(i,j) = adj_b(i,j)**2 + 4*adj_a(i,j)* &
             (a_hamiton(i,i)-a_hamiton(j,j))
      else if (label_diabatic .eq. 3) then
         delta(i,j) = adj_b(i,j)**2 + 4*adj_a(i,j)* &
             (a_hamiton(i,i)-a_hamiton(j,j))
      endif

     
      if (delta(i,j) .lt. 0.d0) then
         current_state = current_state_pre
      else
          if (adj_b(i,j) .lt. 0.d0) then
             adj_r = (adj_b(i,j) + sqrt(delta(i,j))) / &
                     (2*adj_a(i,j))
          else
             adj_r = (adj_b(i,j) - sqrt(delta(i,j))) / &
                     (2*adj_a(i,j))
          endif
           
          do k=1,n_mode
             p_cur(k) = p_cur(k) - adj_r * adj_w(k,i,j)
          enddo

      endif
      
      end subroutine vel_adj

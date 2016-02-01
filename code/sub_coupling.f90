      subroutine coupling(n_state,n_mode,x_cur,p_cur,freq,k_tun,exci_e, &
                          lambda,kinetic,dia_hamiton,a_hamiton,dtoa_u, &
                          delta_v_adiabatic,delta_v_diabatic,nac,nac_vec)

      implicit none

      include "param.def"

      integer :: i_mode,n_mode, n_state, i, j 
      double precision :: kinetic
      double precision, dimension(n_mode) :: x_cur,p_cur,freq
      double precision, dimension(n_state) :: exci_e
      double precision, dimension(n_state,n_state) :: a_hamiton, &
                                                      dia_hamiton, &
                                                      lambda, &
                                                      nac
      double precision, dimension(n_state,n_mode) :: k_tun
      double precision, dimension(n_mode,n_state,n_state) :: delta_v_adiabatic, &
                                                             delta_v_diabatic, & 
                                                             nac_vec

      complex (kind=8), dimension(n_state,n_state) :: dtoa_u

      call gradient(n_state,n_mode,x_cur,p_cur,freq,k_tun,exci_e, &
                    lambda,kinetic,dia_hamiton,a_hamiton,dtoa_u, &
                    delta_v_adiabatic,delta_v_diabatic) 

      do i_mode=1,n_mode
         nac = 0.d0
         do i=1,n_state
            do j=1,n_state
               if (i .ne. j) then
                  nac_vec(i_mode,i,j) = delta_v_adiabatic(i_mode,i,j) &
                  / (a_hamiton(j,j) - a_hamiton(i,i))
               else
                  nac_vec(i_mode,i,j) = 0.d0
               endif 
               
               nac(i,j) = nac(i,j) + nac_vec(i_mode,i,j) * &
                          freq(i_mode)*p_cur(i_mode)
            enddo
         enddo
      enddo

      end subroutine coupling

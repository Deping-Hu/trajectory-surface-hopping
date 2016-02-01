      subroutine gradient(n_state,n_mode,x_cur,p_cur,freq,k_tun,exci_e, &
                          lambda,kinetic,dia_hamiton,a_hamiton,dtoa_u, &
                          delta_v_adiabatic,delta_v_diabatic)

      implicit none

      include "param.def"

      integer :: i_mode,n_mode, n_state, i, j 
      double precision :: kinetic
      double precision, dimension(n_mode) :: x_cur,p_cur,freq
      double precision, dimension(n_state) :: exci_e
      double precision, dimension(n_state,n_mode) :: k_tun
      double precision, dimension(n_state,n_state) :: a_hamiton, &
                                                      dia_hamiton, &
                                                      lambda, &
                                                      delta_v_tmp2, &
                                                      delta_v_tmp1
      double precision, dimension(n_mode,n_state,n_state) :: delta_v_adiabatic, &
                                                             delta_v_diabatic
      complex (kind=8), dimension(n_state,n_state) :: dtoa_u
                                                     
      call read_hamiton(n_state,n_mode,x_cur,p_cur,freq,k_tun,exci_e, &
                        lambda,kinetic,dia_hamiton,a_hamiton,dtoa_u)

!diabatic gradient
      do i_mode = 1,n_mode
         do i=1,n_state
            do j=1,n_state 
               if (i .eq. j) then
                  delta_v_diabatic(i_mode,i,j) = (freq(i_mode)*x_cur(i_mode)+ &
                                     k_tun(i,i_mode))
               else if (i .ne. j) then
                  if (i_mode .eq. 3) then
                     delta_v_diabatic(i_mode,i,j) = lambda(i,j) 
                  else
                     delta_v_diabatic(i_mode,i,j) = 0.d0
                  endif
               endif
            enddo
         enddo
      enddo 

!adiabatic gradient

      do i_mode=1,n_mode

         do i=1,n_state
           do j=1,n_state
              delta_v_tmp1(i,j) = delta_v_diabatic(i_mode,i,j)
           enddo
         enddo

         delta_v_tmp2 = matmul(dtoa_u,delta_v_tmp1)
         delta_v_tmp2 = matmul(delta_v_tmp2,conjg(dtoa_u))
         
         do i=1,n_state
           do j=1,n_state
              delta_v_adiabatic(i_mode,i,j) = delta_v_tmp2(i,j)
           enddo
         enddo

      enddo

      end subroutine gradient

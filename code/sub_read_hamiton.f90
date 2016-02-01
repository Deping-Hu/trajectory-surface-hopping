      subroutine read_hamiton(n_state,n_mode,x_cur,p_cur,freq,k_tun, &
                 exci_e,lambda,kinetic,dia_hamiton,a_hamiton,dtoa_u)

      implicit none

      include "param.def"

      integer :: n_mode, n_state, i, j, k, lwork, info
      double precision :: kinetic
      double precision, dimension(n_mode) :: x_cur,p_cur,freq
      double precision, dimension(n_state) :: exci_e, e_adiabatic
      double precision, dimension(n_state,n_mode) :: k_tun
      double precision, dimension(n_state,n_state) :: dia_hamiton, &
                                                      a_hamiton, &
                                                      lambda
      complex (kind=8), dimension(n_state,n_state) :: dtoa_u

      double precision, allocatable, dimension(:) :: work
      complex (kind=8), allocatable, dimension(:,:) :: rwork
                                                     
      lwork= n_state*(n_state+1)

      allocate (work(lwork))
      allocate (rwork(3*n_state-2,3*n_state-2))

      call read_parameter(n_state,n_mode,freq,k_tun,exci_e,lambda)
      
      kinetic = 0.d0

      do i=1,n_mode
         kinetic = kinetic + 0.5d0*freq(i)*p_cur(i)**2
      enddo

      do i =1,n_state
         do j=1,n_state
            if (i .eq. j) then
               dia_hamiton(i,j) = exci_e(i)
               do k=1,n_mode
                  dia_hamiton(i,j)= dia_hamiton(i,j) &
                   + k_tun(i,k)*x_cur(k) &
                   + 0.5d0* freq(k)*x_cur(k)**2
               enddo
            else if (i .ne. j) then
               dia_hamiton(i,j)=lambda(i,j)*x_cur(3)
            endif
         enddo
      enddo


      dtoa_u = dia_hamiton
      call  ZHEEV('V','U', n_state, dtoa_u, n_state, &
                  e_adiabatic, work, lwork, rwork, info)

      a_hamiton = 0.d0

      do i=1,n_state
         a_hamiton(i,i) = e_adiabatic(i)
      enddo

      end subroutine read_hamiton

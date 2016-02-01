      subroutine read_parameter(n_state,n_mode,freq,k_tun,exci_e,lambda) 

      implicit none

      include "param.def"

      integer :: n_state,i_mode, n_mode, file_freq_input

      double precision, dimension(n_mode) :: freq
      double precision, dimension(n_state) :: exci_e
      double precision, dimension(n_state,n_state) :: lambda
      double precision, dimension(n_state,n_mode) :: k_tun


      file_freq_input=11

      open(unit=file_freq_input, file="freq.input")

      do i_mode = 1 , n_mode
         read (file_freq_input,*) freq(i_mode)
         freq(i_mode)=freq(i_mode)/TOEV
      enddo

      close(11)

      exci_e(1) = 3.94/TOEV
      exci_e(2) = 4.84/TOEV
      k_tun(1,1)=0.037/TOEV
      k_tun(1,2)=-0.254/TOEV
      k_tun(1,3)=0.d0
      k_tun(2,1)=-0.105/TOEV
      k_tun(2,2)=0.149/TOEV
      k_tun(2,3)=0.d0
      lambda(1,2)=0.262/TOEV
      lambda(2,1)=0.262/TOEV
 

      if (n_state .eq. 3) then

      exci_e(3) = 5.54/TOEV
      k_tun(3,1)=0.20/TOEV
      k_tun(3,2)=-0.10/TOEV
      k_tun(3,3)=0.d0
      lambda(1,3)=0.200/TOEV
      lambda(3,1)=0.200/TOEV
      lambda(2,3)=0.300/TOEV
      lambda(3,2)=0.300/TOEV

      endif
 
      end subroutine read_parameter

      subroutine read_ini_trj(label_restart,n_mode,n_state,x_cur,p_cur, &
                 ele_coe) 

      implicit none

      include "param.def"

      integer :: i_mode, n_mode, n_state,label_restart, &
                 file_ini_input, file_restart_input

      double precision :: q
      double precision, dimension(n_mode) :: x_cur,p_cur
      complex (kind=8), dimension(n_state) :: ele_coe

      call RANDOM_NUMBER(q)

      q=q*PI

      ele_coe(1)=cmplx(0.d0,0.d0)
      ele_coe(2)=cmplx(cos(q),sin(q))
!      ele_coe(3)=(0.d0,0.d0)

      if (label_restart .eq. 0 ) then

         file_ini_input=11
         open(unit=file_ini_input, file="trj.input")
         do i_mode=1, n_mode
            read (file_ini_input, *)  x_cur(i_mode), p_cur(i_mode)
         enddo
         close(11)

      endif


      if (label_restart .eq. 1 ) then

         file_restart_input=12
         open(unit=file_restart_input, file="trj_restart.input")
         do i_mode=1, n_mode
            read (file_restart_input, *)  x_cur(i_mode), p_cur(i_mode)
         enddo
         close(12)

      endif

      end subroutine read_ini_trj

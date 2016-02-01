      subroutine hop(i_step,n_state,current_state,current_state_pre,hop_pro) 

      implicit none
      include "param.def"

      integer :: n_state,i,j,i_step,current_state, &
                 current_state_pre,file_random_number
      double precision :: hop_tmp
      double precision, dimension(n_state,n_state) :: hop_pro,random_hop

      call RANDOM_NUMBER(random_hop)

      file_random_number = 11

      if (i_step .eq. 0) then
         open(unit=11,file="random.dat")
      else
         open(unit=11,position='append',file="random.dat")
      endif

      write(11,9990) i_step,random_hop(1,2),random_hop(2,1)

      close(11)

9990  format(i8,1x,10(f20.10, 1x))

      current_state_pre = current_state 

      i = current_state
      hop_tmp = 0.d0
      do j=1,n_state
         if (j .ne. i) then
            hop_tmp = hop_tmp + hop_pro(i,j)
            if (random_hop(i,j) .le. hop_tmp) then
               current_state = j
               exit
            endif
         endif
      enddo

      end subroutine hop

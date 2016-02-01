      subroutine write_trj(i_step,n_state,current_state,n_mode,&
                 x_cur,p_cur,ele_coe,label_diabatic) 

      implicit none

      include "param.def"

      integer :: i_step, i_mode, n_mode,i_state,n_state, current_state, &
                 file_cur_trj, file_restart_trj, file_cur_nac, &
                 file_cur_state,file_cur_energy_a,file_cur_energy_d, &
                 file_cur_population_d,file_cur_population_a,label_diabatic

      double precision :: kinetic,random_pop,pop_excited_state
      double precision, dimension(n_mode) :: x_cur,p_cur,freq
      double precision, dimension(n_state) :: pop_d,pop_a,exci_e
      double precision, dimension(n_state,n_state) :: dia_hamiton, &
                                                      a_hamiton, &
                                                      nac, &
                                                      lambda
      double precision, dimension(n_state,n_mode) :: k_tun
      double precision, dimension(n_mode,n_state,n_state) :: delta_v_adiabatic, &
                                                             delta_v_diabatic, & 
                                                             nac_vec
      complex (kind=8), dimension (n_state) :: ele_coe,ele_coe_a,ele_coe_d
      complex (kind=8), dimension (n_state,n_state) :: dtoa_u


! energy

      file_cur_energy_d=10
      file_cur_energy_a=11

      call coupling(n_state,n_mode,x_cur,p_cur,freq,k_tun, &
                    exci_e,lambda,kinetic,dia_hamiton,a_hamiton,dtoa_u, &
                    delta_v_adiabatic,delta_v_diabatic,nac,nac_vec) 

      
      if (i_step .eq. 0) then
         open(unit=file_cur_energy_d,file="energy_d.dat")
         open(unit=file_cur_energy_a,file="energy_a.dat")
      else
         open(unit=file_cur_energy_d, position='append',file="energy_d.dat")
         open(unit=file_cur_energy_a, position='append',file="energy_a.dat")
      endif

      write (file_cur_energy_d, 9999,advance='no') i_step, kinetic
      write (file_cur_energy_a, 9999,advance='no') i_step, kinetic

      do i_state=1,n_state
          write (file_cur_energy_d, 9998,advance='no')  & 
          dia_hamiton(i_state,i_state)
          write (file_cur_energy_a, 9998,advance='no')  &
          a_hamiton(i_state,i_state)
      enddo

      write (file_cur_energy_d, 9998) &
      kinetic + dia_hamiton(current_state,current_state)

      write (file_cur_energy_a, 9998) &
      kinetic + a_hamiton(current_state,current_state)

      close(10)
      close(11)


! population
      if (i_step .eq. 0) then
         if (label_diabatic .eq. 2) then
            ele_coe = matmul(dtoa_u,ele_coe)
            pop_excited_state=conjg(ele_coe(2))*ele_coe(2)
            call RANDOM_NUMBER(random_pop)
            if (pop_excited_state .lt. random_pop) then
               current_state = 1
            else
               current_state = 2
            endif 
         endif
      endif

      file_cur_population_d=13
      file_cur_population_a=14

      if (label_diabatic .eq. 1) then
         ele_coe_d = ele_coe
         ele_coe_a = matmul(dtoa_u,ele_coe_d)
      else if (label_diabatic .eq. 2) then
         ele_coe_a = ele_coe
         ele_coe_d = matmul(ele_coe_a,conjg(dtoa_u))
      else if (label_diabatic .eq. 3) then
         ele_coe_d = ele_coe
         ele_coe_a = matmul(dtoa_u,ele_coe_d)
      endif

      do i_state=1,n_state
         pop_d(i_state) = real(conjg(ele_coe_d(i_state))*ele_coe_d(i_state))
         pop_a(i_state) = real(conjg(ele_coe_a(i_state))*ele_coe_a(i_state))
      enddo

      if (i_step .eq. 0) then
      open(unit=file_cur_population_d,file="pop_d.dat")
      open(unit=file_cur_population_a,file="pop_a.dat")
      else
      open(unit=file_cur_population_d, position='append',file="pop_d.dat")
      open(unit=file_cur_population_a, position='append',file="pop_a.dat")
      endif
      
      write (file_cur_population_d, 9999) i_step, pop_d
      write (file_cur_population_a, 9999) i_step, pop_a

      close(13)
      close(14)

! state
      file_cur_state=12
      if (i_step .eq. 0) then
         open(unit=file_cur_state,file="state.dat")
      else
         open(unit=file_cur_state, position='append',file="state.dat")
      endif

      write (file_cur_state, *) i_step, current_state

      close(12)


! trajectory

      file_restart_trj=15
      file_cur_trj=16
      if (i_step .eq. 0) then
      open(unit=file_cur_trj,file="trj.dat")
      else
      open(unit=file_cur_trj, position='append',file="trj.dat")
      endif

      open(unit=file_restart_trj,file="trj_restart.input")

      do i_mode=1, n_mode
         write (file_cur_trj, 9999) i_step, x_cur(i_mode), p_cur(i_mode)
         write (file_restart_trj, 9998) x_cur(i_mode), p_cur(i_mode)
      enddo
      
      close(15)
      close(16)

! coupling 
      file_cur_nac=17
      if (i_step .eq. 0) then
         open(unit=file_cur_nac,file="nac.dat")
      else
         open(unit=file_cur_nac,position='append',file="nac.dat")
      endif

      write(file_cur_nac, 9999) i_step, nac(1,2), nac(2,1)

      close(17)
      
     
9998  format(10(f20.10, 1x))
9999  format(i8,1x,10(f20.10, 1x))

      end subroutine write_trj

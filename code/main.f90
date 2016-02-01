      Program main

      implicit none

      include 'param.def'

      integer :: n_state, current_state, &
                 current_state_pre,n_mode, &
                 n_step, i_step, &
                 n_save_trj, &
                 label_restart, &
                 label_diabatic, &
                 label_hop_pro, &
                 file_dyn_in, &
                 n_ele_dt

      double precision :: nuc_dt

      double precision, allocatable, dimension(:) :: x_cur, p_cur, &
                                                     x_pre, p_pre

      double precision, allocatable, dimension(:,:) :: hop_pro

      complex (kind=8), allocatable, dimension(:) :: ele_coe


      namelist /dyn_control/ n_mode, n_step, nuc_dt, n_ele_dt, &
               n_save_trj,label_restart,label_diabatic, label_hop_pro,&
               n_state, current_state

      file_dyn_in=11       
      open(unit=file_dyn_in, file="dyn.input")      
      read(file_dyn_in, nml = dyn_control)
      close(11)
      

      allocate (x_cur(n_mode))
      allocate (p_cur(n_mode))
      allocate (x_pre(n_mode))
      allocate (p_pre(n_mode))
      allocate (hop_pro(n_state,n_state))
      allocate (ele_coe(n_state))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      nuc_dt = nuc_dt / TOFS
      call sub_init_random_seed()
      call read_ini_trj(label_restart,n_mode,n_state,x_cur,p_cur,ele_coe)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do i_step=0, n_step

   
         if (mod(i_step,n_save_trj) .eq. 0 ) then
            call write_trj(i_step,n_state,current_state,n_mode,x_cur, &
                 p_cur,ele_coe,label_diabatic)
         endif
   
   
         x_pre = x_cur
         p_pre = p_cur
   
         call nuc_rk4(nuc_dt,current_state,n_state, &
              x_cur,p_cur,n_mode,label_diabatic) 
   
         call ele_rk4(i_step,n_mode,n_state,&
              x_cur,p_cur,x_pre,p_pre,ele_coe,nuc_dt,n_ele_dt,&
              hop_pro,label_diabatic,label_hop_pro)
 
         call hop(i_step,n_state,current_state, &
              current_state_pre,hop_pro)
   
        if (current_state .ne. current_state_pre) then
            call vel_adj(n_mode,n_state,current_state, &
                 current_state_pre,x_cur,p_cur,label_diabatic)
        endif
   
      enddo

      deallocate (x_cur)
      deallocate (p_cur)
      deallocate (x_pre)
      deallocate (p_pre)
      deallocate (ele_coe)
      deallocate (hop_pro)
      

      end Program main

      subroutine ele_rk4(i_step,n_mode,n_state, x_cur,p_cur,x_pre,p_pre, &
                 ele_coe,nuc_dt,n_ele_dt,hop_pro,label_diabatic,label_hop_pro)

      implicit none

      include "param.def"

      integer :: i_step,n_mode, n_state, n_ele_dt, i_ele_step, i, j, &
                 label_diabatic, label_hop_pro, file_cur_hop_pro

      double precision, dimension(n_state) :: exci_e

      double precision :: ele_dt, nuc_dt, tmp_hop, kinetic

      double precision, dimension(n_mode) :: x_cur,p_cur,x_pre,p_pre,&
                                             x_cur_tmp,p_cur_tmp,freq

      double precision, dimension(n_state,n_state) :: dia_hamiton, &
                                                      a_hamiton, &
                                                      hop_pro, &
                                                      nac, &
                                                      lambda

      double precision, dimension(n_state,n_mode) :: k_tun

      double precision, dimension(n_mode,n_state,n_state) :: delta_v_adiabatic, &
                                                             delta_v_diabatic, & 
                                                             nac_vec

      complex (kind=8), dimension(n_state) :: ele_coe,k_coe_tmp, &
                                              ele_coe_pre,ele_coe_tmp, &
                                              k1_coe,k2_coe, &
                                              k3_coe,k4_coe

      complex (kind=8), dimension(n_state,n_state) :: dtoa_u

      hop_pro = 0.d0

      ele_coe_pre = ele_coe

      ele_dt = nuc_dt / n_ele_dt

      do i_ele_step =0, n_ele_dt

         x_cur_tmp = (x_cur - x_pre) / n_ele_dt * i_ele_step + x_pre
         p_cur_tmp = (p_cur - p_pre) / n_ele_dt * i_ele_step + p_pre

         call coupling(n_state,n_mode,x_cur_tmp,p_cur_tmp,freq,k_tun, &
                      exci_e,lambda,kinetic,dia_hamiton,a_hamiton,dtoa_u, &
                      delta_v_adiabatic,delta_v_diabatic,nac,nac_vec) 

         ele_coe_tmp = ele_coe
   
         call ele_motion_fun(n_state,label_diabatic,dia_hamiton, &
                             a_hamiton,nac,ele_coe,k_coe_tmp)
         k1_coe = k_coe_tmp      
         ele_coe = ele_coe_tmp + 0.5d0*ele_dt*k1_coe
   
         call ele_motion_fun(n_state,label_diabatic,dia_hamiton, &
                             a_hamiton,nac,ele_coe,k_coe_tmp)
         k2_coe = k_coe_tmp      
         ele_coe = ele_coe_tmp + 0.5d0*ele_dt*k2_coe
   
         call ele_motion_fun(n_state,label_diabatic,dia_hamiton, &
                             a_hamiton,nac,ele_coe,k_coe_tmp)
         k3_coe = k_coe_tmp      
         ele_coe = ele_coe_tmp + ele_dt*k3_coe
    
         call ele_motion_fun(n_state,label_diabatic,dia_hamiton, &
                             a_hamiton,nac,ele_coe,k_coe_tmp)
         k4_coe = k_coe_tmp      
   
         ele_coe = ele_coe_tmp + &
                   (1.0d0/6.0d0)*ele_dt*(k1_coe+ &
                   2.0d0*k2_coe+2.0d0*k3_coe+k4_coe)
         
   
         do i=1,n_state
            do j=1,n_state
               if (j .ne. i) then

                  if (label_diabatic .eq. 1) then
                     tmp_hop = -2 * ele_dt * &
                     aimag(conjg(ele_coe(i))*ele_coe(j)*dia_hamiton(i,j))
                  else if (label_diabatic .eq. 2) then
                     tmp_hop = 2 * ele_dt * &
                     real(conjg(ele_coe(i))*ele_coe(j)*nac(i,j)) 
                  else if (label_diabatic .eq. 3) then
                     tmp_hop = 2 * ele_dt * &
                     real(conjg(ele_coe(i))*ele_coe(j)*nac(i,j)) 
                  endif

                  if (label_hop_pro .eq. 1) then
                     hop_pro(i,j) = hop_pro(i,j) + &
                     tmp_hop / (conjg(ele_coe(i))*ele_coe(i))
                  else if (label_hop_pro .eq. 2) then
                     hop_pro(i,j) = hop_pro(i,j) + tmp_hop 
                  endif

               endif
            enddo  
         enddo
         
      enddo 


      do i=1,n_state
         do j=1,n_state
            if (j .ne. i) then
              if (label_hop_pro .eq. 2) then
                 hop_pro(i,j) = hop_pro(i,j) / &
                 (conjg(ele_coe_pre(i))*ele_coe_pre(i)) 
              else if(label_hop_pro .eq. 3) then
                 hop_pro(i,j) = ((conjg(ele_coe_pre(i))*ele_coe_pre(i)) &
                 - (conjg(ele_coe(i))*ele_coe(i))) &
                 / (conjg(ele_coe_pre(i))*ele_coe_pre(i)) &
                 * conjg(ele_coe_pre(j))*ele_coe_pre(j) &
                 / (1-conjg(ele_coe_pre(i))*ele_coe_pre(i))   
              endif
            endif
         enddo
      enddo   

      file_cur_hop_pro = 11
      if (i_step .eq. 0) then
         open(11,file="hop_pro.dat")
      else
         open(11,position='append',file="hop_pro.dat")
      endif

      write(11, 9997) i_step, hop_pro(1,2), hop_pro(2,1)

      close(11)

9997  format(i8,1x,10(f20.10, 1x))

      end subroutine ele_rk4

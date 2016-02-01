       subroutine ele_motion_fun(n_state,label_diabatic,dia_hamiton, &
                                 a_hamiton,nac,ele_coe,k_coe_tmp) 

       implicit none
       include "param.def"

       integer :: n_state,label_diabatic
       complex (kind=8), dimension(n_state) :: ele_coe, k_coe_tmp
       double precision, dimension(n_state,n_state) :: a_hamiton, &
                                                       dia_hamiton,nac

       if (label_diabatic .eq. 1) then
          k_coe_tmp = -c1 * matmul(dia_hamiton,ele_coe)
       else if (label_diabatic .eq. 2) then
          k_coe_tmp = -c1 * matmul(a_hamiton,ele_coe) &
                      - matmul(nac,ele_coe)
       else if (label_diabatic .eq. 3) then
          k_coe_tmp = -c1 * matmul(dia_hamiton,ele_coe)
       endif
       
       end subroutine ele_motion_fun

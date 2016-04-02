      subroutine find_ss_large_dd(ss_m)

* [ P A R A M E T E R S]
* ......................
      integer, parameter :: nss    = 24

* [ C O M M O N ]
* ...............
      common/wblock4/  den_m(nss), den_im(nss), gdot(nss), nnns, nnne

* [ V A R I A B L E S ]
* .....................      
      integer :: ss_m, temp_ss, i, j
      real    :: den_t(nss), temp_dd


      temp_dd = 0.0
      den_t = den_m + den_im
    
      do i = 1, nss
         if (den_t(i) > temp_dd) then
            temp_dd = den_t(i) 
            temp_ss = i
         endif
      end do

      ss_m = temp_ss

      return
      end 
    

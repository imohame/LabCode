      subroutine setary (array,value,npts)
c     implicit double precision (a-h,o-z)                                    dp
      dimension array(*)
      do 10 i=1,npts
   10 array(i)=value
      return
      end

      subroutine setse2(prop,c)
c     implicit double precision (a-h,o-z)                                    dp
c
c     initialize orthotropic material constants
c
      common/bk14/lfna(15),lfnt(6)
      dimension prop(*),c(4,*)
      check=prop(1)*prop(2)*prop(3)*prop(7)
      if (check.gt.1.0e-20) go to 10
      write(lfnt(2),30)
      call bye (2)
   10 c(1,1)=1.0/prop(1)
      c(2,2)=1.0/prop(2)
      c(3,3)=1.0/prop(3)
      c(4,4)=1.0/prop(7)
      c(1,2)=-prop(4)*c(2,2)
      c(1,3)=-prop(5)*c(3,3)
      c(2,3)=-prop(6)*c(3,3)
      c(1,4)=0.0
      c(2,4)=0.0
      c(3,4)=0.0
      do 20 i=1,4
      do 20 j=i,4
   20 c(j,i)=c(i,j)
      call inver (c)
      return
c
   30 format(/' orthotropic material properties are unacceptable ')
      end

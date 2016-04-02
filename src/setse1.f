      subroutine setse1(prop)
c     implicit double precision (a-h,o-z)                                    dp
c
c     initialize elastic material constants
c
      dimension prop(4,*)
      q1=prop(1,1)*prop(2,1)/((1.+prop(2,1))*(1.-2.0*prop(2,1)))
      q2=prop(1,1)*0.5/(1.+prop(2,1))
      q3=q1+2.0*q2
      prop(1,5)=prop(1,1)
      prop(2,5)=prop(2,1)
      do 10 i=1,4
      do 10 j=1,4
   10 prop(i,j)=0.0
      prop(1,1)=q3
      prop(2,2)=q3
      prop(3,3)=q3
      prop(4,4)=q2
      prop(1,2)=q1
      prop(1,3)=q1
      prop(2,1)=q1
      prop(2,3)=q1
      prop(3,1)=q1
      prop(3,2)=q1
ck      write(25,20) q1,q2,q3
ck   20 format(5x,'q1=',e14.4,2x,'q2=',e14.4,2x,'q3=',e14.4/)
      return
      end

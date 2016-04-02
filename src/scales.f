      subroutine scales(s,ne,matp,time,eldeth,elbury)
c     implicit double precision (a-h,o-z)                                    dp
c
c     scale element out of solution leaving trace to prevent
c     instabilities
c
       use mod_parameters
      common/double/iprec,ncpw,unit
      common/vect1/r1(nelemg),r2(nelemg),r3(nelemg),
     > r4(nelemg),r5(nelemg),r6(nelemg),
     1 r7(nelemg),r8(nelemg),mtype(nelemg),mte(nelemg)
      dimension s(*)
c
      scalef=1.0-(time-eldeth)/(elbury-eldeth)
      if(scalef.lt.0.) matp=0
      scalef=max(1.e-8*unit,scalef)
c
      r1(ne)=scalef*r1(ne)
      r2(ne)=scalef*r2(ne)
      r3(ne)=scalef*r3(ne)
      r4(ne)=scalef*r4(ne)
      r5(ne)=scalef*r5(ne)
      r6(ne)=scalef*r6(ne)
      r7(ne)=scalef*r7(ne)
      r8(ne)=scalef*r8(ne)
c
      do 10 n=1,36
   10 s(n)=scalef*s(n)
c
      return
      end

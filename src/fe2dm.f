      subroutine fe2dm(matype,amas,ix,matp,ym,id)
c     implicit double precision (a-h,o-z)                                    dp
c
      common/bk02/ioofc,iphase,imass,lpar(9)
      dimension lm(8),matype(*),amas(*),ix(4,*),matp(*),ym(4,*),id(2,*)
c
      equivalence (lpar(2),numel)
c

        write(7777,*) '-- fe2dm.f'

c
      do 10 n=1,numel
      mtype=matp(n)
      if (mtype.eq.0) go to 10
c
      lm(1)=id(1,ix(1,n))
      lm(2)=id(2,ix(1,n))
      lm(3)=id(1,ix(2,n))
      lm(4)=id(2,ix(2,n))
      lm(5)=id(1,ix(3,n))
      lm(6)=id(2,ix(3,n))
      lm(7)=id(1,ix(4,n))
      lm(8)=id(2,ix(4,n))
c
      call addmas (lm,ym(1,n),amas)
c
   10 continue
      return
      end

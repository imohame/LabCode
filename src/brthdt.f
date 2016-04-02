      subroutine brthdt(freep,mtype,beta,matp,lmm,s,iprec)
c     implicit double precision (a-h,o-z)                                    dp
c
c     element birth and death
c
      common/bk32/nsref,nequit,time,timep,lprint,nprint
      common/range/mft,mlt,lft,llt,nftm1
      dimension freep(5,*),mtype(*),beta(*),matp(*),
     1 lmm(44*iprec,*),s(44,*)
      do 50 i=lft,llt
      elbrth=freep(3,i)
      eldeth=freep(4,i)
      elbury=freep(5,i)
      if (time.lt.elbrth.or.mtype(i).eq.0)
     1call stfset (s(1,i),nftm1+i,mtype(i),beta(i),i)
      if (time.le.eldeth) go to 50
      call scales (s(1,i),i,matp(i),time,eldeth,elbury)
   50 continue
      return
      end

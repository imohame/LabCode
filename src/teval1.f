      subroutine teval1(temp,npc,p,lcc,tt)
c     implicit double precision (a-h,o-z)                                    dp
c
c...  evaluate nodal temperatures from a load curve
c
      common/bk14/lfna(15),lfnt(6)
      common/cn1/numati,numnpi,numeli,nblk1,nslidi,ntslvi,ntmsri,
     1           nnpbi,nepbi,ncnpi
      dimension temp(*),npc(*),p(*)
c
      ierr=0
      xmag=1.0
      loc=npc(lcc)
      npoint=(npc(lcc+1)-loc)/2
      f=0.0
      call interp (p(loc),tt,npoint,f,xmag,ierr)
      if(ierr.eq.1)then
        write(lfnt(2),10)
        call bye(2)
      else
      do 20 i=1,numnpi
      temp(i)=f
   20 continue
      endif
c
      return
   10 format(//' ***Error*** time exceeds temperature load curve')
      end

      subroutine ldcset (fval,npc,p,tt)
c     implicit double precision (a-h,o-z)                                    dp
c
c     compute current value of load curves
c
      common/bk27/nlcur,nptst,nthpy,nthpz,nthps,xmy,xmz,xms,nload,nptm
      common/bks17/fmult(5)
      dimension fval(*),npc(*),p(*)
c
      if (nlcur.eq.0) return
c
      do 10 n=1,nlcur
      fval(n+nlcur)=fval(n)
      ierr=0
      xmag=1.0
      loc=npc(n)
      npoint=(npc(n+1)-loc)/2
      call interp (p(loc),tt,npoint,fvl,xmag,ierr)
      fval(n)=fvl
c.... island (loadset command)
      if(n.le.5)fval(n)=fval(n)*fmult(n)
      if (ierr.eq.1) fval(n)=0.0
   10 continue
c
      return
c
      end

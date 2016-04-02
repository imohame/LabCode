      subroutine lnesch(ui,usi,r,tvc2,u,tvc1,tt,step,g,g0,rlnew,riksf,
     1 alfa)
c     implicit double precision (a-h,o-z)                                    dp
c
      common/bk09/maxref,rhsn,rhsvn,cvtl,iteref,ectl,tolls
      common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
      dimension ui(*),usi(*),tvc2(*),r(*),u(*),tvc1(*)
c
      ga=g
      gb=g0
      sb=0.0
      sa=1.0
      do 20 j=1,10
      step=sa-ga*(sa-sb)/(ga-gb)
      call assmrs (ui,usi,tvc2,r,u,tvc1,tt,step,0,g)
      gb=0.50*gb
      if (g*1.0e+20*ga.gt.0.0) go to 10
      sb=sa
      gb=ga
   10 sa=step
      ga=g
      if (abs(g).gt.  tolls*abs(g0)) goto 20
      if (abs(sb-sa).gt..45*(sa+sb)) go to 20
      return
   20 continue
      return
      end

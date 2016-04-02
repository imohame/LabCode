      subroutine loadpr(id,u,r,npc,p,yn,zn,pmult,nodes,strt,lc,time)
c     implicit double precision (a-h,o-z)                                    dp
c
c     calculates pressure and shear tractions
c
      common/bk18/nummat,ityp2d,ako(31)
      common/bk30/numlp,numpc,h22(2,2),pl2(2,2),h33(3,2),pl3(3,2)
      common/bks17/fmult(5)
c
      dimension id(2,*),u(*),r(*),yn(3,*),zn(3,*),pmult(4,*),npc(*),
     1          nodes(3,*),strt(*),p(*),xx1(3),xx2(3),h(3),q(3),lc(*)
c
      if (numpc.eq.0) return
      do 140 n=1,numpc
      ierr=0
      xmag=1.0
      ishear=0
      lcc=lc(n)
      if (lcc.gt.0) go to 10
      ishear=1
      lcc=-lcc
   10 loc=npc(lcc)
      npoint=(npc(lcc+1)-loc)/2
      tt=time-strt(n)
      if(tt.lt.0.) go to 140
c.... island (loadset command)
      if(lcc.le.5)xmag=xmag*fmult(lcc)
      call interp (p(loc),tt,npoint,f,xmag,ierr)
      if (ierr.eq.1) go to 140
      xx1(1)=yn(1,n)
      xx2(1)=zn(1,n)
      xx1(2)=yn(2,n)
      xx2(2)=zn(2,n)
      k=id(1,nodes(1,n))
      l=id(2,nodes(1,n))
      if (k.ne.0) xx1(1)=xx1(1)+u(k)
      if (l.ne.0) xx2(1)=xx2(1)+u(l)
      k=id(1,nodes(2,n))
      l=id(2,nodes(2,n))
      if (k.ne.0) xx1(2)=xx1(2)+u(k)
      if (l.ne.0) xx2(2)=xx2(2)+u(l)
      if (.5*(xx1(1)+xx1(2)).gt.pmult(3,n)) go to 140
      if (.5*(xx2(1)+xx2(2)).gt.pmult(4,n)) go to 140
      do 130 i=1,2
      dy=pl2(1,i)*xx1(1)+pl2(2,i)*xx1(2)
      dz=pl2(1,i)*xx2(1)+pl2(2,i)*xx2(2)
      pr=h22(1,i)*pmult(1,n)+h22(2,i)*pmult(2,n)
      prs=pr*f
      if (ishear.eq.1) then
      tyy=prs*dy
      tzz=prs*dz
      else
      tyy=-prs*dz
      tzz=prs*dy
      endif
      if (ityp2d.eq.0) then
      yy=h22(1,i)*xx1(1)+h22(2,i)*xx1(2)
      tyy=tyy*yy
      tzz=tzz*yy
      endif
      iy=id(1,nodes(1,n))
      iz=id(2,nodes(1,n))
      if (iy.ne.0) r(iy)=r(iy)+h22(1,i)*tyy
      if (iz.ne.0) r(iz)=r(iz)+h22(1,i)*tzz
      iy=id(1,nodes(2,n))
      iz=id(2,nodes(2,n))
      if (iy.ne.0) r(iy)=r(iy)+h22(2,i)*tyy
      if (iz.ne.0) r(iz)=r(iz)+h22(2,i)*tzz
  130 continue
  140 continue
      return
      end

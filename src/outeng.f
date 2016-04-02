      subroutine outeng(u,udt,udtt,id,numnp,energy,xmass,matp,ix,
     1 matype,ro,temp,knetic,intrnl,total,ener,yz,y,z,thick,ln)
c     implicit double precision (a-h,o-z)                                    dp
c
c     write displacements and velocities into plot file
c
      common/bk05/ifil,iadd,maxsiz,head(12)
      common/bk06/nprnt,mprint,itmpop,numelt,jprint,idump,locstr
      common/bk14/lfna(15),lfnt(6)
      common/bk16/maxint,hgc
      common/bk18/nummat,ityp2d,ako(31)
      common/bk20/ntotal
      common/bk32/nsref,nequit,time,timep,lprint,nprint
      common/bk34/b(1)
      common/cn3/ibbari,intgi,nmbfi,ithopi,ithcri,ithini,iengri,
     1           ijinti
      dimension u(*),udt(*),udtt(*),id(2,*),d(2),energy(ln,*),xmass(4,*)
     1 ,matp(*),ix(4,*),matype(*),ro(*),temp(*),ener(*),y(*),z(*),
     2 yz(2,*),xx(2,4),thick(ln,*)
      real knetic,intrnl
c
      do 50 i=1,nummat
   50 ener(i)=0.
c
      loc    =0
      knetic=0.0
      intrnl=0.0
      lce=1-maxint
      do 80 i=1,numelt
      lce=lce+maxint
      do 66 k=1,4
      xx(1,k)=y(ix(k,i))
      xx(2,k)=z(ix(k,i))
   66 continue
      matpi=matp(i)
      loc=loc+1
      call engcal(energy(1,lce),thick(1,lce),xx,b(loc),
     1 maxint,ln,vol0)
      intrnl=intrnl+b(loc)
      if (matpi.eq.0) go to 70
      ener(matpi)=ener(matpi)+b(loc)
      do 200 j=1,4
      id1=id(1,ix(j,i))
      id2=id(2,ix(j,i))
      if (id1.ne.0) knetic=knetic+xmass(j,i)*udt(id1)**2
      if (id2.ne.0) knetic=knetic+xmass(j,i)*udt(id2)**2
  200 continue
      if(iengri.eq.0)then
      b(loc)=b(loc)/vol0
      elseif(iengri.eq.2)then
      tempe=0.
      do 60 j=1,4
      tempe=tempe+.25*temp(ix(j,i))
   60 continue
      b(loc)=tempe
      endif
   70 if (loc+1.lt.ntotal) go to 80
ck      call wrabsg(lfna(7),b,loc,iadd,b(ntotal+1),0)
ck      call riosta(lfna(7))
ck      iadd=iadd+loc
ck      loc=0
   80 continue
      knetic=knetic/2.
      total  =knetic+intrnl
      if (loc.eq.0) go to 90
ck      call wrabsg(lfna(7),b,loc,iadd,b(ntotal+1),0)
ck      call riosta(lfna(7))
ck      iadd=iadd+loc
   90 do 100 i=1,numnp
      yz(1,i)=y(i)
  100 yz(2,i)=z(i)
ck      call wrabsg(lfna(7),matp,numelt,iadd,b(ntotal+1),1)
ck      call riosta(lfna(7))
ck      iadd=iadd+numelt
ck      call wrabsg(lfna(7),yz,2*numnp,iadd,b(ntotal+1),0)
ck      call riosta(lfna(7))
ck      iadd=iadd+2*numnp
ck      phony=1.e+21
ck      call wrabsg(lfna(7),phony,1,iadd,b(ntotal+1),0)
ck      call riosta(lfna(7))
c
      return
      end

      subroutine printv (freq,phi,id,numnp,neq,yz,y,z,matp)
c     implicit double precision (a-h,o-z)                                    dp
c
      common/bk00/
     1k01,k02,k03,k04,k05,k06,k07,k08,k09,k10,k11,k12,
     2k13,k14,k15,k16,k17,k18,k19,k20,k21,k22,k23,k24,
     3k25,k26,k27,k28,k29,k30,k31,k32,k33,k34,k35,k36,
     4k37,k38,k39,k40,k41,k42,k43,k44,k45,k46,k47,k48,
     5k49,k50,k51,k52,k53,k54,k55,k56,k57,k58,k59,k60,
     6k61,k62,k63,k64,k65,k66,k67,k68,k69,k70,k71,k72,
     7k73,k74,k75,k76,k77,k78,k79,k80,k81,k82,k83,k84,
     8k85,k86,k87,k88,k89,k90,k91,k92,k93,k94,k95,k96
      common/bk05/ifil,iadd,maxsiz,head(12)
      common/bk06/nprnt,mprint,itmpop,numelt,jprint,idump,locstr
      common/bk14/lfna(15),lfnt(6)
      common/bk18/nummat,ityp2d,ako(31)
      common/bk20/ntotal
      common/bk29/numfrq,clengt
      common/bk34/bb(1)
      common/main_block/ a(1)
      dimension freq(*),phi(*),id(2,*),d(2),yz(2,*),y(*),z(*),tim(5),
     1 matp(1)
      data twopi/6.283185307/
c
      do 10 i=1,numfrq
   10 freq(i)=sqrt(1./freq(i))
      call header
      write(lfnt(2),100)
      do 20 i=1,numfrq
      acirc=freq(i)/twopi
      aperd=1./acirc
      write(lfnt(2),110) i,freq(i),acirc,aperd
   20 continue
c
      clengt=clengt/20.0
      nume20=20*numelt
      do 90 jj=1,numfrq
      phimax=0.0
      nprnt=0
      call blkcpy (a(k16+neq*(jj-1)),phi,neq)
c
c     normalize eigenvectors
c
      do 30 ii=1,neq
      if (abs(phi(ii)).gt.phimax) phimax=abs(phi(ii))
   30 continue
      phimax=1.0/phimax
      do 40 ii=1,neq
   40 phi(ii)=phimax*phi(ii)
c
      do 70 ii=1,numnp
      d(1)=0.
      d(2)=0.
      do 50 i=1,2
      kk=id(i,ii)
      if (kk.ne.0) d(i)=phi(kk)
   50 continue
      if (nprnt.gt.0) go to 60
      nprnt=50
      call header
      write(lfnt(2),120) jj,freq(jj)
   60 nprnt=nprnt-1
      write(lfnt(2),130) ii,d
   70 continue
c     call empty (17)                                                   cray1
c
c     scale eigenvectors before writing them into plotfile
c
      do 80 ii=1,neq
   80 phi(ii)=clengt*phi(ii)
      tim(1)=freq(jj)
      ifctr=iadd/maxsiz
      icnt=5+nummat+8*numnp+(21+idump)*numelt
      if (iadd+icnt.gt.(ifctr+1)*maxsiz) iadd=(ifctr+1)*maxsiz
      write (lfnt(2),3000) jj,nume20,iadd
 3000 format (3h xx,3i10)
c     call empty (17)                                                   cray1
ck      call wrabsg(lfna(7),tim,5,iadd,bb(ntotal+1),0)
ck      call riosta(lfna(7))
ck      iadd=iadd+5+nummat+nume20
      call outdv (phi,id,numnp)
ck      call wrabsg(lfna(7),matp,numelt,iadd,bb(ntotal+1),1)
ck      call riosta(lfna(7))
ck      iadd=iadd+(1+idump)*numelt
      do 85 i=1,numnp
      yz(1,i)=y(i)
   85 yz(2,i)=z(i)
ck      call wrabsg(lfna(7),yz,2*numnp,iadd,bb(ntotal+1),0)
ck      call riosta(lfna(7))
ck      iadd=iadd+2*numnp
c
   90 continue
c     tim(1)=0.0
c     call wrabsg(lfna(7),tim,5,iadd,bb(ntotal+1),0)
c     call riosta(lfna(7))
c
      return
c
  100 format(///' r e s u l t s   o f   e i g e n v a l u e',
     1 '   a n a l y s i s'/)
  110 format('frequency number=',i3,' frequency (radians)=',e12.5,
     1 ' frequency (hertz)=',e12.5,' period=',e12.5)
  120 format (///'   m o d e  s h a p e  n o .  ',i3,
     1       55x,'( frequency =  ',e11.4, ' )' //
     2       '    node',13x,'y-displacement',4x,'z-displacement')
  130 format (2x,i5,8x,2e18.6)
      end

      subroutine elemin(numel,ix,betan,freep,matp,matype,id,ifree,
     1 tref,reftem,temmat,nit,nt)
c     implicit double precision (a-h,o-z)                                    dp
c
c     read and store element data

      common /wblock20/ mat_type(40000)
c
      common/bk14/lfna(15),lfnt(6)
      common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
      common/cn0/iconv
      common/cn3/ibbari,intgi,nmbfi,ithopi,ithcri,ithini,iengri,
     1           ijinti
      dimension ix(4,*),betan(*),freep(5,*),matp(*),matype(*),id(2,*),
     1 ifree(*),ic(9),iy(9),tref(4,*),reftem(*),temmat(*)
      character*80 txts,mssg
      character*2 nt
      logical nit
      data ic/0,0,0,0,0,0,0,0,0/
!!!!!!!!      common /WMLthermal/thermalflag
!!!!!!!!!!!      common /WMLthermali/connect(40000,4)
!!!!!!!!!!!      integer connect
!!!!!!!      integer thermalflag

c
!!!!    ----------------------------testing  
!!        write(*,*)id(1,1:25)
!!        write(*,*)id(2,1:25)
        
        
      nprnt=0
      mxnods=0
      ips=0
      i=0
      do 10 l=1,numnp
   10 ifree(l)=0
   20 inc=inc0
      imt=ic(5)
      bfmgy=0.0
      bfmgz=0.0
      elbrth=0.0
      eldeth=1.e20
      elbury=1.e20
      if(nt.ne.'91')then
      iconsv=iconv
      iconv=0
      call gttxsg (txts,lcount)
      if(.not.nit)then
      read(unit=txts,fmt=190,err=165) m,(ic(j),j=1,9),inc0,beta,irf
      else
      read(unit=txts,fmt=250,err=165) m,(ic(j),j=1,5),inc0,beta,irf
      endif
      if (irf.eq.1) then
      call gttxsg (txts,lcount)
      read(unit=txts,fmt=230,err=165) bfmgy,bfmgz,elbrth,eldeth,elbury
      endif
      if(iconsv.eq.1)then
ck      write(lfnt(4),251)m,(ic(j),j=1,5),inc0,beta,elbrth,eldeth,elbury
      endif
      iconv=iconsv
      else
      irf=2
      call gttxsg (txts,lcount)
      read(unit=txts,fmt=252,err=165)
     1  m,(ic(j),j=1,5),inc0,beta,elbrth,eldeth,elbury
      endif
      if (elbrth.le.0.)     elbrth=-1.e+10
      if (eldeth.le.elbrth.or.eldeth.le.0.) eldeth=1.e+20
      if (elbury.le.eldeth) elbury=1.e-09+eldeth
      if (ic(5).eq.0) go to 30
      model=matype(ic(5))
      if (beta.eq.0.0.and.model.eq.5)  beta=1.0
      if (beta.eq.0.0.and.model.eq.14) beta=1.0
      if (beta.eq.0.0.and.model.eq.21) beta=1.0
   30 if (inc.eq.0) inc=1
   40 i=i+1
      if (m-i) 140,70,50
   50 do 60 l=1,4
   60 iy(l)=iy(l)+inc
      iy(5)=imt
      go to 90
   70 do 80 j=1,9
   80 iy(j)=ic(j)
      btta=beta
   90 n=i
      ix(1,n)=iy(1)
      ix(2,n)=iy(2)
      ix(3,n)=iy(3)
      ix(4,n)=iy(4)
      matp(n)=iy(5)
  
!!!!!!c     passing connectivity on to thermal calcuations WML 82010
!!!!!!      if (thermalflag > 0) then
!!!!!!c      write(*,*)'thermalflag=',thermalflag
!!!!!!      connect(n,1)=ix(1,n)
!!!!!!      connect(n,2)=ix(2,n)
!!!!!!      connect(n,3)=ix(3,n)
!!!!!!      connect(n,4)=ix(4,n)
!!!!!!      endif


      mat_type(n) = iy(5)

      betan(n)=btta
      freep(1,n)=bfmgy
      freep(2,n)=bfmgz
      freep(3,n)=elbrth
      freep(4,n)=eldeth
      freep(5,n)=elbury
      if (nprnt.gt.0) go to 100
      nprnt=50
      call header
      write(lfnt(2),210) (j,j=1,4)
  100 nprnt=nprnt-1
      write(lfnt(2),220) n,btta,iy(5),(iy(j),j=1,4)
c      write(29,270) n,(iy(j),j=1,4)
c.... define element reference temperature from nodal values
      if(nt.ne.'91')then
        tref(1,n)=.25*(reftem(iy(1))+reftem(iy(2))+reftem(iy(3))
     2         +reftem(iy(4)))
      else
        tref(1,n)=temmat(matp(n))
      endif
      if(iy(5).eq.0) go to 110
      ifree(iy(1))=1
      ifree(iy(2))=1
      ifree(iy(3))=1
      ifree(iy(4))=1
  110 continue
      if(irf.eq.1)then
      nprnt=nprnt-1
cwail  write(lfnt(2),240) bfmgy,bfmgz,elbrth,eldeth,elbury
      elseif(irf.eq.2)then
      nprnt=nprnt-1
cwail  write(lfnt(2),241) elbrth,eldeth,elbury
      endif
  120 if (m-i) 140,130,40
  130 if (numel-i) 150,150,20
  140 write(lfnt(2),200) i
      write (*,200) i
      call bye (2)
c
  150 inun=0
      do 160 i=1,numnp
      if (ifree(i).ne.0) go to 160
      if (id(1,i).ne.0.and.id(2,i).ne.0) go to 160
      id(1,i)=1
      id(2,i)=1
      write(lfnt(2),170) i
      inun=inun+1
  160 continue
      if (inun.ne.0) write (*,180) inun
!!!!    ----------------------------testing  
!!      write(*,*)matp(1:1+15)
      
      return
c
  165 j=i+1
      write (unit=mssg,fmt=260) j
      call termin(txts,mssg,lcount,1)
c
  170 format(/'node ',i5,' is not connected to any element'/
     1 'its boundary condition code has been reset to 3')
  180 format(///'warning ',i4,' disconnected nodes found in mesh'/
     1 'they have been constrained- see output file for node numbers')
  190 format(11i5,e15.0,i5)
  200 format(/' error in element data',i5)
  210 format(///4x,'e l e m e n t   i n f o r m a t i o n '/
     1 /'     m          bet        mtyp',5x,4('node',i1,3x))
  220 format(1x,i5,8x,e10.2,i7,4(4x,i4))
  230 format(6e10.0)
  240 format(' magnetic body forces   x=',e11.3,'   z=',e11.3,5x,
     1 'birth=',e11.3,'   death=',e11.3,'   burial=',e11.3)
  241 format( 'birth=',e11.3,'   death=',e11.3,'   burial=',e11.3)
  250 format(7i5,e10.0,i5)
  251 format(7i5,4e10.3)
  252 format(7i5,4e10.0)
  260 format(' error reading element cards, probably element#',i6)
  270 format(i5,4(4x,i4))
      end

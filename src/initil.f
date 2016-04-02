      subroutine initil
c     implicit double precision (a-h,o-z)                                    dp
c
      use mod_parameters
      real*8 hed                                                        
      common/bk00/
     1k01,k02,k03,k04,k05,k06,k07,k08,k09,k10,k11,k12,
     2k13,k14,k15,k16,k17,k18,k19,k20,k21,k22,k23,k24,
     3k25,k26,k27,k28,k29,k30,k31,k32,k33,k34,k35,k36,
     4k37,k38,k39,k40,k41,k42,k43,k44,k45,k46,k47,k48,
     5k49,k50,k51,k52,k53,k54,k55,k56,k57,k58,k59,k60,
     6k61,k62,k63,k64,k65,k66,k67,k68,k69,k70,k71,k72,
     7k73,k74,k75,k76,k77,k78,k79,k80,k81,k82,k83,k84,
     8k85,k86,k87,k88,k89,k90,k91,k92,k93,k94,k95,k96
      common/bk01/h4(4,5),p14(4,5),p24(4,5)
      common/bk02/ioofc,iphase,imass,lpar(9)
      common/bk04/nprm(8)
      common/bk06/nprnt,mprint,itmpop,numelt,jprint,idump,locstr
      common/bk07/mbfc,nelpg,hed(12)
      common/bk08/kprint,nstep,ite,ilimit,newstf
      common/bk11/cnwmk(2),iequit,iprint,isref
      common/bk14/lfna(15),lfnt(6)
      common/bk15/cpuio(36),cpuip(36)
      common/bk16/maxint
      common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
      common/bk18/nummat,ityp2d,ako(31)
      common/bk23/itemp,itherm,irtin
      common/bk26/dt,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10
      common/bk27/nlcur,nptst,nthpy,nthpz,nthps,xmy,xmz,xms,nload,nptm
      common/bk29/numfrq,clengt
      common/bk30/numlp,numpc,h22(2,2),pl2(2,2),h33(3,2),pl3(3,2)
      common/bk32/nsref,nequit,time,timep,lprint,nprint
      common/bk34/bb(1)
      common/bk44/hh(8),pp(2,8),xjj(4)
      common/bk45/gps(5),gpt(5)
      common/slar3/nsl,nsntl,nmntl,nslnmx,sltol,slhrd
      common/fissn0/maxneq,mwspac,ntpe0,ntpe1,nfissl(3)
      common/fissn1/melemt,nnns,ntpe2,n2g,llls
      common/fissn2/nwpblk,numblk,mwsusd,mxnepb,maxch,matpr,mench,ifa(2)
      common/cn3/ibbari,intgi,nmbfi,ithopi,ithcri,ithini,iengri,
     1           ijinti
      common/main_block/ a(1)
	  common/mbsize/numelt2, numnp2, neq2
c
      dimension gpst(10)
      equivalence (gpst,gps)
c
      data gpst/.5773502691896,-.5773502691896,-.5773502691896,
     1   .5773502691896, .0000000000000, .5773502691896, .5773502691896,
     2  -.5773502691896,-.5773502691896, .0000000000000/
c
cw      call timin (cpuio,cpuip,1,3)
cw      call timin (cpuio,cpuip,3,2)
      if (itemp.eq.0) go to 10
      numnp1=numnp+1
   10 lpar(5)=ityp2d
      loc=1
      lpar(3)=nelpg
      lpar(4)=numelt
      ndm=8
c
      nprm(1)=ndm
      nprm(2)=(ndm*(ndm+1))/2
      nprm(3)=nprm(2)+ndm
      nprm(4)=0
      nprm(5)=4*lpar(9)
      nprm(6)=48
      nprm(7)=8
      nprm(8)=1
c
c     compute and store basis functions for 4 node elements
c
      do 30 lst=1,5
      e1=gps(lst)
      e2=gpt(lst)
      call basis1 (e1,e2,hh,pp)
      do 20 i=1,4
      h4(i,lst)=hh(i)
      p14(i,lst)=pp(1,i)
   20 p24(i,lst)=pp(2,i)
   30 continue
c
c----------------------------------------------------------------
c----------------------------------------------------------------
ck      write(25,12) maxint
ck   12 format(////5x,'Gaussian integeration points =',i2)  
ck      write(25,35)
ck   35 format(////5x,'Shape function for 4 integration points :'//)
ck      do 37 j=1,4
ck      write(25,36) j,(h4(i,j),i=1,4)
ck   36 format(5x,i1,' point : ',4(f8.3)/)
ck   37 continue
c
ck      write(25,38)
ck   38 format(///5x,'Derivative of shape function respect to local'
ck     .       ,' coordinates :'//)
ck      do 42 j=1,4
ck      write(25,41) j,(p14(i,j),i=1,4),j,(p24(i,j),i=1,4)
ck   41 format(5x,i1,' point respect to yy : ',4(f8.3)/5x,
ck     .          i1,' point respect to zz : ',4(f8.3)//)
ck   42 continue
c----------------------------------------------------------------   
c----------------------------------------------------------------
c.... set load curves at time zero
      call ldcset (a(k57+2*numnp2),a(k63),a(k64),0.0)
c
c.... set initialization step for reference temperatures
      if(imass.eq.2.and.ithopi.ne.0)then
        write(lfnt(2),199)
        call bye(2)
      endif
      if(ithini.ne.0.and.ithopi.ne.-3)then
      time  =-dt
      nstep =-2
      elseif(imass.eq.2)then
      time=-dt
      nstep=-2
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c     else
c.... set initialization step if any nonzero loads at t=0
c     if(mthsol.lt.6) then
c     do 50 i=1,nlcur
c     if (ithcri.eq.i) go to 50
c     if (a(k57+2*numnp-1+i).eq.0.0) go to 50
c     time  =-dt
c     nstep =-2
c     go to 60
c  50 continue
c     endif
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      endif
c
c.... initialize temperatures
   60 continue
      call intmp(a(k81),a(k81+1),a(k82),a(k82+1))
      lpar(2)=numelt
c
      call intelt
c
c     set variables in block solver's common blocks
c
      maxneq=neq
      nnns=ndm
      llls=nnns+nnns*(nnns+1)/2
      n2g=128
c
c     initialize a array
c
      call blkcpy (a(k13),a(k16),nwebuf)
      call azero (a(k20),k21-k20)
c
cw      call timin (cpuio,cpuip,3,3)
cw      call timin (cpuio,cpuip,5,2)
c
c      if (nsl.ne.0) call intsld
c
c.... j-integral initialization
c     if(ijinti.ne.0)call jconsr
c
      return
  199 format(' **Error** cannot use static initialization option'/
     1 ' with thermal option')
      end

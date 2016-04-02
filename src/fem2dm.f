      subroutine fem2dm
c     implicit double precision (a-h,o-z)                                    dp
c
c     ovrlay for 2-dimensional solid elements
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
      common/bk02/ioofc,iphase,imass,lpar(9)
      common/bk06/nprnt,mprint,itmpop,numelt,jprint,idump,locstr
      common/bk12/ntlen
      common/bk16/maxint,hgc
      common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
      common/bk18/nummat,ityp2d,ako(31)
      common/main_block/ a(1)
	  common/mbsize/numelt2, numnp2, neq2
c
cwq      m1=k06+lpar(2)
      m1=k06+numelt2
      m2=k18		!u
      m3=k19
      m4=k57+2*numnp2
      iequit=0
      if (iphase.ge.3)then
      iequit=1
      m2=k22		!tvc1
      m3=k21		!
      endif
c
      call f2d (a(k17),a(m1),a(m2),a,a(m3),a(k20),a(k02),a(k03),a(k04),
     1a(k06),a(k07),a(k08),a(k10),a(k11),a(k09),a(k01),
     2 a(k55),a(k56),a(m4),iequit,a(ntlen),a(k57))                           nk
c    2 a(k55),a(k56),a(m4),iequit,a(ntlen),a(k57),a(k07+5*numelt))           pl
c
c.... phase change
      if(iphase.eq.4)iphase=2
c
      return
      end

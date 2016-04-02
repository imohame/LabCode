      subroutine assmas
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
      common/bk02/ioofc,iphase,imass,lpar(9)
      common/bk04/ndm,nmt,nms,nd5dim,idwa,ncon,idw,itdim
      common/main_block/ a(1)
      equivalence (lpar(7),mxnods)
c

c         write(7777,*) '-- assmas.f'

c
      call fe2dm (a(k10),a(k63),a(k02),a(k08),a(k09),a(k57))
c
      return
      end

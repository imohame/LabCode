      subroutine eigen
c     implicit double precision (a-h,o-z)                                    dp
c
c     out-of-core eigenvalue solution by subspace iteration
c     is based on incore eigenvalue routine provided by
c     prof. r. l. taylor university of california, berkeley
c     extended and implemented in crystal2d by j.o. hallquist
c
      logical conv,plt
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
      common/bk12/ntlen
      common/bk14/lfna(15),lfnt(6)
      common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
      common/bk29/numfrq,clengt
      common/main_block/ a(1)

      data its/50/
c
c     assemble lumped mass vector
c
      call blkcpy (a(k16),a(k13),nwebuf)
      call lumps (a(k18),a(k57),a(k76),a(k77),neq)
      call setpnt (63,54,2*numnp)
      call blkcpy (a(k18),a(k63),neq)
c
      call assmas
c
      mm=2*numnp
      call blkcpy (a(k57),a(k54),mm)
      nv=min(neq,numfrq+8)
      neqnv=neq*nv
      call setpnt (18,17,neq )
      call setpnt (19,18,neq )
      call setpnt (20,19,2*numnp )
      call setpnt (21,20,neq )
      call setpnt (22,21,nv*(nv-1)/2+nv)
      call setpnt (55,22,nv*(nv-1)/2+nv)
      call setpnt (56,55,nv  )
      call setpnt (57,56,nv*nv   )
      call setpnt (63,57,nv*nv   )
      call setpnt (64,63,nv  )
      call setpnt (65,64,nv+1)
      call setpnt (66,65,nv  )
c
      call blkcpy (a(k54),a(k19),mm+neq)
c
      call pzero (a(k18),k19-k18)
      call pzero (a(k22),k66-k22)
c
c     assemble and factor stiffness matrix
c
      call bsolvr (a(k18),a(ntlen),2,6)
c
c     compute the initial iteration vectors
c
      call intvec (a(k20),a(k18),neq,nv)
c
      call header
      write(lfnt(2),40)
      conv=.false.
      do 20 iter=1,its
c
c     setup reduced eigenvalue problem g*p=h*p*d  (scott pencil)
c
      call formm (a(k17),a(k18),a(k22),a(k20),neq,nv)
c
      do 10 i=1,nv
      j=(i-1)*neq
      call blkcpy (a(k16+j),a(k18),neq)
      call bsolvr (a(k18),a(ntlen),3,6)
   10 call blkcpy (a(k18),a(k54+j),neq)
c
      call forms (a(k17),a(k18),a(k21),neq,nv)
c
c     solve the reduced general eigenproblem
c
      call geig (a(k21),a(k22),a(k55),a(k56),a(k57),a(k63),a(k64),nv)
c
c     check for convergence
c
      call convck (a(k55),a(k65),a(k64),conv,nv,iter)
c
c     multiply eigenvectors by eigenvalue to prevent overflows
c
      call multev (a(k55),a(k56),nv)
c
c     compute new iteration vectors
c
      call newvec (a(k18),a(k56),a(k17),nv,neq)
      if (conv) go to 30
   20 continue
c
c     print the frequencies and mode shapes
c
   30 call printv (a(k55),a(k18),a(k19),numnp,neq,a(k13),a(k03),a(k04),
     1 a(k08))
c
      call bye (1)
c
   40 format(///' e i g e n v a l u e   a n a l y s i s   b y   s u b s
     1p a c e   i t e r a t i o n'/)
      end

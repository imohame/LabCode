      subroutine addem(je,ke,le,ia,a,ld,s,ns,ls,ls2,neg)
c     implicit double precision (a-h,o-z)                                    dp
c
c     add contributions from a group of elements to two blocks of eqs.
c
       use mod_parameters
       dimension ia(*),a(*),ld(ls2,1),s(ls,1)
c
c     input arguments
c           je        eqn. number of the first eqn. in the block -1
c           ke        eqn. number of the last eqn. in the block
c           le        eqn. number of the last eqn. in the next block
c           ia        array of diagonal pointers for the blocks
c            a        array containing blocks of equations
c           ld        equation numbers for each element's degrees of
c                     freedom (a zero equation number implies that the
c                     degree of freedom is inactive)
c            s        upper triangle of each element coefficient matrix,
c                     stored columnwise
c           ns        number of degrees of freedom per element
c           ls        length of an element record
c           neg       number of elements in the group
c
c     element (i,j) of the element coefficient matrix is added (mapped)
c     to element (ld(i),ld(j)) of the system coefficient matrix. this
c     element is stored in a(k), where k=ia(ld(j)-je)-(ld(j)-ld(i)),
c     and ld(j) is greater than or equal to ld(i). if ld(j) is greater
c     than ke, the equation is in the second block and k is set equal
c     to ia(ld(j)-ke+nwb)-(ld(j)-ld(i))+nwb.
c
c     output arguments
c            a        updated blocks of equations (system coef. matrix)
c
      common/fissn2/nwb,nblk,nwu,meb,mxc,mpr,mch,ia1,ia2
c
c     note- zero origin addresses are passed to this subroutine for
c     the arrays ia and a. thus ia(i) must be referenced as ia(i+1)
c     and a(i) must be referenced as a(i+1).
c
c.... scratch arrays for vectorization on groups of 128 elements
      common/fissl6/k(nelemg),l(nelemg),m(nelemg),ss(nelemg)
      common/double/iprec,ncpw,unit

c
      nw=nwb*iprec
      me=ke-je
      ne=le-je
      ij=0
      do 200 j=1,ns
      do 180 i=1,j
      do 100 n=1,neg
      k(n)=max(ld(i,n),ld(j,n))-je
      l(n)=abs(ld(i,n)-ld(j,n))
c     m(n)=or(or(k(n),ld(i,n)-1),ld(j,n)-1)                             cray1
c     k(n)=cvmgm(0,k(n),or(m(n),ne-k(n)))                               cray1
      if((k(n).lt.0).or.(ld(i,n).eq.0).or.(ld(j,n).eq.0)                vax75
     1 .or.(k(n).gt.ne))k(n)=0                                          vax75
  100 continue
      if(i.ne.j) go to 121
      do 120 n=1,neg
      ss(n)=s(i+ij,n)
  120 continue
      go to 141
  121 do 130 n=1,neg
c     ss(n)=cvmgz(s(i+ij,n)+s(i+ij,n),s(i+ij,n),l(n))                   cray1
      ss(n)=s(i+ij,n)                                                   vax75
      if(l(n).eq.0)ss(n)=ss(n)+ss(n)                                    vax75
  130 continue
  141 if(me.eq.ne) go to 151
      do 150 n=1,neg
      m(n)=me-k(n)
c     k(n)=cvmgp(k(n),nwb-m(n),m(n))                                    cray1
c     l(n)=cvmgp(l(n),l(n)-nwb,m(n))                                    cray1
      if(m(n).lt.0)then                                                 vax75
      k(n)=nw-m(n)                                                      vax75
      l(n)=l(n)-nwb                                                     vax75
      endif                                                             vax75
  150 continue
  151 do 160 n=1,neg
      if(k(n).eq.0) go to 160
      kl = ia(k(n)) - l(n)
      a(kl) = a(kl) + ss(n)
ck----------------------------------------------
ck----------------------------------------------
ck      write(*,*) kl,a(kl)
ck  152 format(5x,' k( ',i3,' )= ',e14.6)
ck----------------------------------------------
ck----------------------------------------------  
  160 continue
  180 continue
  200 ij=ij+j
      return
      end










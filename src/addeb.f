       subroutine addeb(ia,a,ld,s,ns,ls,ls2,neg)
c     implicit double precision (a-h,o-z)                                    dp
c
c      add contributions from a group of elements to global equations
c
       use mod_parameters
       dimension ia(1),a(1),ld(ls2,1),s(ls,1)
c
c      input arguments
c            ia        array of diagonal pointers for the equations
c             a        array containing system coefficient matrix
c            ld        equation numbers for each element's degrees of
c                      freedom (a zero equation number implies that the
c                      degree of freedom is inactive)
c             s        upper triangle of each element coefficient matrix
c                      stored columnwise
c            ns        number of degrees of freedom per element
c            ls        length of an element record
c            neg       number of elements in the group
c
c      element (i,j) of the element coefficient matrix is added (mapped)
c      to element (ld(i),ld(j)) of the system coefficient matrix. this
c      element is stored in a(k), where k=ia(ld(j))-(ld(j)-ld(i)), and
c      ld(j) is greater than or equal to ld(i).
c
c      output arguments
c             a        updated equations (system coef. matrix)
c
c      note- zero origin addresses are passed to this subroutine for
c      the arrays ia and a. thus ia(i) must be referenced as ia(i+1)
c      and a(i) must be referenced as a(i+1).
c
c      note- this routine has been vectorized for the cray-1. statements
c      flagged with 'cray' in columns 74-77 apply to the cray only.
c      statements flagged with 'ansi' in columns 74-77 are ansi standard
c      fortran equivalent to the cray coding.
c
cc.... scratch arrays for vectorization on groups of 128 elements
       common/fissl6/k(nelemg),l(nelemg),ss(nelemg)
c
       ij=0
       do 200 j=1,ns
       do 180 i=1,j
       do 100 n=1,neg
       k(n)=max(ld(i,n),ld(j,n))
       l(n)=abs(ld(i,n)-ld(j,n))
       if((ld(i,n).eq.0).or.(ld(j,n).eq.0)) k(n)=0
  100 continue
       if(i.ne.j) go to 121
       do 120 n=1,neg
       ss(n)=s(i+ij,n)
  120 continue
       go to 141
  121 do 130 n=1,neg
       ss(n)=s(i+ij,n)
       if(l(n).eq.0) ss(n)=ss(n)+ss(n)
  130 continue
  141 continue
       do 160 n=1,neg
       if(k(n).eq.0) go to 160
       kl=ia(k(n))-l(n)
       a(kl+1)=a(kl+1)+ss(n)
  160 continue
  180 continue
  200 ij=ij+j
       return
       end

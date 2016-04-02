      subroutine coorot(y,z,yz,numnp)
c     implicit double precision (a-h,o-z)                                    dp
c
c     write nodal coordinates on plot tape
c
      common/bk05/ifil,iadd,maxsiz,head(12)
      common/bk14/lfna(15),lfnt(6)
      common/bk20/ntotal
      common/bk34/bb(1)
      dimension y(*),z(*),yz(*)
c
      do 10 j=1,numnp
      yz(2*j-1)=y(j)
      yz(2*j)=z(j)
   10 continue
      nn2=2*numnp
ck      call wrabsg(lfna(7),yz(1),nn2,iadd,bb(ntotal+1),0)
ck      call riosta(lfna(7))
ck      iadd=iadd+nn2
c
      return
c
      end

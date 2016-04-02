      subroutine geig(g,h,d,p,s,t,id,nv)
c     implicit double precision (a-h,o-z)                                    dp
      dimension g(*),h(*),d(*),p(nv,1),s(nv,1),t(*),id(*)
      call jacobi (h,d,s,t,id,nv,ir)
      do 10 i=1,nv
      c=1./sqrt(d(i))
      id(i)=i*(i+1)/2
      do 10 j=1,nv
   10 s(j,i)=s(j,i)*c
      call pzero (p,nv*nv)
      do 20 i=1,nv
   20 call promul (g,s(1,i),p(1,i),id,nv)
      k=0
      do 30 j=1,nv
      do 30 i=1,j
      k=k+1
      h(k)=fdot(s(1,i),p(1,j),nv)
   30 continue
      call jacobi (h,d,p,t,id,nv,ir)
      call psorte (d,p,nv,.true.)
      do 60 j=1,nv
      do 40 i=1,nv
   40 t(i)=0.
      do 50 k=1,nv
      do 50 i=1,nv
   50 t(i)=t(i)+s(i,k)*p(k,j)
      do 60 i=1,nv
   60 p(i,j)=t(i)
      return
      end

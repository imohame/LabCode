      subroutine rezzap (numnp,nr,nrv,nd,ka)
c     implicit double precision (a-h,o-z)                                    dp
c***********************************************************************
c     this routine removes any zero degree nodes,
c     then sorts and packs the connectivity array
c***********************************************************************
      dimension nr(*),nrv(2,*), nd(2,*), ka(*), nodcon(300)
      common /sabr0/ maxnd
      common /sabr1/ mind, maxd, ibw, nbw, ipr, npr, nnc, nzn, nnp
      nzn=0
      do 10 i=1,numnp
      if (nd(1,i).ne.0) go to 10
      nzn=nzn+1
      nrv(1,nzn)=0
      nrv(2,nzn)=i
   10 continue
      ibw=0
      ipr=0
      nnc=0
      j=0
      nk=-maxnd
      mind=-nk
      maxd=0
      do 70 i=1,numnp
      nk=nk+maxnd
      if (nd(1,i).eq.0) go to 70
      l=nd(1,i)
      mind=min(mind,l)
      maxd=max(maxd,l)
      j=j+1
      nrv(1,j+nzn)=l
      nrv(2,j+nzn)=nnc
      do 15 k=1,l
      nodcon(k)=ka(nk+k)
   15 continue
      if (l.eq.1) go to 30
      do 20 k=2,l
      ndk=nd(1,nodcon(k-1))
      do 20 m=k,l
      ndm=nd(1,nodcon(m))
      if (ndk.le.ndm) go to 20
      ndk=ndm
      nt=nodcon(k-1)
      nodcon(k-1)=nodcon(m)
      nodcon(m)=nt
   20 continue
   30 continue
      if (nzn.eq.0) go to 50
      do 40 m=1,nzn
      do 40 k=1,l
      if (nodcon(k).lt.nrv(2,nzn-m+1)) go to 40
      nodcon(k)=nodcon(k)-1
   40 continue
   50 continue
      do 55 k=1,l
      ka(nnc+k)=nodcon(k)
   55 continue
      nnc=nnc+l
      kkh=0
      do 60 k=1,l
      kkh=max(kkh,j-nodcon(k))
   60 continue
      ipr=ipr+kkh
      ibw=max(ibw,kkh)
   70 continue
      nnc=nnc+1
      nnp=j
      nbw=ibw
      npr=ipr
      n=2*nnp
c     call blkcpy(nrv(1,nzn+1),nd,n)                                    cray1
      call blicpy(nrv(1,nzn+1),nd,n)                                    vax750
      if (nzn.eq.0) return
      do 80 i=1,nzn
   80 nr(i)=nrv(2,i)
      return
      end

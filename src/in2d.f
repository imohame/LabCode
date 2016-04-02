      subroutine in2d (ix,y,z,betap,matp,matype,den,nbck,ym,thick)
c     implicit double precision (a-h,o-z)                                    dp
c
      common/bk02/ioofc,iphase,imass,model,lpar(8)
      common/bk06/nprnt,mprint,itmpop,numelt,jprint,idump,locstr
      common/bk34/bb(1)
      common/bk36/beta,bfmgy,bfmgz,elbrth,eldeth
      common/bk46/nodm(8),reftem,thicknes
      common/bk48/stress(4),strain(4),d(4,4),ipt,n,nstate
      common/newz1/ibase,iadd5,iadd6,icnt6,locst6
      dimension xx(2,4),y(*),z(*),ix(4,*),betap(*),matp(*),
     1 den(*),nbck(*),matype(*),ym(4,*),thick(*)
c


        write(7777,*) '-- in2d.f'

c
      do 30 n=1,numelt
      nbck(n)=0
      mtype=matp(n)
      thicknes=thick(mtype)
      if (mtype.ne.0) model=matype(mtype)
      do 10 i=1,4
   10 nodm(i)=ix(i,n)
c
      do 20 i=1,4
      xx(1,i)=y(nodm(i))
   20 xx(2,i)=z(nodm(i))
      if (mtype.eq.0) go to 30
      beta  =betap(n)
      reftem=ym(1,n)
c
      call strint
c
c     compute and store lumped mass vector
c
      rho=den(mtype)
c
      call massl (ym(1,n),xx,rho,thicknes)
c
   30 continue
c&&&&&&&&&&&&&&&&&&
      if(ibase.eq.1)then
      call strind
      endif
c&&&&&&&&&&&&&&&&&&&&
      return
      end

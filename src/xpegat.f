      subroutine xpegat(iopt)
c     implicit double precision (a-h,o-z)                                    dp
       use mod_parameters
      common/bk01/h4(4,5),p14(4,5),p24(4,5)
      common/bk02/ioofc,iphase,imass,lpar(9)
      common/range/mft,mlt,lft,llt,nftm1
      common/intgrt/ipt
      common/xcom25/ipcmp(7,10),ipmcmp(10,10),ipecmp(10,10),
     1              xpcmp(2,10),npcmp,ipcm
      common/vect6/
     1 cg11(nelemg),cg21(nelemg),cg12(nelemg),cg22(nelemg),
     2 cg13(nelemg),cg23(nelemg),cg14(nelemg),cg24(nelemg),
     3 ch11(nelemg),ch21(nelemg),ch12(nelemg),ch22(nelemg),
     4 ch13(nelemg),ch23(nelemg),ch14(nelemg),ch24(nelemg),
     5 ci11(nelemg),ci21(nelemg),ci12(nelemg),ci22(nelemg),
     6 ci13(nelemg),ci23(nelemg),ci14(nelemg),ci24(nelemg)
      common/vect12/svvol(nelemg,4),
     1 spy1(nelemg,4),spy2(nelemg,4),spy3(nelemg,4),spy4(nelemg,4),
     2 spz1(nelemg,4),spz2(nelemg,4),spz3(nelemg,4),spz4(nelemg,4),
     3 sph1(nelemg,4),sph2(nelemg,4),sph3(nelemg,4),sph4(nelemg,4)
      common/vect91/xargat(nelemg),exx(nelemg,4),eps(nelemg,5),
     > sa(nelemg),
     1 sb(nelemg),sc(nelemg),ang1(nelemg),ang2(nelemg)
      if(npcmp.eq.0)return
      if(iopt.eq.1)then
      do 10 i=1,npcmp
      nortyp=ipcmp(4,i)
      if(nortyp.eq.1.or.nortyp.eq.2)then
      xpcmp(1,i)=0.
      xpcmp(2,i)=0.
      elseif(nortyp.eq.3)then
      xpcmp(1,i)=-1.e20
      elseif(nortyp.eq.4)then
      xpcmp(1,i)=1.e20
      endif
   10 continue
      elseif(iopt.lt.0)then
      do 30 i=1,npcmp
      nmat=ipcmp(5,i)
      if(nmat.ne.0)then
      do 40 j=1,nmat
      if(ipmcmp(j,i).eq.-iopt)then
      call xecomp(ipcmp(3,i),mft,mlt)
      nortyp=ipcmp(4,i)
      if(nortyp.eq.1.and.lpar(5).eq.0)then
      do 50 k=mft,mlt
      radius=h4(1,ipt)*ch11(k)+h4(2,ipt)*ch12(k)+
     1       h4(3,ipt)*ch13(k)+h4(4,ipt)*ch14(k)
      vol=radius*svvol(k,ipt)
      xpcmp(1,i)=xpcmp(1,i)+xargat(k)*vol
      xpcmp(2,i)=xpcmp(2,i)+vol
   50 continue
      elseif(nortyp.le.2)then
      do 55 k=mft,mlt
      xpcmp(1,i)=xpcmp(1,i)+xargat(k)*svvol(k,ipt)
      xpcmp(2,i)=xpcmp(2,i)+svvol(k,ipt)
   55 continue
      elseif(nortyp.eq.3)then
      do 60 k=mft,mlt
      if(xargat(k).gt.xpcmp(1,i))xpcmp(1,i)=xargat(k)
   60 continue
      elseif(nortyp.eq.4)then
      do 65 k=mft,mlt
      if(xargat(k).lt.xpcmp(1,i))xpcmp(1,i)=xargat(k)
   65 continue
      endif
      goto 30
      endif
   40 continue
      endif
   30 continue
      elseif(iopt.eq.3)then
      do 100 i=1,npcmp
      if(ipcmp(5,i).eq.0.and.ipcmp(6,i).eq.0)then
      call xecomp(ipcmp(3,i),lft,llt)
      nortyp=ipcmp(4,i)
      if(nortyp.eq.1.and.lpar(5).eq.0)then
      do 80 k=lft,llt
      radius=h4(1,ipt)*ch11(k)+h4(2,ipt)*ch12(k)+
     1       h4(3,ipt)*ch13(k)+h4(4,ipt)*ch14(k)
      vol=radius*svvol(k,ipt)
      xpcmp(1,i)=xpcmp(1,i)+xargat(k)*vol
      xpcmp(2,i)=xpcmp(2,i)+vol
   80 continue
      elseif(nortyp.le.2)then
      do 85 k=lft,llt
      xpcmp(1,i)=xpcmp(1,i)+xargat(k)*svvol(k,ipt)
      xpcmp(2,i)=xpcmp(2,i)+svvol(k,ipt)
   85 continue
      elseif(nortyp.eq.3)then
      do 90 k=lft,llt
      if(xargat(k).gt.xpcmp(1,i))xpcmp(1,i)=xargat(k)
   90 continue
      elseif(nortyp.eq.4)then
      do 95 k=lft,llt
      if(xargat(k).lt.xpcmp(1,i))xpcmp(1,i)=xargat(k)
   95 continue
      endif
      elseif(ipcmp(6,i).ne.0)then
      nele=ipcmp(6,i)
      do 110 j=1,nele
      ielem=ipecmp(j,i)-nftm1
      if(ielem.ge.lft.and.ielem.le.llt)then
      call xecomp(ipcmp(3,i),ielem,ielem)
      nortyp=ipcmp(4,i)
      if(nortyp.eq.1.and.lpar(5).eq.0)then
      radius=h4(1,ipt)*ch11(ielem)+h4(2,ipt)*ch12(ielem)+
     1       h4(3,ipt)*ch13(ielem)+h4(4,ipt)*ch14(ielem)
      vol=radius*svvol(ielem,ipt)
      xpcmp(1,i)=xpcmp(1,i)+xargat(ielem)*vol
      xpcmp(2,i)=xpcmp(2,i)+vol
      elseif(nortyp.le.2)then
      xpcmp(1,i)=xpcmp(1,i)+xargat(ielem)*svvol(ielem,ipt)
      xpcmp(2,i)=xpcmp(2,i)+svvol(ielem,ipt)
      elseif(nortyp.eq.3)then
      if(xargat(ielem).gt.xpcmp(1,i))xpcmp(1,i)=xargat(ielem)
      elseif(nortyp.eq.4)then
      if(xargat(ielem).lt.xpcmp(1,i))xpcmp(1,i)=xargat(ielem)
      endif
      endif
  110 continue
      endif
  100 continue

      elseif(iopt.eq.4)then
      do 20 i=1,npcmp
      nortyp=ipcmp(4,i)
      if(nortyp.eq.1.or.nortyp.eq.2)then
      xpcmp(1,i)=xpcmp(1,i)/xpcmp(2,i)
      endif
   20 continue
      endif
      return
      end

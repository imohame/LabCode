      subroutine engcal (energy,thick,xx,ener,maxint,ln,vol0)
c     implicit double precision (a-h,o-z)                                    dp
c
      common/bk02/ioofc,iphase,imass,lpar(9)
      common/cn3/ibbari,intgi,nmbfi,ithopi,ithcri,ithini,iengri,
     1           ijinti
      dimension energy(ln,*),thick(ln,*),xx(2,4),h(4)
      equivalence (lpar(5),ityp2d)
c
      ener=0.
      wgt =1.
      vol0=0.
      if (maxint.eq.1) wgt =4.
      if(iengri.eq.0)then
      do 60 lst=1,maxint
      ipt=lst
      if (maxint.eq.1) ipt=5
      call bass2r (h,det,xx,ipt)
      if (ityp2d.eq.0) go to 20
      radius=1.0
      go to 40
   20 radius=h(1)*xx(1,1)+h(2)*xx(1,2)+h(3)*xx(1,3)+h(4)*xx(1,4)
   40 fac=radius*det*wgt
      vol0=vol0+fac
      ener=ener+energy(1,lst)*fac
   60 continue
      elseif(iengri.eq.1)then
      do 70 lst=1,maxint
      ener=ener+wgt*.25*thick(1,lst)
   70 continue
      endif
c
      return
      end

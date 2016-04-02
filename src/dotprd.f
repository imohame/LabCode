      function dotprd (v1,v2)
c     implicit double precision (a-h,o-z)                                    dp
      common/bk03/numdc,imassn,idampn,irller,penstf
      common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
      dimension v1(*),v2(*)
      common/main_block/ a(1)
c
      dotprd=fdot(v1,v2,neq)
      if (numdc.eq.0) return
      k73=igtpnt(73)
      k57=igtpnt(57)
      k72=igtpnt(72)
      change=0.0
      call moddot(a(k57),a(k72),a(k73),numdc,v1,v2,change)
      dotprd=dotprd-change
      return
c
      end

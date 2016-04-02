      subroutine sslcs (a,model)

c     implicit double precision (a-h,o-z)                                    dp
c
      common/bk00/k01,k02,k03,k04,k05,k06,k07,k08,k09,k10,k11,k12
      common/bk02/ioofc,iphase,imass,lpar(9)
      common/bk06/nprnt,mprint,itmpop,numelt,jprint,idump,locstr
      common/bk14/lfna(15),lfnt(6)
      common/bk16/maxint,hgc
      common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
      common/bk18/nummat,ityp2d,ako(31)
      common/intgp/d(4,4),ipt,nel,nelsub
      dimension a(*)
c
      ln   = maxint*lpar(9)
      k08  = igtpnt(8)
      mtp  = locdbl(a(k08),nelsub-1)
      nm   = igtpnt(12)+48*(mtp-1)
      ii   = igtpnt(14)+ln*(nel-1)
      kk   = igtpnt(6) +nel-1
      nn   = ii+lpar(9)*(ipt-1)
      ne   = nn+lpar(9)-1
      nt   = ne-idump
      nu   = igtpnt(11)+nummat+mtp-1
      nray = igtpnt(11)+3*nummat+(mtp-1)*3
      mm   = igtpnt(2) +4*(nel-1)
      k63  = igtpnt(63)
      k64  = igtpnt(64)
      k81  = igtpnt(81)
      k82  = igtpnt(82)
      mtb  = igtpnt(12)+96*nummat
      mtb1 = locdbl(a(mtb),mtp-1)
      mtb2 = mtb+nummat+mtb1

      if(iphase.eq.2)then
         ktm = k81+1
      else
         ktm = k82+1
      endif
c
      if (model.eq.0) then
         call s0mn(a(nn),a(ne),a(nt),ln)
         return
      endif
c
ck----Elasticity :
      if(model.eq.1)then
      call s1mn (a(nm),a(nn),a(ne),a(nt),ln)
ck
ck----Orthotropic Elasticity :
      elseif(model.eq.2)then
      call s2mn (a(nm),a(nm+29),a(nn),a(ne),a(nt),ln)
ck
ck----Elastoplasticity :
      elseif(model.eq.3)then
      if (ityp2d.le.1) then
      call s3mn(a(nm),a(nm+29),a(nn),a(nn+4),a(ne),ln)
      else
      call s3mnp(a(nm),a(nm+29),a(nn),a(nn+4),a(ne),a(nt),ln)
      endif

c---------------------------------
c.....add up matslip model by Waeil Ashamwi
c---------------------------------

*-----Single Crystal (Double-Slip) :
       elseif (model.eq.4)then
          if (ityp2d.le.1) then
          call matpoly (a(k08),a(nn),a(k12),ln)
       endif
*
* ... Polycrystal with dislocation eveol.  (Double-Slip) :
       elseif(model.eq.5)then
          if (ityp2d.le.1) then
*          call matpoly_dislocation_2s (a(k08),a(nn),a(k12),ln)
       endif
*
* ... Polycrystal with dislocation eveol.  (4-Slip) :
       elseif(model.eq.6)then
          if (ityp2d.le.1) then
*          call matpoly_dislocation_4s (a(k08),a(nn),a(k12),ln)
       endif

*
*-----Bicrystal (Double-slip) :
       elseif(model.eq.7) then
          if (ityp2d.le.1) then
*          call matslipbxds (a(nn),a(k12),ln)
       endif
*
* ... Polycrystal  (Double-Slip) :
       elseif(model.eq.8)then
          if (ityp2d.le.1) then
*          call matpoly_ds (a(k08),a(nn),a(k12),ln)
       endif



c.... no material specified

      else
         write(lfnt(2),250)model
         call bye(2)
      endif

c
c.... stiffness proportional damping
      if(imass.eq.1)then
      call raydmp(a(nray))
      endif
c

  231 format(' **fatal error** '//
     1       ' model',i5,'is not implemented for plane stress'//)
  250 format(' **fatal error** mat model',i5,' not coded (sslcs)')


      return
      end

      subroutine fissln(ist,b,a,s)
c     implicit double precision (a-h,o-z)                                    dp
      common/cn4/itypei,ianali,ibwmni,icorei,ipcori,isolti,gamai,
     1           betai,raydmi
      common/double/iprec,ncpw,unit
      dimension b(1),a(1),s(1)
c
c.... solve with fissle
      if(icorei.ge.0)then
         call fissla(ist,b,a,s)
c
c.... solve with actol
      else
         call fisslb(ist,b,a,s)
      endif
c
      return
      end

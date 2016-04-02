      subroutine strset
c     implicit double precision (a-h,o-z)                                    dp
      common/bk14/lfna(15),lfnt(6)
      common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
      common/main_block/ a(1)
c
      k13=igtpnt(13)
      k16=igtpnt(16)
      call getnam  (lfna(11),namef)
      write (*,10) namef
      call rdabsf (lfna(11),lendr,1,0,ioerr)
      call riosta (lfna(11))
      iadd=lendr+k16
      call rdabsf (lfna(11),a(k13),nwebuf,iadd,ioerr)
      call riosta  (lfna(11))
      call intemp
      call blkcpy (a(k13),a(k16),nwebuf)
      return
c
   10 format(//'  file',2x,a8,2x,'opened for stress initialization')
      end

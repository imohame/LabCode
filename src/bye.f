      subroutine bye (n)
c     implicit double precision (a-h,o-z)                                    dp
c
      common/bk14/lfna(15),lfnt(6)
      common/bk15/cpuio(3,12),cpuip(3,12)
      common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
      common/bk29/numfrq,clengt
      common/fissn0/maxneq,mwspac,ntpe0,ntpe1,nfissl
      common/total/itrlas,irflas,irhlas,itrtot,irftot,irhtot
	  integer:: n
      common/iobufs/iobuf(700,13)                                       vax750
      character*8 keep                                                  vax750
      keep='keep'                                                       vax750
c

       call prtout
	   write(*,*) '==>>> bye.f'
       write(*,*) 'In routine bye, n = ',  n

      call rwabsf(iobuf(1,9) ,keep,0,0,0)             
      call rwabsf(iobuf(1,10),keep,0,0,0)            
      call rwabsf(iobuf(1,11),keep,0,0,0)
                      
c

cw      call timin (cpuio,cpuip,5,3)
cw      call timin (cpuio,cpuip,11,3)

c
c     calculate and print solution time log
c
      call header
      write(lfnt(2),30) ((cpuio(i,j),i=1,3),j=1,11)
      tctm=cpuio(1,11)+cpuio(2,11)+cpuio(3,11)
c
      write(lfnt(2),60) tctm,cpuio(1,11)
c*waeil*      write (*,60) tctm,cpuio(1,11)
c*waeil*      write (*,120) itrtot,irftot,irhtot
      write(lfnt(2), 120) itrtot,irftot,irhtot
c
      if (n.ne.1) go to 10
      write(lfnt(2),40)
c*waeil*      write (*,40)
      go to 20
   10 if (n.ne.2) go to 20
      write(lfnt(2),50)
c*waeil*      write (*,50)
c
   20 continue
 
	  call exita(n)
c       call exita (cdim(n,1))
c
c
   30 format(///1x,'t i m i n g   i n f o r m a t i o n ',/
     1  58x,'cpu',7x,'sys',7x,'i/o'
     2///5x,'i n p u t   p h a s e                             ',3f10.3
     3 //5x,'     portion due to bandwidth minimizer           ',3f10.3
     4///5x,'i n i t i a l i z a t i o n   p h a s e           ',3f10.3
     5 //5x,'     portion due to fissle                        ',3f10.3
     6///5x,'s o l u t i o n   p h a s e                       ',3f10.3
     7 //5x,'     portion due to fissle excluding eq. iter.    ',3f10.3
     8 //5x,'     portion due to equilibrium iterations        ',3f10.3
     9  /5x,'          portion due to fissle in eq. iter       ',3f10.3
     $ //5x,'     high speed printer output                    ',3f10.3
     $ //5x,'     plotfile generation                          ',3f10.3
     $///5x,'t o t a l    s o l u t i o n   t i m e            ',3f10.3)
   40 format(///1x,'n o r m a l   t e r m i n a t i o n',//)
   50 format(///1x,'e r r o r    t e r m i n a t i o n ',//)
   60 format(/' total time charged for CPU+IO =',e14.5,
     1       /' total time charged for CPU    =',e14.5,/)
  120 format('     total iterations             = ',i5,/,
     1       '     total stiffness reformations = ',i5,/,
     2       '     total rhs evaluations        = ',i5//)
      end

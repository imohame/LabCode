       subroutine addres(jdiag,mht,nwk,ma,mb)
c     implicit double precision (a-h,o-z)                                    dp
c
       common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
       dimension jdiag(*),mht(*)
c
       jdiag(1)=1
       ma=0
       mb=0
       do 10 i=2,neq
       ma=max(ma,mht(i))
       mb=mb+mht(i)
   10 jdiag(i)=jdiag(i-1)+1+mht(i)
       ma=ma+1
       mb=mb/neq
       nwk=jdiag(neq)
c
       return
       end

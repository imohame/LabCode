       subroutine fwdbak(b,jdiag,a,neq)
c     implicit double precision (a-h,o-z)                                    dp
c
       dimension a(*),b(*),jdiag(*)
       jr=0
       do 20 j=1,neq
       jd=jdiag(j)
       jh=jd-jr
       is=j-jh+2
       if(jh-2) 20,10,10
   10 b(j)=b(j)-fdot(a(jr+1),b(is-1),jh-1)
   20 jr=jd
       do 30 i=1,neq
   30 b(i)=b(i)/a(jdiag(i))
       j=neq
       jd=jdiag(j)
   40 d=b(j)
       j=j-1
       if(j.le.0) return
       jr=jdiag(j)
       if(jd-jr.le.1) go to 60
       is=j-jd+jr+2
       k=jr-is+1
       do 50 i=is,j
   50 b(i)=b(i)-a(i+k)*d
   60 jd=jr
       go to 40
       end

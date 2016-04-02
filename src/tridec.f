       subroutine tridec(b,jdiag,a,neq)
c     implicit double precision (a-h,o-z)                                    dp
c
       dimension a(*),b(*),jdiag(*)
       jr=0
       do 50 j=1,neq
       jd=jdiag(j)
       jh=jd-jr
       is=j-jh+2
       if(jh-2) 50,30,10
   10 ie=j-1
       k=jr+2
       id=jdiag(is-1)
       do 20 i=is,ie
       ir=id
       id=jdiag(i)
       ih=min(id-ir-1,i-is+1)
       if(ih.gt.0) a(k)=a(k)-fdot(a(k-ih),a(id-ih),ih)
   20 k=k+1
   30 ir=jr+1
       ie=jd-1
       k=j-jd
       do 40 i=ir,ie
       d=a(i)
       a(i)=a(i)/a(jdiag(k+i))
       a(jd)=a(jd)-d*a(i)
   40 continue
   50 jr=jd
       return
       end

      subroutine clht (mht,lm,n,iprec)
c     implicit double precision (a-h,o-z)                                    dp
c
       use mod_parameters
      common/vect16/
     1 ls(nelemg),u12(nelemg),u14(nelemg),u21(nelemg),
     > u22(nelemg),u24(nelemg),
     2 u41(nelemg),u42(nelemg),u44(nelemg),q1(nelemg),
     > q2(nelemg),q3(nelemg),q4(nelemg),
     3 t11(nelemg),t12(nelemg),t14(nelemg),t41(nelemg),
     > t44(nelemg)
      dimension lm(44*iprec,*),mht(*)
c
      do 10 i=1,n
      ls(i)=100000
      ls(i)=ls(i)*(1-sign(1,lm(1,i)-1))/2+min(ls(i),lm(1,i))
      ls(i)=ls(i)*(1-sign(1,lm(2,i)-1))/2+min(ls(i),lm(2,i))
      ls(i)=ls(i)*(1-sign(1,lm(3,i)-1))/2+min(ls(i),lm(3,i))
      ls(i)=ls(i)*(1-sign(1,lm(4,i)-1))/2+min(ls(i),lm(4,i))
      ls(i)=ls(i)*(1-sign(1,lm(5,i)-1))/2+min(ls(i),lm(5,i))
      ls(i)=ls(i)*(1-sign(1,lm(6,i)-1))/2+min(ls(i),lm(6,i))
      ls(i)=ls(i)*(1-sign(1,lm(7,i)-1))/2+min(ls(i),lm(7,i))
      ls(i)=ls(i)*(1-sign(1,lm(8,i)-1))/2+min(ls(i),lm(8,i))
   10 continue
      do 20 i=1,n
      if (lm(1,i).ne.0) mht(lm(1,i))=max(mht(lm(1,i)),lm(1,i)-ls(i))
      if (lm(2,i).ne.0) mht(lm(2,i))=max(mht(lm(2,i)),lm(2,i)-ls(i))
      if (lm(3,i).ne.0) mht(lm(3,i))=max(mht(lm(3,i)),lm(3,i)-ls(i))
      if (lm(4,i).ne.0) mht(lm(4,i))=max(mht(lm(4,i)),lm(4,i)-ls(i))
      if (lm(5,i).ne.0) mht(lm(5,i))=max(mht(lm(5,i)),lm(5,i)-ls(i))
      if (lm(6,i).ne.0) mht(lm(6,i))=max(mht(lm(6,i)),lm(6,i)-ls(i))
      if (lm(7,i).ne.0) mht(lm(7,i))=max(mht(lm(7,i)),lm(7,i)-ls(i))
      if (lm(8,i).ne.0) mht(lm(8,i))=max(mht(lm(8,i)),lm(8,i)-ls(i))
   20 continue
c
      return
      end

      subroutine prtstr (ipst,matp,matype)
c     implicit double precision (a-h,o-z)                                    dp
c
c.... print stresses in result file and write stresses to tecplot
c

* [ P A R A M E T E R S]
* ......................
      integer, parameter :: print_intv = 10000
      

      common/bk02/  ioofc,iphase,imass,model,numel,lpar(7)
      common/bk05/  ifil,iadd,maxsiz,head(12)
      common/bk06/  nprnt,mprint,itmpop,numelt,jprint,idump,locstr
      common/bk08/  kprint,nstep,ite,ilimit,newstf
      common/bk14/  lfna(15),lfnt(6)
      common/bk18/  nummat,ityp2d,ako(31)
      common/bk20/  ntotal
      common/bk32/  nsref,nequit,time,timep,lprint,nprint
      common/bk34/  stf(1)
      common/bk48/  stress(4),strain(4),d(4,4),ipt,n,nstate
      common/newz1/ ibase,iadd5,iadd6,icnt6,locst6
      common/poutt/ nst(print_intv)
      common/main_block/  a(1)

      common /wblock1/ iplotDirxy, x_area, yield_stress


      dimension ipst(*),matp(*),matype(*),
     1 ast(4,26),last(4)
      character*8 astate(3)
      character*9 fname(8)
      integer :: ipr

      data astate/'elastic','plastic','failure'/,
     1fname/'sig-yy','sig-zz','sig-hoop','sig-yz','max-sig',
     2 'min-sig','angle','yield fn'/
cw
c      numel=numelt
c
      if(ityp2d.eq.1)then
        fname(3)='sig-xx'
      endif
c
c.... loop over elements
      do 70 n = 1,numel
         mtype = matp(n)
         if (mtype.eq.0) go to 60
c
c.... find model
         model = matype(mtype)
         if (mprint.le.0) go to 10
         if (ipst(n).eq.1) go to 70
c
c.... loop over integration points
   10    continue
         do 50 lst=1,4
            ipt = lst
c.... extract stresses
            call getstr(a)
            if (mprint.le.0) go to 40
            if (nprnt.gt.0) go to 20
c.... print header if necessary
            nprnt = 40
!!!!!            call header
            write(lfnt(2),80) nstep,timep
   20       nprnt = nprnt - 1

c.... compute quantities other than stresses
            call prinst
            ast(lst,1) = stress(1)
            ast(lst,2) = stress(2)
            ast(lst,3) = stress(3)
            ast(lst,4) = stress(4)
            ast(lst,5) = strain(2)
            ast(lst,6) = strain(3)
            ast(lst,7) = strain(4)
            ast(lst,8) = strain(1)
            last(lst)  = nstate

c.... print results
            if (ipt.lt.4) go to 50
c            write(lfnt(2),200)n
c            write(lfnt(2),210)mtype,(astate(last(i)),i=1,4)
c            do 220 i=1,8
c  220       write(lfnt(2),230) fname(i),(ast(j,i),j=1,4)

*  ____________________________________________________________________
* |                                                                    |
* |... Print out stresses for each element                             |
* |...      ast(1,1) : stress_xx                                       |
* |...      ast(1,2) : stress_xy                                       |
* |...      ast(1,4) : stress_yy                                       |
* |____________________________________________________________________|
*

c      do ipr = 1, print_intv   
         if (mod(nstep,(nst(2)-nst(1))) == 0) then
            if (nstep .ne. 0) then
			write(500,235) ast(1,1)/yield_stress
            write(501,235) ast(1,2)/yield_stress
            write(502,235) ast(1,4)/yield_stress
            write(503,235) ast(1,3)/yield_stress
            exit
			endif
         endif
c      end do 
         
           go to 50
c
   40      kloc = 5*(lst-1)
   50    continue
   60    locstr = locstr+20
   70 continue
c
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      if(kprint.le.0) return
      if (ibase.ge.1) then
         lmark  = icnt6 - 2
         ipoint = igtpnt(14)
         call wrabsg(lfna(11),a(ipoint),lmark,iadd6,stf(ntotal+1),0)
!!!!!!!!!!!!!         call riosta(lfna(11))
         iadd6 = iadd6 + lmark
      endif
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
      return

   80 format(/' element stress results for step ',i5,
     1 ' at time',3x,e14.4)
c  200 format(/,' element',i5,' int pt:',t25,'1',t40,'2',t55,
c     1 '3',t70,'4')
c  210 format(' material',i4,' state:',t22,a7,t37,a7,t52,a7,
c     1 t67,a7)
  223 format(t19,e11.3,t34,e11.3,t49,e11.3,t64,e11.3)
c  230 format(3x,a9,t19,e11.3,t34,e11.3,t49,e11.3,t64,e11.3)

  235 format(5x,e20.10)

      end

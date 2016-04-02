      subroutine odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad)

      parameter ( maxstp = 100,
     >            nmax   = 24,
     >            two    = 2.0,
     >            zero   = 0.00,
     >            tiny   = 1.0e-20 )
      common /path/ kmax,kount,dxsav,xp(200),yp(10,200)
      dimension ystart(nmax),yscal(nmax),y(nmax),dydx(nmax)

C     Adaptive step method based on fifth order Runge-Kutta
C     Switch is done to Implicit method if problem is stiff
C     Based on algorithm developed by M.A.  Zikry
c
c
       x    = x1
       h    = sign(h1,x2-x1)
       nok  = 0
       nbad = 0
       kmax = 0
	   !write(*,*) 'beginning of odeint'
	   !write(*,*) nvar
       do i = 1,nvar
          y(i)=ystart(i)
       end do

* [Assure storage of first step]
* ..............................
           xsav=x-dxsav*two

* [Take at most maxstp steps]
* ...........................
       do 16 nstp=1,maxstp
          call derivs(x,y,dydx,nvar)
		  !write(*,*) 'in odeint after derivs'
		  !write(*,*) nvar

* [Scaling used to monitor accuracy. this can be modified as needed]
* ..................................................................
          do 12 i = 1, nvar
             yscal(i)=abs(y(i))+abs(h*dydx(i))+tiny
12        continue

          if((x+h-x2)*(x+h-x1) .gt. zero) h = x2-x

* [If step can overshoot end, cut down stepsize]
* ..............................................
          call RKQC(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,hmin)
          !write(*,*) 'in odeint after rkqc'
		  !write(*,*) nvar
          if ((x-x2)*(x2-x1) .ge. zero) then

* [Are we done?]
* ..............
             do 14 i=1,nvar
                ystart(i)=y(i)
14           continue
             return
          endif

16     continue
       if(abs(hnext) .lt. hmin) h = hnext

c      write(27,31)
31     format(' stepsize smaller than minimum')

       return
       end


      subroutine printm (n,den,thick,temmat,raystf,prop,prp48,
     >                   headng)


c     implicit double precision (a-h,o-z)                                    dp
c
c.... print out material properties and define other material
c     constants (e.g. elastic matrix coeffients)
c

      common /bk02/    ioofc,iphase,imass,model,lpar(8)
      common /bk14/    lfna(15),lfnt(6)
      common /meichu/  modeltype


      character*6 headng
      dimension den(*),prop(*),prp48(*),headng(12)


      write (lfnt(2),160) n,model,den(1),thick,temmat,
     >                    raystf

      if(den(1).eq.0.0) den(1)=1.
      modeltype = model
c--------------------------------
c.....Elastic-Isotropic Materials
c--------------------------------
      if(model.eq.1) then
      write (lfnt(2),170) (prop(i),i=1,3)
      call setse1 (prop)

c----------------------------------
c.....Elastic-Orthotropic Materials
c----------------------------------
      elseif(model.eq.2)then
      write (lfnt(2),180) (prop(i),i=1,11)
      call setse2 (prop,prop(30))

c----------------------------------------
c.....Elastoplastic Materials (Von Mises)
c----------------------------------------
      elseif(model.eq.3)then
      write (lfnt(2),190) (prop(i),i=1,21)
      prop(30)=prop(1)
      prop(31)=prop(2)
      call setse1 (prop(30))

c---------------------------------
c.....Single Crystal (Double-Slip)
c---------------------------------
      elseif(model.eq.4)then
*         write (lfnt(2),200) (prop(i),i=1,9)
         write (lfnt(2),210) (prop(i),i=1,4), 
     >                       (prop(i), i=7,10)
         write (lfnt(2),300) n

         prop(10)=prop(1)
         prop(11)=prop(2)
         call setse1 (prop(10))

c------------------------------------------------
c.....Polycrystal (G.B.) --- Double-Slip
c------------------------------------------------
      elseif(model.eq.5)then
      write (lfnt(2),200) (prop(i),i=1,9)
      write (lfnt(2),300) n
      prop(10)=prop(1)
      prop(11)=prop(2)
      call setse1 (prop(10))

c----------------------------------------------
c.....Bicrystal (Sigma9 G.B.) --- Multiple-Slip
c----------------------------------------------
      elseif(model.eq.6)then
         write (lfnt(2),210) (prop(i),i=1,4), 
     >                       (prop(i), i=7,9)
*        write (lfnt(2),220) (slip_angle(i,n), i=1,4) 
         write (lfnt(2),300) n

         prop(10)=prop(1)
         prop(11)=prop(2)
         call setse1 (prop(10))

c-----------------------------------------------
c.....Bicrystal ---  Double-Slip
c-----------------------------------------------
      ELSEIF(model.eq.7) THEN
         WRITE (lfnt(2),210) (prop(i),i=1,8)
         prop(10) = prop(1)
         prop(11) = prop(2)
         CALL setse1 (prop(10))

c------------------------------------------------
c.....Polycrystal --- Double-Slip
c------------------------------------------------
      elseif(model.eq.8)then
         write (lfnt(2),200) (prop(i),i=1,9)
         write (lfnt(2),300) n
         prop(10)=prop(1)
         prop(11)=prop(2)
         call setse1 (prop(10))
c
      else
         write(lfnt(2),355)model
         call bye(2)
      endif
c
      return

  160 format(//' material constants set number .... ',i5,
     1        4x,' material model .... ',i5//
     2 5x,'den .............................. =', e12.4/
     3 5x,'thickness  (plane stress) ........ =', e12.4/
     4 5x,'material reference temperature ... =', e12.4/
     5 5x,'stiffness prop Rayleigh dmp coeff. =', e12.4/)
  170 format(
     1 5x,'e ................................ =', e12.4/
     2 5x,'vnu .............................. =', e12.4/
     3 5x,'ipc .............................. =', e12.4)
  180 format(
     1 5x,'e(a) ............................. =', e12.4/
     2 5x,'e(b) ............................. =', e12.4/
     3 5x,'e(c) ............................. =', e12.4/
     4 5x,'vnu(ab) .......................... =', e12.4/
     5 5x,'vnu(ac) .......................... =', e12.4/
     6 5x,'vnu(bc) .......................... =', e12.4/
     7 5x,'g(ab) ............................ =', e12.4/
     8 5x,'material axes option ............. =', e12.4/
     9 5x,'yc    (option =1.0) ............. =', e12.4/
     $ 5x,'zc    (option =1.0) ............. =', e12.4/
     $ 5x,'gamma (option =2.0) ............. =', e12.4)
  190 format(
     1 5x,'e ................................ =', e12.4/
     2 5x,'vnu .............................. =', e12.4/
     3 5x,'yield ............................ =', e12.4/
     4 5x,'e (harden) ....................... =', e12.4/
     5 5x,'hardening parmeter ............... =', e12.4/
     6 5x,'effective plastic strain ......... =',8(e9.2,1x)/
     7 5x,'effective stress ................. =',8(e9.2,1x))

  200 format(
     1 5x,'e (Youngs modulus)....................... =', e12.4/
     2 5x,'vnu (Poissons ratio)..................... =', e12.4/
     3 5x,'tau y (Static yield stress).............. =', e12.4/
     4 5x,'1/m (Mater. rate sensitivity parameter).. =', e12.4/
     5 5x,'phi1 (Initial angle of slip dir. (1)).... =', e12.4/
     6 5x,'phi2 (Initial angle of slip dir. (2)).... =', e12.4/
     7 5x,'rr (Reference shear strain rate)......... =', e12.4/
     8 5x,'rcr (Critical shear strain rate)......... =', e12.4/
     9 5x,'ipc (Interval of printing for contours).. =', e12.4)

  210 format(
     1 5x,'e (Youngs modulus)....................... =', e12.4/
     2 5x,'vnu (Poissons ratio)..................... =', e12.4/
     3 5x,'tau y (Static yield stress).............. =', e12.4/
     4 5x,'1/m (Mater. rate sensitivity parameter).. =', e12.4/
     5 5x,'rr (Reference shear strain rate)......... =', e12.4/
     6 5x,'rcr (Critical shear strain rate)......... =', e12.4/
     7 5x,'ipc (Interval of printing for contours).. =', e12.4)

 220  format(
     > 5x,'phi1 (Initial angle of slip dir. (1)).... =', e12.4/ 
     > 5x,'phi2 (Initial angle of slip dir. (2)).... =', e12.4/
     > 5x,'phi3 (Initial angle of slip dir. (3)).... =', e12.4/
     > 5x,'phi4 (Initial angle of slip dir. (4)).... =', e12.4)

 300  format(
     > 5x,'Grain Number (in the polycrystalline).... =', i7)


  355 format(' **fatal error** printout of properties and/or',/,
     1       ' setup of coeff. for mat',i5,' not coded (printm)')


      end










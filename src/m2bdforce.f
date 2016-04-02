
      subroutine m2bdforce(ix,nodes,idir,u,nn)

      parameter (nodet = 40000)

      common/bk03/numdc,imassn,idampn,irller,penstf
      common/bk06/nprnt,mprint,itmpop,numelt,jprint,idump,locstr
      common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
      common/forces/fff(8,nodet)
      common/pvari/iaa,ibb
      common/wblock1/ iplot, x_area, yield_stress
 
      dimension sf(2,nodet),ix(4,*),nodes(*),idir(*),u(*)

      iaa = iaa + 1
      sf=0.0
!!!!      do 99 ii = 1, numnp
!!!!        sf(1,ii) = 0.
!!!!        sf(2,ii) = 0.
!!!!   99 continue
c
      dispp = 0.
      do 10 ii = 1, numdc
        ij    = idir(ii)
        dispp = dispp + u(ij)
   10 continue

      dispp = dispp/numdc
c ------------------------------------- 
      do 333 ii = 1 ,numdc
      ndd = nodes(ii)
      do 222 jj = 1, numelt
        do 111 kk = 1, 4
           nod = ix(kk,jj)
           if (nod.eq.ndd) then
           sf(1,ndd) = sf(1,ndd) + fff(2*kk-1,jj)
           sf(2,ndd) = sf(2,ndd) + fff(2*kk,jj)
           endif
  111   continue
  222 continue
  333 continue
c----------------------------------------
      tf1 = 0.
      tf2 = 0.
      do 666 ii=1,numdc
        ndd = nodes(ii)
        tf1 = tf1+sf(1,ndd)
        tf2 = tf2+sf(2,ndd)
  666 continue 
c
      tf3 = tf1/x_area
      tf4 = tf2/x_area

*     divide tf4 by the yield stress (.001 for non-dim.)
      tf5 = tf3/(yield_stress)
      tf6 = tf4/(yield_stress)
      
      if(iplot.eq.1) then 
* >>  {
        if(iaa.eq.ibb) then
!!           if(nn.ne.51) then
              write(126,*) dispp,tf5
              flush(126)
!!           endif
           ibb = iaa+100
        endif
* <<  }
      elseif(iplot.eq.2) then
* >>  {
        if(iaa.eq.ibb) then
!!           if(nn.ne.51) then
              write(126,*) dispp,tf6
              flush(126)
!!           endif
           ibb = iaa+100
        endif
* <<  { 
      endif

c
* ----------- F O R M A T  S T A T E M E N T S --------------
* -----------------------------------------------------------
  557   format(5x,f12.3,3x,f20.6)
c  
      return
      end










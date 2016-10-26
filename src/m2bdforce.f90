
    subroutine m2bdforce(ix,nodes,idir,u,nn)

    parameter (nodet = 40000)
    common/bk03/numdc,imassn,idampn,irller,penstf
    common/bk06/nprnt,mprint,itmpop,numelt,jprint,idump,locstr
    common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
    common/forces/fff(8,nodet)
    common/pvari/iaa,ibb
    common/wblock1/ iplotDirxy, x_area, yield_stress

    dimension sf(2,nodet),ix(4,*),nodes(*),idir(*),u(*)

    integer iaa,ndd,ii,jj,ij,nod
    real sf,dispp,ForceX,ForceY,StressX,StressY,StressXNor,StressYNor
    
    real x_area, yield_stress,u,fff
    integer numdc,ix,nodes,idir,ibb,iplotDirxy
    
    iaa = iaa + 1
    !-- only go through if iaa == ibb
    if(iaa.ne.ibb) then
        return
    endif
    
    sf=0.0
!!!!      do 99 ii = 1, numnp
!!!!        sf(1,ii) = 0.
!!!!        sf(2,ii) = 0.
!!!!   99 continue

    dispp = 0.
    do  ii = 1, numdc
        ij    = idir(ii)
        dispp = dispp + u(ij)
    enddo 
    !- average displ
    dispp = dispp/numdc
    !-- strain
    dispp = dispp/x_area
!!!c ------------------------------------- 
    do  ii = 1 ,numdc
        ndd = nodes(ii)
        do  jj = 1, numelt
            do  kk = 1, 4
               nod = ix(kk,jj)
               if (nod.eq.ndd) then
                sf(1,ndd) = sf(1,ndd) + fff(2*kk-1,jj)
                sf(2,ndd) = sf(2,ndd) + fff(2*kk,jj)
               endif
            enddo
        enddo
    enddo
!!!c----------------------------------------
    ForceX = 0.
    ForceY = 0.
    do  ii=1,numdc
        ndd = nodes(ii)
        ForceX = ForceX+sf(1,ndd)
        ForceY = ForceY+sf(2,ndd)
    enddo 

    StressX = ForceX/x_area
    StressY = ForceY/x_area

!!!!*     divide StressY by the yield stress (.001 for non-dim.)
    StressXNor = StressX/(yield_stress)
    StressYNor = StressY/(yield_stress)

    if(iplotDirxy.eq.1) then !x-dir
        if(iaa.eq.ibb) then
            write(126,*) dispp,StressXNor
            flush(126)
            ibb = iaa+100
        endif
    elseif(iplotDirxy.eq.2) then !y-dir
        if(iaa.eq.ibb) then
            write(126,*) dispp,StressYNor
            flush(126)
            ibb = iaa+100
        endif
    endif

!!!* ----------- F O R M A T  S T A T E M E N T S --------------
!!!* -----------------------------------------------------------
!!!  557   format(5x,f12.3,3x,f20.6)
      return
      end












  subroutine LineSegIntersection(p,p2,q,q2,pRatio,qRatio,rPtc,bIntersect)
    !-- source http://www.codeproject.com/Tips/862988/Find-the-Intersection-Point-of-Two-Line-Segments
    use mod_parameters
    implicit none
    real*8 , intent(in)::p(2),p2(2),q(2),q2(2) !-input lines p-p2, q-q2
    real*8 , intent(out)::rPtc(2),pRatio,qRatio !- output intersection point
    integer ,intent(out)::bIntersect !- =0 for no intersection
    integer i
    real*8 r(2),s(2),qp(2),rxs,qpxr,qpxs,CorssProduct2d

    bIntersect=0
    r=p2-p
    s=q2-q
    qp=q-p
    rxs=CorssProduct2d(r,s)
    qpxr=CorssProduct2d(qp,r)
    !- If r x s = 0 and (q - p) x r = 0, then the two lines are collinear.
    if ((abs(rxs) < TolZero16) .and. (abs(qpxr) <TolZero16 )) then
      return
    endif
    !- If r x s = 0 and (q - p) x r != 0, then the two lines are parallel and non-intersecting.
    if ((abs(rxs) < TolZero16) .and. (qpxr /=TolZero16 )) then
      return
    endif

    qpxs=CorssProduct2d(qp,s)
    !- pRatio = (q - p) x s / (r x s)
    pRatio = qpxs/rxs

    !- qRatio = (q - p) x r / (r x s)
    qRatio = qpxr/rxs

    !- If r x s != 0 and 0 <= t <= 1 and 0 <= u <= 1
    !- the two line segments meet at the point p + t r = q + u s.
    if ((0 <= pRatio .and. pRatio <= 1) .and. (0 <= qRatio .and. qRatio <= 1)) then
        !- We can calculate the intersection point using either t or u.
        rPtc = p + pRatio*r
        !-An intersection was found.
        bIntersect=1
    endif

    !- Otherwise, the two line segments are not parallel but do not intersect.
  end subroutine LineSegIntersection
!##############################################################################
!##############################################################################
  function CorssProduct2d(v1,v2)result(v1xv2)
    implicit none
    real*8 , intent(in)::v1(2),v2(2) !-input lines p-p2, q-q2
    real*8 ::v1xv2 !- output intersection point
    v1xv2=v1(1)*v2(2)-v1(2)*v2(1) != v1x*v2y-v1y*v2x

  end function CorssProduct2d
  !##############################################################################
  !##############################################################################
    subroutine GetLineStartEndPts(rPt,Vnx,Vny,rDelta,p1,p2)
      implicit none
      real*8 , intent(in)::Vnx,Vny,rDelta,rPt(2) !- input line direction
      real*8 , intent(out)::p1(2),p2(2) !-output points
      real*8 ::v1xv2 !- output intersection point
      p1(1)=rPt(1)-rDelta*Vnx
      p1(2)=rPt(2)-rDelta*Vny

      p2(1)=rPt(1)+rDelta*Vnx
      p2(2)=rPt(2)+rDelta*Vny

    end subroutine GetLineStartEndPts
    !##############################################################################
    !##############################################################################
    function CalcPolygonArea(xs,ys,nPts)result(rArea)
        implicit none
        real*8 , intent(in)::xs(nPts),ys(nPts) !-input lines p-p2, q-q2
        real*8 ::rArea !- output area
        real*8 ::v1(2),v2(2),CorssProduct2d !- output area
        integer i,nPts

        rArea=0.0
        do i=2,nPts
            v1(1)=xs(i+1)-xs(i)
            v1(2)=ys(i+1)-ys(i)

            v2(1)=xs(i)-xs(i-1)
            v2(2)=ys(i)-ys(i-1)
            rArea=rArea+CorssProduct2d(v1,v2)
        end do
        v1(1)=xs(1)-xs(4)
        v1(2)=ys(1)-ys(4)

        v2(1)=xs(4)-xs(3)
        v2(2)=ys(4)-ys(3)
        rArea=rArea+CorssProduct2d(v1,v2)

        rArea=0.5*abs(rArea)

    end function CalcPolygonArea
    !##############################################################################
    !##############################################################################
    subroutine GetPtOnLine(rP1,rP2,DistRatio,rPtc)
        implicit none
        real*8 , intent(in)::rP1(2),rP2(2),DistRatio !- input line direction
        real*8 , intent(out)::rPtc(2) !-output points
        real*8 ::v21(2) !- output intersection point

        v21=rP2-rP1
        rPtc=rP1+DistRatio*v21

    end subroutine GetPtOnLine
    !##############################################################################
    !##############################################################################

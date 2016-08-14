subroutine print_str(iele, ink, matp)

    use mod_parameters
    common/custr/ sign1(nelemg), sign2(nelemg), sign3(nelemg), sign4(nelemg)
    common /wblock1/ iplotDirxy, x_area, yield_stress
    common/wblock8/ abc(573, nume, 4), his(573, nume, 4)
    common /wblock10/ng, grain_mo(1000,3),bv(no_mat, 87), nssmat(1000), nssimmat(1000) !changed to accomodate multiple bv lengths WML
    common/slipplane110/ nsp110(6)

    real sig(4), cleave(3, 3), cleave110(6, 3)
    real str, str_cleavage, normc, yield_stress
    integer iele, ink, iHE, jHE, nsp110, nsp110fcc(6)
    real str110, strmax110
!    real vectn(2)

    dimension matp(*)

    data nsp110fcc/1, 2, 3, 4, 6, 8/
    
!!--------------------------------------debugging
!      write(969, *) 'print_str'
!      do i = 1, 400
!          matid= matp(i)!a(k08+i)
!      write(969, *) i,matid,nssmat(matid),nssimmat(matid)
!      enddo
!      stop
!!--------------------------------------debugging

    str = 0.0
    str_cleavage = 0.0

    str110 = 0.0
    strmax110 = 0.0

    sig(1) = sign1(iele)
    sig(2) = sign2(iele)
    sig(3) = sign3(iele)
    sig(4) = sign4(iele)

    !c     cleavage planes {100}
    do i = 1, 3
        cleave(1, i) = abc(362 + i, ink, 1)
        cleave(2, i) = abc(365 + i, ink, 1)
        cleave(3, i) = abc(368 + i, ink, 1)
    end do

    do i = 1, 3
        str = sig(1) * cleave(i, 2)**2.0 + sig(2) * cleave(i, 3)**2.0 +sig(4) * cleave(i, 2) * cleave(i, 3) * 2.0
        if (abs(str) > str_cleavage) then
            str_cleavage = abs(str)
        end if
    end do
    
    write(935, 1001) str_cleavage/yield_stress !	  'stress001.out'

    !c     hydrogen embrittlement planes {110} for bcc
!    write(969, *) '------------------------------',ink,iele,matp(ink),nssmat(matp(ink))
    if (nssmat(matp(ink)) == 24) then
        do iHE = 1, 6
            jHE = nsp110(iHE)
            cleave110(iHE, 1) = abc(105 + jHE, ink, 1)
            cleave110(iHE, 2) = abc(129 + jHE, ink, 1)
            cleave110(iHE, 3) = abc(153 + jHE, ink, 1)
!            write(969, *) '24', cleave110(iHE, 1), cleave110(iHE, 2), cleave110(iHE, 3)
        end do
    else if (nssmat(matp(ink)) == 12) then ! fcc
        do iHE = 1, 6
            jHE = nsp110fcc(iHE)
            cleave110(iHE, 1) = abc(177 + jHE, ink, 1)
            cleave110(iHE, 2) = abc(201 + jHE, ink, 1)
            cleave110(iHE, 3) = abc(225 + jHE, ink, 1)
!            write(969, *) '12', cleave110(iHE, 1), cleave110(iHE, 2), cleave110(iHE, 3)
        end do
    end if

    do iHE = 1, 6
        str110 = sig(1)*cleave110(iHE, 2)**2.0 +sig(2)*cleave110(iHE, 3)**2.0 +sig(4)*cleave110(iHE, 2)*cleave110(iHE, 3)*2.0
        if (abs(str110) > strmax110) then
!            vectn(1) = cleave110(iHE, 2)
!            vectn(2) = cleave110(iHE, 3)
            strmax110 = abs(str110)
        end if
    end do
    
    write(936, 1001) strmax110/yield_stress !       stress110.out
    !     stresscomp.out
!    write(969, *) sig(1), sig(2), sig(3), sig(4), vectn(1), vectn(2)

    1001 format(5x, e20.10)

end
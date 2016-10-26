    subroutine crack_overlap
!!!!!!
!!!!!!    parameter (nume=40000)
!!!!!!    parameter (nume2=20000)
!!!!!!    common/bk06/nprnt,mprint,itmpop,numelt,jprint,idump,locstr
!!!!!!    common /stressflag/ ElemFractCode(nume),ElemDecayCount(nume)
!!!!!!    common /pcracktip/ connect(4,nume2), node(2,nume2),penta(nume2),ndflag(2,nume2), numnpt, numeltu, ndc
!!!!!!	  
!!!!!!    integer numelt, numeltu, ndnum(4), connect, penta
!!!!!!    integer crflag, ele, ne, ele_int
!!!!!!    integer ele2, ne2, ele_int2
!!!!!!    integer ElemFractCode
!!!!!!
!!!!!!    crflag=0
!!!!!!    do ele=1,numeltu
!!!!!!        if((ele<=numelt .and. ElemFractCode(ele)==2) .or.(ele>numelt .and. ElemFractCode(penta(ele-numelt))==2)) then
!!!!!!            do ne=1,4
!!!!!!                if(ne<=3) then
!!!!!!                    ndnum(1)=connect(ne,ele)
!!!!!!                    ndnum(2)=connect(ne+1,ele)
!!!!!!                else if(ne==4) then
!!!!!!                    ndnum(1)=connect(4,ele)
!!!!!!                    ndnum(2)=connect(1,ele)
!!!!!!                end if
!!!!!!				  
!!!!!!                do ele2=1, numeltu
!!!!!!                    if(ele2/=ele) then
!!!!!!                        if((ele2<=numelt .and. ElemFractCode(ele2)==2).or. &
!!!!!!                           (ele2>numelt  .and. ElemFractCode(penta(ele2-numelt))==2)) then	 
!!!!!!                            do ne2=1,4
!!!!!!                                if(ne2<=3) then
!!!!!!                                    ndnum(3)=connect(ne2,ele2)
!!!!!!                                    ndnum(4)=connect(ne2+1,ele2)
!!!!!!                                else if(ne2==4) then
!!!!!!                                    ndnum(3)=connect(4,ele2)
!!!!!!                                    ndnum(4)=connect(1,ele2)
!!!!!!                                end if
!!!!!!
!!!!!!                                call edge_cross(ndnum, crflag)
!!!!!!
!!!!!!                                if(crflag==1) then
!!!!!!                                    if (ele>numelt) then
!!!!!!                                        ele_int=penta(ele-numelt)
!!!!!!                                    else
!!!!!!                                        ele_int=ele
!!!!!!                                    end if
!!!!!!
!!!!!!                                    if (ele2>numelt) then
!!!!!!                                        ele_int2=penta(ele2-numelt)
!!!!!!                                    else
!!!!!!                                        ele_int2=ele2
!!!!!!                                    end if
!!!!!!
!!!!!!                                    write(*,*) 'contact element:',ele_int, ele_int2
!!!!!!                                    stop
!!!!!!
!!!!!!                                end if
!!!!!!								
!!!!!!                            end do
!!!!!!                        end if
!!!!!!                    end if
!!!!!!                end do				  
!!!!!!            end do			  
!!!!!!        end if
!!!!!!    end do
	  
end
								  
									 			  
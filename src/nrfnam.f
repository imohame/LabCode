      subroutine nrfnam(frn,i,nfn)                                      vax750
c     implicit double precision (a-h,o-z)                                    dp
c
c     form the file name for member i+1 of random family kfn
c
c     input arguments
c           frn       family root name (name of the first family member)
c            i        family member index for member i+1
c
c     output arguments
c           nfn       file name for member i+1
c
c     this version written by robert whirley to handle naming probs.
c
      character*8 frn,nfn                                               vax750
      character*42 msg                                                  vax750
      character*1 ni(10)                                                vax750
      data ni/'0','1','2','3','4','5','6','7','8','9'/                  vax750
      data msg/'family member index exceeds 99 for file - '/            vax750
c
      if(i.ne.0)goto 11                                                 vax750
      nfn=frn                                                           vax750
      return                                                            vax750
   11 if(i.lt.100)goto 21                                               vax750
c      call abort(msg//frn)                                              vax750
   21 do 30 k=1,6                                                       vax750
      if(frn(k:k).eq.' ')goto 40                                        vax750
   30 continue                                                          vax750
   40 k=k-1                                                             vax750
      j=i/10                                                            vax750
      nfn=frn(1:k)//ni(j+1)//ni(i-10*j+1)                               vax750
      return                                                            vax750
      end                                                               vax750

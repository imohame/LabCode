      subroutine asgrfm(da,fit)                                         
c
c     assign next family member for random i/o
c
      implicit integer(a-z)                                             
      dimension fit(8)                                                  
c
      common/frfcm1/mxfrf,ifrf,buflen,fcsize,dskloc ,curlen,kop,ier     
      character*8 frfn,frn,kfn                                          
      common/frfcm2/frfn(2,16),frn,kfn                                
      common/double/iprec,ncpw,unit                                    
c
      logical fxist                                                   
      character nfn*8,msg*48                                           
c.... parameter giving number of record units per integer word
c     for most systems a single character is used as a record unit
      data msg/'read requested from nonexistent family member - '/    
c

c.... compute family member index & bias disk address for correct access
      i=da/(fcsize)                               
      da=da-i*fcsize                                                  

c.... get the name of the requested family member
      call nrfnam(frn,i,nfn)                                          
c.... return if current family member is the desired one
      if(kfn.eq.nfn) return                                            
c.... flush the buffer if data is present which is not on disk
      if (fit(3).lt.0) then                                                   
      fit(3)=-fit(3)                                                          
      call wdiska (fit(1),fit(8),buflen,dskloc )                              
      fit(7)=max(fit(7),dskloc +buflen)                                       
      endif                                                                   
c.... determine if requested family member exists
      if (kop.ne.0) go to 20                                                  
      inquire (file=nfn,exist=fxist)                                          
      if (fxist) go to 20                                                     
      if (ier.ne.0) return                                                    
c      call abort (msg)                                                   
c.... close the current family member
   20 continue                                                                
      if (kfn.ne.'        ') then                                             
      close (fit(1),status='keep')                                            
      endif                                                                   
c.... open/create the requested family member
c     lrecl=buflen                                                   
      lrecl=ncpw*buflen                                               
      open (fit(1),file=nfn,access='direct',form='unformatted'                
     1,recl=lrecl,status='unknown')                                           
      kfn=nfn                                                                 
      frfn(2,ifrf)=kfn                                                        
      fit(5)=fit(4)                                                           
      dskloc =fit(5)                                                          
      fit(6)=0                                                                
      curlen=fit(6)                                                           
      fit(7)=0                                                                
      if (fxist) fit(7)=fit(4)                                                
      ier=0                                                                   
      return                                                                  
      end                                                                     

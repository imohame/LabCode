      subroutine rdarg                                                  unix
c
c.... subroutine to parse input command line for unix
c
      character*80 lnkarg,arg                                           unix
      dimension arg(14)                                                 unix
      integer numargs 
      integer ilen, ierror
      common/args/lnkarg,numargs                                        unix
      logical comflg                                                    unix
      data comflg /.false./ 

                                           
      do 10 i=1,80                                                    
         lnkarg(i:i)=' '                                                  
 10   continue                 
                                       
      numargs=iargc()                  
                               
      if (numargs.eq.0) return  
                                       
      if (numargs.gt.8) then                                            
         write(*,1001)                                                     
         stop 'rdarg'                                                      
      endif  
                                                          
      do 20 i=1,numargs                                                 
         call getarg(i,arg(i)) 
* THis is for flyer
*         call pxfgetarg(i,arg(i),ilen,ierror)
 20   continue                                                         

      if (index(arg(1),',') .gt. 0) comflg=.true.         
      if (comflg) then                                   
        do 21 i=1,80                 
           if (arg(1)(i:i).eq.',') then   
              lnkarg(i:i)=' '                
           else                                
              lnkarg(i:i)=arg(1)(i:i)  
           endif                                
 21     continue                     
        return                                  
      endif
                                
      do 40 i=1,numargs                             
         ic=11*(i-1)+1                                      
         il=ic+10                                
         lnkarg(ic:il)=arg(i)                    
 40   continue
                         
      return                                  
 1001 format(//5x,'Too many arguments in rdarg, max is 8')       
      end                                              

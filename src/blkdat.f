      block data blkdat
c     implicit double precision (a-h,o-z)                                    dp
      common/irdmp1/lendr,lenhr,irt,trt,ityprs
      common/bk00/ioff(96)
      common/bk02/ioofc,iphase,imass,lpar(9)
      common/bk05/ifil,iadd,maxsiz,head(12)
      common/bk06/nprnt,mprint,itmpop,numelt,jprint,idump,locstr
      common/bk08/kprint,nstep,ite,ilimit,newstf
      common/bk09/maxref,rhsn,rhsvn,cvtl,iteref,ectl,tolls
      common/bk11/cnwmk(2),iequit,iprint,isref
      common/bk13/xnorm0(6),xnormc(6) !,xnk2d(20)
      common/bk14/lfna(15),lfnt(6)
      common/bk15/cpuio(36),cpuip(36)
      common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
      common/bk18/nummat,ityp2d,ako(31)
      common/bk20/ntotal
      common/bk23/itemp,itherm,irtin
      character*4 mess
      common/bk25/mess
      common/bk32/nsref,nequit,time,timep,lprint,nprint
      common/bk33/irfreq,krfreq,iress
      common/bk49/bulkmx,ncon(30)
      logical rezone
      common/rezone/rezone,nrzn,nctr,irzcnt,nrezon
      common/slar3/nsl,nsntl,nmntl,nslnmx,sltol,slhrd
      common/automt/dtmin,dtmax,mxback,termtm
      common/fissn0/maxneq,mwspac,ntpe0,ntpe1,nfissl(3)
      common/fissn1/melemt,nnns,ntpe2,n2g,llls
      common/fissn3/ifissl,kfissl(3)
      common/total/itrlas,irflas,irhlas,itrtot,irftot,irhtot
c     common/effort/number                                              cray1
      common/bks17/fmult(5)
      common/fmeml/fl
      character*8 names                                                 
      common/filen/names(25)                                            
      common/array/maxa,maxadd,ifield
      common/excute/execut(256)
      character*1 label
      common /alphab/ label(26)
      common/frfcm1/mxfrf,ifrf,buflen,fcsize,dskloc ,curlen,kop,ier     
      character*8 frfn,frn,kfn                                          
      common/frfcm2/frfn(2,16),frn,kfn                                  
c
c
c     version number and compile date
c
      character*8 vn,cdate
      character*72 tx1,tx2
      common/vrsn/vn,cdate,tx1(10),tx2(10)

      data vn/'3.1.3'/,cdate/'09/25/91'/
      data n2g/0/,llls/0/
c     data number/3r680/                                                cray1
      data label /'a','b','c','d','e','f','g','h','i','j',
     1'k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z'/
      data cpuio/36*0.0/
      data cpuip/36*0.0/
c
      data names / 'rstxyz', 'lkjunk', 'disk6' , 'thrxyz','      ',     
     1             '      ', 'disk7' , 'disk3' , 'n2plot','n2dump',     
     2             'n2$$rr', 'rz$$rz', 'strxyz', 'cctxyz','      ',     
     3             '      ', 'result', 'xyz123', 'island','newfle',     
     4             'cct123','       ','       ','       ','      '/     

c     data maxa,maxadd/500000,500000/                                   vms
*      data maxa,maxadd/1500000,1500000/                                 wkstn
*      data maxa,maxadd/5000000,5000000/                                 unics
      data maxa,maxadd/90.0E+06, 90.0E+06/                                     
c     data maxa,maxadd/4000000,4000000/                                 nersc
      data ifield/1/
      data mxfrf/16/                                                    
      data frfn/32*'        '/                                          
c      data ntotal,maxsiz,nwebuf,ite,kprint/4001,1600000,0,0,1/          cray1
      data ntotal,maxsiz,nwebuf,ite,kprint/4001,  262144  ,0,0,1/       
      data mprint,lprint,nprint,nsref,nequit,nstep/1,0,0,0,0,-1/
      data time,timep/0.0,0.0/
      data xnorm0,xnormc/12*0./
      data fmult/5*1./
c      data ntpe0,ntpe1,ntpe2/36,10,0/                                        pl
c     data lfna/9,10,11,12,13,14,15,16,17,18,19,70,70,70,70/                 pl
c     data lfnt/35,36,37,71,71,71/                                           pl
      data ntpe0,ntpe1,ntpe2/17,2,0/                                         
      data lfna/1,2,3,4,7,8,9,10,11,12,13,70,70,70,70/                       
      data lfnt/14,17,18,20,19,21/                                           
      data ifil/0/
      data mess/'  '/
      data dn1,dn2/0.,0./
      data iteref,iress/0,10000000/
      data ncon/4,4,8,8,8,8,8,9,4,5,10,8,6,9,9,16,8,10,6,8,
     1          8,4,9,11,7,5*0/
      data iadd,itherm/0,0/
      data fl/0.90/
      data ioofc/1/
      data ioff/96*0/
      data rezone/.false./
      data ako/1.,1.,0.,1.,1.,1.,1.,0.,1.,0.,1.,1.,1.,1.,0.,1.,1.,1.,1.,
     1         1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1./
      data nctr/10000000/
      data sltol,slhrd/.0010,.0100/
      data irzcnt/100000/
      data nrezon/0/
      data nfissl/3*20/
      data kfissl/3*0/
      data itrtot,irftot,irhtot/0,0,0/
      data dtmin,dtmax/0.0,0.0/								 
      end

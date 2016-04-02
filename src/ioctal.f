      integer function ioctal (l)
      l1=l-(l/10)*10
      l2=(l-(l/100)*100-l1)/10
      l3=(l-(l/1000)*1000-l1-l2*10)/100
      l4=(l-(l/10000)*10000-l1-l2*10-l3*100)/1000
      ioctal=l1+8*l2+64*l3+512*l4
      return
      end

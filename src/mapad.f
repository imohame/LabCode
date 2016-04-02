      function mapad(idiag,n,mxc)
c
c     compute coefficient matrix profile and set diagonal pointers
c
      dimension idiag(*)
c
      mxc=0
      lk=mxc
      do 100 i=1,n
      mxc=max(mxc,idiag(i))
      lk=lk+idiag(i)+1
      idiag(i)=lk
  100 continue
      mxc=mxc+1
      mapad=lk
      return
      end

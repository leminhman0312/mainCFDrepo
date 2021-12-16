c program to read a grid in Charley Swanson's format and
c output it for ISAAC (plot3d binary)
c
      parameter (idim = 225, jdim = 33)
      implicit real*8 (a-h,o-z)
      dimension x(idim,jdim), y(idim,jdim)
c
c read in Charley's grid
c
      open(8,file='gr224n20',status='old')
      read(8,*)
      read(8,*)
      read(8,*)
      read(8,*)xi,xj,cord
      read(8,*)ni,nj,nte1,nte2,nbot
      if(ni.ne.idim.or.nj.ne.jdim)then
        write(*,'(a,24i5)')'Input grid size differs from specified ',
     1     idim, jdim, ni, nj
         stop
      endif
      read(8,1000)((x(i,j),i=1,ni),j=1,nj)
      read(8,1000)((y(i,j),i=1,ni),j=1,nj)
 1000 format(10E13.7)
c
c output grid
c
      inc=1
      niout=idim
      njout=jdim
      open(9,file='n12_225_33.fmt',status='unknown',
     1      form='formatted')
      write(9,*)niout,njout
      write(9,*)((x(i,j),i=1,idim,inc),j=1,jdim,inc),
     1          ((y(i,j),i=1,idim,inc),j=1,jdim,inc)
      close(9)
      stop
      end

c Program to generate a grid for Klebanoff's flatplate
      parameter(idim=65,jdim=65)
      implicit real*8 (a-h,o-z)
      dimension x(idim,jdim),
     1          y(idim,jdim),
     2          s1(jdim)
c
c ytop  +----------------------------+
c       |                            |
c       |                            |
c       |                            |
c       |                            |
c       |                            |
c ybot  +----------------------------+
c       |                            |
c       xinit                        xfinal
c
c Flatplate data:
c
      xinit=-0.5D0
      xfinal=1.5D0
      ybot=0.0D0
      ytop=0.18D0
c
c Set up equal spaced x distribution
c
      ds=(xfinal-xinit)/float(idim-1)
c
c Y distribution:
c
      strtch=1.00025
      call stretch(jdim,strtch,s1)
c
c Calculate i=1 distribution
c
      do j=1,jdim
        x(1,j)=xinit
        y(1,j)=ybot+s1(j)*(ytop-ybot)
      enddo
c
c Copy i=1 distribution throughout domain
c
      do i=2,idim
        xsta=xinit+float(i-1) * ds
        do j=1,jdim
          x(i,j)=xsta
          y(i,j)=y(1,j)
        enddo
      enddo
c
c Output grid
c
      open(7,file='flatplate.grd',status='unknown')
      call wrp3d(idim,jdim,x,y)
      close(7)
      stop
      end
      subroutine stretch(jdim,strtch,s)
c
c Subroutine to set up a non-dimensionalized stretching function
c
      implicit real*8(a-h,o-z)
      dimension s(jdim)
      rh(b,psi)=1D0-b+2D0*b/(1D0+((b+1D0)/(b-1D0))**psi)
      s(1)=0D0
      do j=2,jdim
        ps=float(jdim-j)/float(jdim-1)
        s(j)=rh(strtch,ps)
      enddo
      return
      end
      subroutine wrp3d(idim,jdim,x,y)
c
c write a binary plot3d grid file
c
      implicit real*8(a-h,o-z)
      dimension x(idim,jdim),
     1          y(idim,jdim)
      write(7,*)idim,jdim
      write(7,*)((x(i,j),i=1,idim),j=1,jdim),
     1         ((y(i,j),i=1,idim),j=1,jdim)
      return
      end

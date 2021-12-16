      subroutine aires
c     -----------------------------------------------------------------
c     computation of the triangle and control volume areas
c     -----------------------------------------------------------------
      include 'nsc2ke.inc'
c     -----------------------------------------------------------------
      real    a(2), b(2), c(2)
c
      do 5 is=1,ns
         airs(is)  = 0.0
5     continue
c
      do 10 jt=1,nt
c
         is1       = nu(1,jt)
         is2       = nu(2,jt)
         is3       = nu(3,jt)
c
c        computing the current triangle area
c
         a(1)      = coor(1,is1)
         a(2)      = coor(2,is1)
         b(1)      = coor(1,is2)
         b(2)      = coor(2,is2)
         c(1)      = coor(1,is3)
         c(2)      = coor(2,is3)
         u         = a(2)-c(2)
         v         = c(1)-a(1)
         w         = a(1)*c(2)-a(2)*c(1)
         vol       = u*b(1)+v*b(2)+w
         airt(jt)  = 0.5*abs(vol)
c
c        gathering the current triangle area into the
c        control volume one
c
         airs(is1) = airs(is1)+airt(jt)/3.0
         airs(is2) = airs(is2)+airt(jt)/3.0
         airs(is3) = airs(is3)+airt(jt)/3.0
c
10    continue
c
        som_airt=0.
        do it=1,nt
         som_airt=som_airt+airt(it)
        enddo
c
        som_airs=0.
        do is=1,ns
         som_airs=som_airs+airs(is)
        enddo
        print *,'test airt , som airt =',som_airt
        print *,'test airs , som airs =',som_airs
c
      return
      end
      subroutine axigeo
      include 'nsc2ke.inc'
      real vnt(nn)
      if(iaxi.eq.0) then
       do is=1,ns
        airsa(is)=airs(is)
        rrs(is)=1.0
       enddo
       do it=1,nt
        airta(it)=1.0/airt(it)
       enddo
       return
      endif
c
       do it=1,nt
         is1=nu(1,it)
         is2=nu(2,it)
         is3=nu(3,it)
         y1=coor(2,is1)
         y2=coor(2,is2)
         y3=coor(2,is3)
         rrt=(y1+y2+y3)/3.0
         rrs(is1)=rrs(is1)+airt(it)*rrt
         rrs(is2)=rrs(is2)+airt(it)*rrt
         rrs(is3)=rrs(is3)+airt(it)*rrt
         airta(it)=rrt/airt(it)
       enddo
c
       do  is=1,ns
         if(coor(2,is).lt.1.e-6) then
           rrs(is)=rrs(is)/(3.*airs(is))
                                else
           rrs(is)=coor(2,is)
         endif
        airsa(is)=airs(is)*rrs(is)
        vnt(is)=0.
       enddo
c
c redefinition de vno pour l'axi
c
      do ifr=1,nfr
         iseg          = nufr(ifr)
         nuor          = nubo(1,iseg)
         nuex          = nubo(2,iseg)
         r             =(coor(2,nuor)+coor(2,nuex))*0.5
         vnt(nuex)=vnt(nuex)+0.5*r
         vnt(nuor)=vnt(nuor)+0.5*r
      enddo
c
      do is=1,ns
         vno(is)   = vno(is)*vnt(is)
      enddo
c
       do iseg=1,nseg
          is1=nubo(1,iseg)
          is2=nubo(2,iseg)
          r=(coor(2,is1)+coor(2,is2))*0.5
          vnocl(3,iseg)=vnocl(3,iseg)*r
       enddo
c
         return
         end
      subroutine caldtl(dtmin,dt)
      include 'nsc2ke.inc'
c
      cste              = 110.*tinf/tinfd
cc
cc       sutherland   viscosity
cc
         if (ivis .ne. 0) then
         do is=1,ns
         temp = pres(is)/(gam1*ua(1,is))
         xnum=tinf+cste
         xden=temp+cste
         reylam(is)=reynolds*(temp/tinf)**1.5*(xnum/xden)
         enddo
         endif
c
            if(iloc.eq.2) then
            xcoef=2.*gam/pr
            else
            xcoef=0.0
            endif
            dtmin              = 1.0e+12
            do is=1,ns
               reytot=reylam(is)+reyturb(is)
               dsrey= xcoef*reytot/ua(1,is)
               unorm=sqrt(ua(2,is)*ua(2,is) +
     &                    ua(3,is)*ua(3,is))
               sigmax=gam*pres(is)/ua(1,is)
               sigmax=max(1.e-2,sqrt(sigmax)+unorm)
               dtl(is)       = cfl*dthaut(is)*dthaut(is)/
     &                        (dthaut(is)*sigmax + ivis*dsrey)
               dtmin           = min(dtmin, dtl(is))
            enddo
c
            if (iloc .eq. 0) then
               do is=1,ns
                  dtl(is)    = dtmin
               enddo
            endif
c
            dt               = min(dtmin, tmax-t)
            t                = t + dt
c
           return
           end
      SUBROUTINE calprc
      include 'nsc2ke.inc'
c
       XRHO = 1.E-6
       DO 100 is=1,ns
       if(ua(1,is).lt.xrho) then
              print *,logfr(is),coor(1,is),coor(2,is),pres(is)
              stop 'density < 0 , calprc '
       endif
       pres(is)=gam1*(ua(4,is) - 0.5*(ua(2,is)*ua(2,is)+
     &               ua(3,is)*ua(3,is))/ua(1,is))
        IF(pres(is).lt.xrho) then
              print *,logfr(is),coor(1,is),coor(2,is),pres(is)
              stop 'pressure < 0 calprc '
       endif
100    CONTINUE
C
      RETURN
      END
      subroutine cdl
c     -----------------------------------------------------------------
c     computation of the convective fluxes at  boundaries
c     -----------------------------------------------------------------
      include 'nsc2ke.inc'
c     -----------------------------------------------------------------
c     local variables definition
      integer  is  ,  i
      real     cap , capd1, capd2, uq41, uq42,
     &         u, v, c    , usc  , vp1 , vp3 , vp4,
     &         fdc1, fdc2 , fdc3 , fdc4,
     &         fgp1, fgp2 , fgp3 , fgp4
c
c     vertices on the body, slip boundary condition
c
      do 700 i=1,nslog2
         is=log2(i)
         ce(2,is)  = ce(2,is) - vno(is)*vnox(is)*pres(is)
         ce(3,is)  = ce(3,is) - vno(is)*vnoy(is)*pres(is)
700   continue
c
c solid bodies
c
      if(ivis.eq.1.and.ilaw.eq.0) then
c
      do 100 i=1,nslog3
         is=log3(i)
         ce(2,is)  = 0.0
         ce(3,is)  = 0.0
         ce(5,is)  = 0.0
         ce(6,is)  = 0.0
100   continue
c
      if(iecc.eq.2) then
      do 101 i=1,nslog3
         is=log3(i)
         ce(4,is)  = ce(1,is)*tbrd2
101   continue
      endif
      endif
c
c     vertices on the upstream boundary
c
      do 850 i=1,nslog5
         is=log5(i)
         cap       = vno(is)
         capd1     = vnox(is)
         capd2     = vnoy(is)
         u         = ua(2,is)/ua(1,is)
         v         = ua(3,is)/ua(1,is)
         c         = sqrt(gam*pres(is)/ua(1,is))
         usc       = 1.0/c
         uq41      = 0.5*(u*u + v*v)*usc + gam4*c
         uq42      = capd1*u + capd2*v
c
c     computation of a-(wi).winf
c
         vp1       = (capd1*u + capd2*v)*cap
         vp3       = min(vp1 + c*cap, 0.0)
         vp4       = min(vp1 - c*cap, 0.0)
         vp1       = min(vp1        , 0.0)
         fdc1  = vp1*(roin + gam1*usc*usc*(-0.5*(u*u + v*v)*roin   +
     &                                      u*ruxin + v*ruyin - ein))
         fdc2  = vp1*((capd1*v - capd2*u)*roin +
     &                    capd2*ruxin - capd1*ruyin)
         fdc3  = vp3*(0.5*(-uq42*roin + capd1*ruxin + capd2*ruyin) +
     &                    0.5*gam1*usc*(0.5*(u*u + v*v)*roin -
     &                                  u*ruxin - v*ruyin + ein))
         fdc4  = vp4*(0.5*(uq42*roin - capd1*ruxin - capd2*ruyin)  +
     &                    0.5*gam1*usc*(0.5*(u*u + v*v)*roin -
     &                                  u*ruxin - v*ruyin + ein))
         fgp1  = fdc1 + usc*(fdc3 + fdc4)
         fgp2  = u*fdc1 + capd2*fdc2  +
     &               (u*usc + capd1)*fdc3 + (u*usc - capd1)*fdc4
         fgp3  = v*fdc1 - capd1*fdc2  +
     &               (v*usc + capd2)*fdc3 + (v*usc - capd2)*fdc4
         fgp4  = 0.5*(u*u + v*v)*fdc1 + (capd2*u - capd1*v)*fdc2   +
     &               (uq41 + uq42)*fdc3   + (uq41 - uq42)*fdc4
c
         ce(1,is)  = ce(1,is) - fgp1
         ce(2,is)  = ce(2,is) - fgp2
         ce(3,is)  = ce(3,is) - fgp3
         ce(4,is)  = ce(4,is) - fgp4
         ce(5,is)  = ce(5,is) - fgp1*ua(5,is)/ua(1,is)
         ce(6,is)  = ce(6,is) - fgp1*ua(6,is)/ua(1,is)
c
c     computation of a+(wi).wi
c
         vp1   = (capd1*u + capd2*v)*cap
         vp3   = max(vp1 + c*cap, 0.0)
         vp4   = max(vp1 - c*cap, 0.0)
         vp1   = max(vp1        , 0.0)
         fdc1  = vp1*(ua(1,is) +
     &               gam1*usc*usc*(-0.5*(u*u + v*v)*ua(1,is)    +
     &                    u*ua(2,is) + v*ua(3,is) - ua(4,is)))
         fdc2  = vp1*((capd1*v - capd2*u)*ua(1,is) +
     &                    capd2*ua(2,is) - capd1*ua(3,is))
         fdc3  = vp3*(0.5*(-uq42*ua(1,is) +
     &                capd1*ua(2,is) + capd2*ua(3,is)) +
     &             0.5*gam1*usc*(0.5*(u*u + v*v)*ua(1,is)     -
     &                           u*ua(2,is) - v*ua(3,is)      +
     &                                  ua(4,is)))
         fdc4  = vp4*(0.5*(uq42*ua(1,is)  -
     &                   capd1*ua(2,is) - capd2*ua(3,is))      +
     &              0.5*gam1*usc*(0.5*(u*u + v*v)*ua(1,is)     -
     &                            u*ua(2,is) - v*ua(3,is)      +
     &                                ua(4,is)))
         fgp1  = fdc1 + usc*(fdc3 + fdc4)
         fgp2  = u*fdc1 + capd2*fdc2  +
     &               (u*usc + capd1)*fdc3 + (u*usc - capd1)*fdc4
         fgp3  = v*fdc1 - capd1*fdc2  +
     &               (v*usc + capd2)*fdc3 + (v*usc - capd2)*fdc4
         fgp4  = 0.5*(u*u + v*v)*fdc1 + (capd2*u - capd1*v)*fdc2 +
     &               (uq41 + uq42)*fdc3   + (uq41 - uq42)*fdc4
c
         ce(1,is)  = ce(1,is) - fgp1
         ce(2,is)  = ce(2,is) - fgp2
         ce(3,is)  = ce(3,is) - fgp3
         ce(4,is)  = ce(4,is) - fgp4
         ce(5,is)  = ce(5,is) - fgp1*ua(5,is)/ua(1,is)
         ce(6,is)  = ce(6,is) - fgp1*ua(6,is)/ua(1,is)
850   continue
c
c     vertices on the downstream boundary
c
      do 950 i=1,nslog4
         is=log4(i)
         cap       = vno(is)
         capd1     = vnox(is)
         capd2     = vnoy(is)
         u         = ua(2,is)/ua(1,is)
         v         = ua(3,is)/ua(1,is)
         c         = sqrt(gam*pres(is)/ua(1,is))
         usc       = 1.0/c
         uq41      = 0.5*(u*u + v*v)*usc + gam4*c
         uq42      = capd1*u + capd2*v
c
c  constant invariants of riemann for u-c
c
         pstar=pout
         xroout=ua(1,is)*(pstar/pres(is))**(1./gam)
         unorm=sqrt(u*u+v*v)
         unorm=unorm+2./gam1*(c-sqrt(gam*pstar/xroout))
         xruxout=ruxout*unorm*xroout
         xruyout=ruyout*unorm*xroout
         xeout=pstar/gam1+xroout*unorm*unorm*0.5
c
c     computation of a-(wi).winf
c
         vp1       = (capd1*u + capd2*v)*cap
         vp3       = min(vp1 + c*cap, 0.0)
         vp4       = min(vp1 - c*cap, 0.0)
         vp1       = min(vp1        , 0.0)
         fdc1      = vp1*(xroout +
     &                    gam1*usc*usc*(-0.5*(u*u + v*v)*xroout +
     &                                  u*xruxout + v*xruyout - xeout))
         fdc2      = vp1*((capd1*v - capd2*u)*xroout +
     &                    capd2*xruxout - capd1*xruyout)
         fdc3      = vp3*(0.5*(-uq42*xroout  +
     &                         capd1*xruxout + capd2*xruyout) +
     &                    0.5*gam1*usc*(0.5*(u*u + v*v)*xroout  -
     &                                  u*xruxout - v*xruyout + xeout))
         fdc4      = vp4*(0.5*(uq42*xroout   -
     &                         capd1*xruxout - capd2*xruyout) +
     &                    0.5*gam1*usc*(0.5*(u*u + v*v)*xroout  -
     &                                  u*xruxout - v*xruyout + xeout))
         fgp1      = fdc1 + usc*(fdc3 + fdc4)
         fgp2      = u*fdc1 + capd2*fdc2  +
     &               (u*usc + capd1)*fdc3 + (u*usc - capd1)*fdc4
         fgp3      = v*fdc1 - capd1*fdc2  +
     &               (v*usc + capd2)*fdc3 + (v*usc - capd2)*fdc4
         fgp4    = 0.5*(u*u + v*v)*fdc1 + (capd2*u - capd1*v)*fdc2 +
     &               (uq41 + uq42)*fdc3   + (uq41 - uq42)*fdc4
c
         ce(1,is)  = ce(1,is) - fgp1
         ce(2,is)  = ce(2,is) - fgp2
         ce(3,is)  = ce(3,is) - fgp3
         ce(4,is)  = ce(4,is) - fgp4
         ce(5,is)  = ce(5,is) - fgp1*ua(5,is)/ua(1,is)
         ce(6,is)  = ce(6,is) - fgp1*ua(6,is)/ua(1,is)
c
c     computation of a+(wi).wi
c
         vp1       = (capd1*u + capd2*v)*cap
         vp3       = max(vp1 + c*cap, 0.0)
         vp4       = max(vp1 - c*cap, 0.0)
         vp1       = max(vp1        , 0.0)
         fdc1      = vp1*(ua(1,is) +
     &                    gam1*usc*usc*(-0.5*(u*u + v*v)*ua(1,is)    +
     &                    u*ua(2,is) + v*ua(3,is) - ua(4,is)))
         fdc2      = vp1*((capd1*v - capd2*u)*ua(1,is) +
     &                    capd2*ua(2,is) - capd1*ua(3,is))
         fdc3      = vp3*(0.5*(-uq42*ua(1,is) +
     &                    capd1*ua(2,is) + capd2*ua(3,is)) +
     &                    0.5*gam1*usc*(0.5*(u*u + v*v)*ua(1,is)     -
     &                                  u*ua(2,is) - v*ua(3,is)      +
     &                                  ua(4,is)))
         fdc4   = vp4*(0.5*(uq42*ua(1,is)  -
     &                         capd1*ua(2,is) - capd2*ua(3,is))      +
     &                    0.5*gam1*usc*(0.5*(u*u + v*v)*ua(1,is)     -
     &                                  u*ua(2,is) - v*ua(3,is)      +
     &                                  ua(4,is)))
         fgp1      = fdc1 + usc*(fdc3 + fdc4)
         fgp2      = u*fdc1 + capd2*fdc2  +
     &               (u*usc + capd1)*fdc3 + (u*usc - capd1)*fdc4
         fgp3      = v*fdc1 - capd1*fdc2  +
     &               (v*usc + capd2)*fdc3 + (v*usc - capd2)*fdc4
         fgp4    = 0.5*(u*u + v*v)*fdc1 + (capd2*u - capd1*v)*fdc2 +
     &               (uq41 + uq42)*fdc3   + (uq41 - uq42)*fdc4
c
         ce(1,is)  = ce(1,is) - fgp1
         ce(2,is)  = ce(2,is) - fgp2
         ce(3,is)  = ce(3,is) - fgp3
         ce(4,is)  = ce(4,is) - fgp4
         ce(5,is)  = ce(5,is) - fgp1*ua(5,is)/ua(1,is)
         ce(6,is)  = ce(6,is) - fgp1*ua(6,is)/ua(1,is)
950   continue
c
      do 1000 i=1,nslog6
         is=log6(i)
         ce(1,is)  = 0.0
         ce(2,is)  = 0.0
         ce(3,is)  = 0.0
         ce(4,is)  = 0.0
         ce(5,is)  = 0.0
         ce(6,is)  = 0.0
1000   continue
      return
      end
      subroutine clhaut
c     -----------------------------------------------------------------
c     computation of the control volume altitudes
c     -----------------------------------------------------------------
      include 'nsc2ke.inc'
c     -----------------------------------------------------------------
c     local variables definition
      integer nubo1  , nubo2  , nubo3  , is , jt
      real    dbxx(3), dbyy(3), haut(3), htt, ait,
     &        x1, x2 , x3, y1 , y2, y3
c
c     initializing the vertex altitudes
c
      do 77 is=1,ns
         dthaut(is)     = 1.0e+12
77    continue
c
      do 1000 jt=1,nt
c
         nubo1          = nu(1,jt)
         nubo2          = nu(2,jt)
         nubo3          = nu(3,jt)
         x1             = coor(1,nubo1)
         y1             = coor(2,nubo1)
         x2             = coor(1,nubo2)
         y2             = coor(2,nubo2)
         x3             = coor(1,nubo3)
         y3             = coor(2,nubo3)
         dbxx(1)        = y2 - y3
         dbxx(2)        = y3 - y1
         dbxx(3)        = y1 - y2
         dbyy(1)        = x3 - x2
         dbyy(2)        = x1 - x3
         dbyy(3)        = x2 - x1
c
         ait            = 1.0/(2.0*airt(jt))
         haut(1)        = ait*(sqrt(dbxx(1)*dbxx(1) + dbyy(1)*dbyy(1)))
         haut(2)        = ait*(sqrt(dbxx(2)*dbxx(2) + dbyy(2)*dbyy(2)))
         haut(3)        = ait*(sqrt(dbxx(3)*dbxx(3) + dbyy(3)*dbyy(3)))
c
         htt            = min(1.0/haut(1), 1.0/haut(2))
         htt            = min(htt        , 1.0/haut(3))
         dthaut(nubo1)  = min(dthaut(nubo1), htt)
         dthaut(nubo2)  = min(dthaut(nubo2), htt)
         dthaut(nubo3)  = min(dthaut(nubo3), htt)
c
1000  continue
c
      return
      end
      subroutine config
c     -----------------------------------------------------------------
c     reading the application dependent parameters
c     -----------------------------------------------------------------
      include 'nsc2ke.inc'
c     -----------------------------------------------------------------
c
c     perferct gaz constant
c
      gam     = 1.4
      gam1    = gam - 1.0
      gam4    = 1.0/gam1
c
c     prandtl's number
c
      pr      = 0.72
      prt     = 0.9
c
cccccccccccccccccccccccccccccccccccccccccccc
      open(1,file="data.inp",status='old')
      read(1, *) iaxi
      read(1, *) ivis
      read(1, *) reynolds
      read(1, *) froud
      read(1, *) xmach
      read(1, *) xpcoef
      read(1, *) iecc
      read(1, *) tinfd
      read(1, *) tbrd1
      read(1, *) tetadeg
      read(1, *) iflux
      read(1, *) nordre
      read(1, *) iloc
      read(1, *) cfl
      read(1, *) ktmax
      read(1, *) ifre
      read(1, *) tmax
      read(1, *) resf
      read(1, *) ncont
      read(1, *)
      read(1, *) iturb
      read(1, *) ilaw
      read(1, *) delta
      read(1, *) ncont1
      read(1, *) xtmin,xtmax,ytmin,ytmax
      close(1)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      print *,'*********************************************'
      print *,'        NSC2KE : RELEASE 1.0, 1994           '
      print *,'*********************************************'
      if(iaxi.eq.1) print *,'AXI SYMMETRIC        '
      if(iaxi.eq.0) print *,'2D              '
      if(ivis.eq.1) print *,'Navier-Stokes computation '
      if(ivis.eq.0) print *,'Euler computation      '
      print *,'*********************************************'
      if(iturb.eq.1) then
          print *,'k-epsilon turbulence model      '
      if(ilaw.eq.1) print *,'classical wall-law  technique  '
      if(ilaw.eq.0) print *,'two-layer technique            '
      endif
      print *,'*********************************************'
c
c     mach number, angle of attack
c
      if(iecc.eq.1) print *,'adiabatic walls '
      if(iecc.eq.2) print *,'isothermal walls '
      print*,'free stream mach number  : ', xmach
      print*,'angle of attack          : ', tetadeg
c
      if (ivis .eq. 1) then
      print*,'reynolds number          : ', reynolds
      reynolds = 1./reynolds
      endif
      print*,'Froude number            : ', froud
c
        pin     = 1.0/(gam*xmach**2)
        tinf    = 1.0/(gam*gam1*xmach*xmach)
c
      tbrd2=tinf*tbrd1/tinfd
c
c     free stream uniform solution
c     external flow around a body
c
      roin    = 1.0
      teta    = tetadeg*3.141592654/180.0
      uxin    = cos(teta)
      uyin    = sin(teta)
      ruxin   = roin*uxin
      ruyin   = roin*uyin
      ein     = 0.5*roin*(uxin*uxin+uyin*uyin)+pin/gam1
      roout   = roin
      uxout   = uxin
      uyout   = uyin
      ruxout  = ruxin
      ruyout  = ruyin
      pout    = pin*xpcoef
      eout    = 0.5*roout*(uxout**2+uyout**2)+pout/gam1
      xkin    = 1.e-5
      xein    = 1.e-5
c
      print*,'courant number           : ', cfl
c
      if(iloc.ne.0) print *,'local time stepping'
c
c     upwinding parameter
c
      beta=0.33333
      beta1   = 1.0 - beta
c
      if (iflux .eq. 1) print*,'Roe''s scheme'
      if (iflux .eq. 2) print*,'Osher''s scheme'
      if (iflux .eq. 3) print*,'Kinetic''s scheme'
c
c     runge-kutta time integration process parameters
c
c      irk=3
c      alpha(1)=0.3333
c      alpha(2)=0.5
c      alpha(3)=1.0
ccccccccccccccccccccccccccccccccc
       irk=4
       alpha(1)=0.11
       alpha(2)=0.2766
       alpha(3)=0.5
       alpha(4)=1.0
cccccccccccccccccccccccccccccccc
      print *,'*********************************************'
      close(1)
      return
      end
      subroutine flucin
c     -----------------------------------------------------------------
c     edgewise convective fluxes computation using the kinetic boltzmann
c     flux.
c     the nodal gradients are computed using a beta-combination
c     of centered and hermitian (half-upwind) gradients
c
c     bijan mohammadi, inria-menusin
c     -----------------------------------------------------------------
      include 'nsc2ke.inc'
c     -----------------------------------------------------------------
c     local variables definition
      integer  nsg    , nubo1  , nubo2
      real dpm(4),dpor(4),dpex(4),aux1(4),aux2(4),
     1     gradi(4),gradj(4)
c
      gam4m5       = (2.-gam)/gam1
      ra3          = sqrt(3.0)
      alp          = 1./2./ra3
c
      beta2 = beta
      beta3 = 0.5*(1.0 - 2.0*beta)
      if(nordre.eq.2) then
      beta2 = 2.0*beta
      beta3 = 1.0 - 2.0*beta
      endif
      e2           = 1.e-6
c
c     loop on global list of edges
c
      do 500 nsg=1,nseg
c
         nubo1     = nubo(1,nsg)
         nubo2     = nubo(2,nsg)
c
         xnn       =-vnocl(1,nsg)
         ynn       =-vnocl(2,nsg)
         rnn       = vnocl(3,nsg)
c
         aix       = coor(1,nubo2)-coor(1,nubo1)
         aiy       = coor(2,nubo2)-coor(2,nubo1)
c
         uas11     = ua(1,nubo1)
         uas12     = ua(1,nubo2)
         uas21     = ua(2,nubo1)
         uas22     = ua(2,nubo2)
         uas31     = ua(3,nubo1)
         uas32     = ua(3,nubo2)
         uas41     = pres(nubo1)
         uas42     = pres(nubo2)
c
         if(nordre.eq.1)  goto 1234
c
         gradI(1)  = beta2*(aix*dx(1,nubo1) + aiy*dy(1,nubo1)) +
     &               beta3*(uas12 - uas11)
c
         gradI(2)  = beta2*(aix*dx(2,nubo1) + aiy*dy(2,nubo1)) +
     &               beta3*(uas22 - uas21)
c
         gradI(3)  = beta2*(aix*dx(3,nubo1) + aiy*dy(3,nubo1)) +
     &               beta3*(uas32 - uas31)
c
         gradI(4)  = beta2*(aix*dx(4,nubo1) + aiy*dy(4,nubo1)) +
     &               beta3*(uas42 - uas41)
c
         gradJ(1)  = beta2*(aix*dx(1,nubo2) + aiy*dy(1,nubo2)) +
     &               beta3*(uas12 - uas11)
c
         gradJ(2)  = beta2*(aix*dx(2,nubo2) + aiy*dy(2,nubo2)) +
     &               beta3*(uas22 - uas21)
c
         gradJ(3)  = beta2*(aix*dx(3,nubo2) + aiy*dy(3,nubo2)) +
     &               beta3*(uas32 - uas31)
c
         gradJ(4)  = beta2*(aix*dx(4,nubo2) + aiy*dy(4,nubo2)) +
     &               beta3*(uas42 - uas41)
c
         if(nordre.eq.2) then
c
         uas11     = uas11 + 0.5*gradI(1)
         uas21     = uas21 + 0.5*gradI(2)
         uas31     = uas31 + 0.5*gradI(3)
         uas41     = uas41 + 0.5*gradI(4)
c
         uas12     = uas12 - 0.5*gradJ(1)
         uas22     = uas22 - 0.5*gradJ(2)
         uas32     = uas32 - 0.5*gradJ(3)
         uas42     = uas42 - 0.5*gradJ(4)
c
         elseif(nordre.eq.3) then
c
         dpm(1)    =-(uas12 - uas11)
         dpex(1)   =-4.0*gradI(1) - dpm(1)
         aux1(1)   = 0.25*(1.0+  SIGN(1.0, dpex(1)*dpm(1)))
         dpor(1)   =-4.0*gradJ(1) - dpm(1)
         aux2(1)   = 0.25*(1.0 + SIGN(1.0, dpor(1)*dpm(1)))
c
         dpm(2)    =-(uas22 - uas21)
         dpex(2)   =-4.0*gradI(2) - dpm(2)
         aux1(2)   = 0.25*(1.0 + SIGN(1.0, dpex(2)*dpm(2)))
         dpor(2)   =-4.0*gradJ(2) - dpm(2)
         aux2(2)   = 0.25*(1.0 + SIGN(1.0, dpor(2)*dpm(2)))
c
         dpm(3)    =-(uas32 - uas31)
         dpex(3)   =-4.0*gradI(3) - dpm(3)
         aux1(3)   = 0.25*(1.0 + SIGN(1.0, dpex(3)*dpm(3)))
         dpor(3)   =-4.0*gradJ(3) - dpm(3)
         aux2(3)   = 0.25*(1.0 + SIGN(1.0, dpor(3)*dpm(3)))
c
         dpm(4)    =-(uas42 - uas41)
         dpex(4)   =-4.0*gradI(4) - dpm(4)
         aux1(4)   = 0.25*(1.0 + SIGN(1.0, dpex(4)*dpm(4)))
         dpor(4)   =-4.0*gradJ(4) - dpm(4)
         aux2(4)   = 0.25*(1.0 + SIGN(1.0, dpor(4)*dpm(4)))
c
         gradI(1)  = aux1(1)*
     &               ((dpex(1)*dpex(1) + e2)*dpm(1)+
     &                (dpm(1)*dpm(1)   + e2)*dpex(1))/
     &               (dpex(1)*dpex(1)  + dpm(1)*dpm(1) + 2.0*e2)
         gradJ(1)  = aux2(1)*
     &               ((dpor(1)*dpor(1) + e2)*dpm(1)+
     &                (dpm(1)*dpm(1)   + e2)*dpor(1))/
     &               (dpor(1)*dpor(1)  + dpm(1)*dpm(1) + 2.0*e2)
c
         gradI(2)  = aux1(2)*
     &               ((dpex(2)*dpex(2) + e2)*dpm(2)+
     &                (dpm(2)*dpm(2)   + e2)*dpex(2))/
     &               (dpex(2)*dpex(2)  + dpm(2)*dpm(2) + 2.0*e2)
         gradJ(2)  = aux2(2)*
     &               ((dpor(2)*dpor(2) + e2)*dpm(2)+
     &                (dpm(2)*dpm(2)   + e2)*dpor(2))/
     &               (dpor(2)*dpor(2)  + dpm(2)*dpm(2) + 2.0*e2)
c
         gradI(3)  = aux1(3)*
     &               ((dpex(3)*dpex(3) + e2)*dpm(3)+
     &                (dpm(3)*dpm(3)   + e2)*dpex(3))/
     &               (dpex(3)*dpex(3)  + dpm(3)*dpm(3) + 2.0*e2)
         gradJ(3)  = aux2(3)*
     &               ((dpor(3)*dpor(3) + e2)*dpm(3)+
     &                (dpm(3)*dpm(3)   + e2)*dpor(3))/
     &               (dpor(3)*dpor(3)  + dpm(3)*dpm(3) + 2.0*e2)
c
         gradI(4)  = aux1(4)*
     &               ((dpex(4)*dpex(4) + e2)*dpm(4)+
     &                (dpm(4)*dpm(4)   + e2)*dpex(4))/
     &               (dpex(4)*dpex(4)  + dpm(4)*dpm(4) + 2.0*e2)
         gradJ(4)  = aux2(4)*
     &               ((dpor(4)*dpor(4) + e2)*dpm(4)+
     &                (dpm(4)*dpm(4)   + e2)*dpor(4))/
     &               (dpor(4)*dpor(4)  + dpm(4)*dpm(4) + 2.0*e2)
c
         uas11     = uas11 - gradI(1)
         uas21     = uas21 - gradI(2)
         uas31     = uas31 - gradI(3)
         uas41     = uas41 - gradI(4)
c
         uas12     = uas12 + gradJ(1)
         uas22     = uas22 + gradJ(2)
         uas32     = uas32 + gradJ(3)
         uas42     = uas42 + gradJ(4)
c
         endif
c
1234   continue
c
c rotation
c
         uas210 = uas21
         uas21  = xnn*uas210+ynn*uas31
         uas31  =-ynn*uas210+xnn*uas31
         uas41  = uas41/uas11
         uas220 = uas22
         uas22  = xnn*uas220+ynn*uas32
         uas32  =-ynn*uas220+xnn*uas32
         uas42  = uas42/uas12
c
         rat1 = sqrt(uas41)
         ext1 = min(ra3,max(-ra3,-uas21/rat1))
c
         a01  = alp*(ra3-ext1)
         a11  = alp*(ra3**2-ext1**2)/2.
         a21  = alp*(ra3**3-ext1**3)/3.
         a31  = alp*(ra3**4-ext1**4)/4.
c
         flu11= uas11*(uas21*a01+rat1*a11)
         flu21= uas11*(uas21*uas21*a01+2.*uas21*rat1*a11
     1                 +uas41*a21)
         flu31= uas31*flu11
         flu41= uas11*(uas21*(uas21*uas21+uas31*uas31+uas41)*a01
     1         +rat1*(3.*uas21*uas21+uas31*uas31+uas41)*a11
     1         +3.*uas21*uas41*a21+rat1*uas41*a31)*0.5
     1         +gam4m5*uas41*flu11
c
         rat2 = sqrt(uas42)
         ext2 = min(ra3,max(-ra3,-uas22/rat2))
c
         a01  = alp*(ra3-ext2)
         a11  = alp*(ra3**2-ext2**2)/2.
         a21  = alp*(ra3**3-ext2**3)/3.
         a31  = alp*(ra3**4-ext2**4)/4.
c
         a02  =1.0 -a01
         a12  =-a11
         a22  =1.0 -a21
         a32  =-a31
c
         flu12= uas12*(uas22*a02+rat2*a12)
         flu22= uas12*(uas22*uas22*a02+2.*uas22*rat2*a12
     1                 +uas42*a22)
         flu32= uas32*flu12
         flu42= uas12*(uas22*(uas22*uas22+uas32*uas32+uas42)*a02
     1         +rat2*(3.*uas22*uas22+uas32*uas32+uas42)*a12
     1         +3.*uas22*uas42*a22+rat2*uas42*a32)*0.5
     1         +gam4m5*uas42*flu12
c
c rotation inverse
c
         flu210 = flu21
         flu21  = xnn*flu210-ynn*flu31
         flu31  = ynn*flu210+xnn*flu31
         flu220 = flu22
         flu22  = xnn*flu220-ynn*flu32
         flu32  = ynn*flu220+xnn*flu32
c
         flu11=(flu11 + flu12)*rnn
         flu12=(flu21 + flu22)*rnn
         flu13=(flu31 + flu32)*rnn
         flu14=(flu41 + flu42)*rnn
c
         fluro(nsg)=-flu11
c
         ce(1,nubo1) = ce(1,nubo1) - flu11
         ce(2,nubo1) = ce(2,nubo1) - flu12
         ce(3,nubo1) = ce(3,nubo1) - flu13
         ce(4,nubo1) = ce(4,nubo1) - flu14
c
         ce(1,nubo2) = ce(1,nubo2) + flu11
         ce(2,nubo2) = ce(2,nubo2) + flu12
         ce(3,nubo2) = ce(3,nubo2) + flu13
         ce(4,nubo2) = ce(4,nubo2) + flu14
c
500   continue
c
      return
      end
      subroutine fluosh
c     ---------------------------------------------------------
c     edgewise convective fluxes computation using osher's
c     approximate riemann solver
c     the nodal gradients are computed using a beta-combination
c     of centered and hermitian (half-upwind) gradients
c     bijan mohammadi, INRIA
c     ---------------------------------------------------------
      include 'nsc2ke.inc'
c     ---------------------------------------------------------
c
           real dpm(4),dpor(4),dpex(4),aux1(4),aux2(4),
     1     gradi(4),gradj(4)
      e2       =1.e-16
      beta2 = beta
      beta3 = 0.5*(1.0 - 2.0*beta)
      if(nordre.eq.2) then
      beta2 = 2.0*beta
      beta3 = 1.0 - 2.0*beta
      endif
      gg1      = gam/gam1
      GAMO     = (GAM-1.)*.5
      USG0     = 1./GAMO
      POW      = 1./(2.*GAM)
      COEFF    = GAM1/(GAM+1.)
c
c     loop on global list of edges
c
      do 500 nsg=1,nseg
c
         xnn       =- vnocl(1,nsg)
         ynn       =- vnocl(2,nsg)
         rnn       =  vnocl(3,nsg) * 0.5
c
         nubo1     = nubo(1,nsg)
         nubo2     = nubo(2,nsg)
c
         aix       = coor(1,nubo2)-coor(1,nubo1)
         aiy       = coor(2,nubo2)-coor(2,nubo1)
c
         uas11     = ua(1,nubo1)
         uas12     = ua(1,nubo2)
         uas21     = ua(2,nubo1)
         uas22     = ua(2,nubo2)
         uas31     = ua(3,nubo1)
         uas32     = ua(3,nubo2)
         uas41     = pres(nubo1)
         uas42     = pres(nubo2)
c
         if(nordre.eq.1) goto 1234
c
         gradI(1)  = beta2*(aix*dx(1,nubo1) + aiy*dy(1,nubo1)) +
     &               beta3*(uas12 - uas11)
c
         gradI(2)  = beta2*(aix*dx(2,nubo1) + aiy*dy(2,nubo1)) +
     &               beta3*(uas22 - uas21)
c
         gradI(3)  = beta2*(aix*dx(3,nubo1) + aiy*dy(3,nubo1)) +
     &               beta3*(uas32 - uas31)
c
         gradI(4)  = beta2*(aix*dx(4,nubo1) + aiy*dy(4,nubo1)) +
     &               beta3*(uas42 - uas41)
c
         gradJ(1)  = beta2*(aix*dx(1,nubo2) + aiy*dy(1,nubo2)) +
     &               beta3*(uas12 - uas11)
c
         gradJ(2)  = beta2*(aix*dx(2,nubo2) + aiy*dy(2,nubo2)) +
     &               beta3*(uas22 - uas21)
c
         gradJ(3)  = beta2*(aix*dx(3,nubo2) + aiy*dy(3,nubo2)) +
     &               beta3*(uas32 - uas31)
c
         gradJ(4)  = beta2*(aix*dx(4,nubo2) + aiy*dy(4,nubo2)) +
     &               beta3*(uas42 - uas41)
c
         if(nordre.eq.2) then
c
         uas11     = uas11 + 0.5*gradI(1)
         uas21     = uas21 + 0.5*gradI(2)
         uas31     = uas31 + 0.5*gradI(3)
         uas41     = uas41 + 0.5*gradI(4)
c
         uas12     = uas12 - 0.5*gradJ(1)
         uas22     = uas22 - 0.5*gradJ(2)
         uas32     = uas32 - 0.5*gradJ(3)
         uas42     = uas42 - 0.5*gradJ(4)
c
         elseif(nordre.eq.3) then
c
         dpm(1)    =-(uas12 - uas11)
         dpex(1)   =-4.0*gradI(1) - dpm(1)
         aux1(1)   = 0.25*(1.0+  SIGN(1.0, dpex(1)*dpm(1)))
         dpor(1)   =-4.0*gradJ(1) - dpm(1)
         aux2(1)   = 0.25*(1.0 + SIGN(1.0, dpor(1)*dpm(1)))
c
         dpm(2)    =-(uas22 - uas21)
         dpex(2)   =-4.0*gradI(2) - dpm(2)
         aux1(2)   = 0.25*(1.0 + SIGN(1.0, dpex(2)*dpm(2)))
         dpor(2)   =-4.0*gradJ(2) - dpm(2)
         aux2(2)   = 0.25*(1.0 + SIGN(1.0, dpor(2)*dpm(2)))
c
         dpm(3)    =-(uas32 - uas31)
         dpex(3)   =-4.0*gradI(3) - dpm(3)
         aux1(3)   = 0.25*(1.0 + SIGN(1.0, dpex(3)*dpm(3)))
         dpor(3)   =-4.0*gradJ(3) - dpm(3)
         aux2(3)   = 0.25*(1.0 + SIGN(1.0, dpor(3)*dpm(3)))
c
         dpm(4)    =-(uas42 - uas41)
         dpex(4)   =-4.0*gradI(4) - dpm(4)
         aux1(4)   = 0.25*(1.0 + SIGN(1.0, dpex(4)*dpm(4)))
         dpor(4)   =-4.0*gradJ(4) - dpm(4)
         aux2(4)   = 0.25*(1.0 + SIGN(1.0, dpor(4)*dpm(4)))
c
         gradI(1)  = aux1(1)*
     &               ((dpex(1)*dpex(1) + e2)*dpm(1)+
     &                (dpm(1)*dpm(1)   + e2)*dpex(1))/
     &               (dpex(1)*dpex(1)  + dpm(1)*dpm(1) + 2.0*e2)
         gradJ(1)  = aux2(1)*
     &               ((dpor(1)*dpor(1) + e2)*dpm(1)+
     &                (dpm(1)*dpm(1)   + e2)*dpor(1))/
     &               (dpor(1)*dpor(1)  + dpm(1)*dpm(1) + 2.0*e2)
c
         gradI(2)  = aux1(2)*
     &               ((dpex(2)*dpex(2) + e2)*dpm(2)+
     &                (dpm(2)*dpm(2)   + e2)*dpex(2))/
     &               (dpex(2)*dpex(2)  + dpm(2)*dpm(2) + 2.0*e2)
         gradJ(2)  = aux2(2)*
     &               ((dpor(2)*dpor(2) + e2)*dpm(2)+
     &                (dpm(2)*dpm(2)   + e2)*dpor(2))/
     &               (dpor(2)*dpor(2)  + dpm(2)*dpm(2) + 2.0*e2)
c
         gradI(3)  = aux1(3)*
     &               ((dpex(3)*dpex(3) + e2)*dpm(3)+
     &                (dpm(3)*dpm(3)   + e2)*dpex(3))/
     &               (dpex(3)*dpex(3)  + dpm(3)*dpm(3) + 2.0*e2)
         gradJ(3)  = aux2(3)*
     &               ((dpor(3)*dpor(3) + e2)*dpm(3)+
     &                (dpm(3)*dpm(3)   + e2)*dpor(3))/
     &               (dpor(3)*dpor(3)  + dpm(3)*dpm(3) + 2.0*e2)
c
         gradI(4)  = aux1(4)*
     &               ((dpex(4)*dpex(4) + e2)*dpm(4)+
     &                (dpm(4)*dpm(4)   + e2)*dpex(4))/
     &               (dpex(4)*dpex(4)  + dpm(4)*dpm(4) + 2.0*e2)
         gradJ(4)  = aux2(4)*
     &               ((dpor(4)*dpor(4) + e2)*dpm(4)+
     &                (dpm(4)*dpm(4)   + e2)*dpor(4))/
     &               (dpor(4)*dpor(4)  + dpm(4)*dpm(4) + 2.0*e2)
c
         uas11     = uas11 - gradI(1)
         uas21     = uas21 - gradI(2)
         uas31     = uas31 - gradI(3)
         uas41     = uas41 - gradI(4)
c
         uas12     = uas12 + gradJ(1)
         uas22     = uas22 + gradJ(2)
         uas32     = uas32 + gradJ(3)
         uas42     = uas42 + gradJ(4)
c
         endif
c
1234   continue
c
      PROD1 = UAS21*xnn+UAS31*ynn
      C1    = SQRT(GAM*UAS41/UAS11)
      PROD2 = UAS22*xnn+UAS32*ynn
      C2    = SQRT(GAM*UAS42/UAS12)
C
C    -- FLUXOSHER = F(NUBO1)+F(NUBO2)
C
      FM1 = PROD1*UAS11+PROD2*UAS12
      FM2 = UAS11*UAS21*PROD1+UAS41*xnn+
     &      UAS12*UAS22*PROD2+UAS42*xnn
      FM3 = UAS11*UAS31*PROD1+UAS41*ynn+
     &      UAS12*UAS32*PROD2+UAS42*ynn
      FM4 = (GG1*UAS41+.5*UAS11*
     & (UAS21*UAS21+UAS31*UAS31))*PROD1+
     &      (GG1*UAS42+.5*UAS12*
     & (UAS22*UAS22+UAS32*UAS32))*PROD2
C
C  II/ CALCUL DES QUANTITES AUX POINTS LIMITES POUR LE C.L.D
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DNUM       = GAMO*(PROD2-PROD1)+C1+C2
      DEN        = ((UAS41/UAS42)**POW)*SQRT(UAS12/UAS11)
      ULD11      = ((DNUM/((1.+1./DEN)*C1))**USG0)*UAS11
      ULD21      = ((DNUM/((1.+   DEN)*C2))**USG0)*UAS12
      ULD4       = ((ULD11/UAS11)**GAM)*UAS41
C
      AUX        = USG0*(C1-SQRT(GAM*ULD4/ULD11))
      ULD12      = UAS21-AUX*xnn
      ULD13      = UAS31-AUX*ynn
      ULD22      = UAS22-(PROD2-PROD1+AUX)*xnn
      ULD23      = UAS32-(PROD2-PROD1+AUX)*ynn
C
C  III/ CALCUL DES FLUX AUX POINTS LIMITES DU C.L.D
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     -- FLULD = FLUX LINEAIREMENT DEGENERE
c
      PRODLD       = ULD12*xnn+ULD13*ynn
c
      FLULD11      = ULD11*PRODLD
      FLULD12      = ULD11*ULD12*PRODLD+ULD4*xnn
      FLULD13      = ULD11*ULD13*PRODLD+ULD4*ynn
      FLULD14      = (gg1*ULD4+.5*(ULD12*ULD12+ulD13*ULD13)*ULD11)
     & *PRODLD
c
      FLULD21      = ULD21*PRODLD
      FLULD22      = ULD21*ULD22*PRODLD+ULD4*xnn
      FLULD23      = ULD21*ULD23*PRODLD+ULD4*ynn
      FLULD24      = (gg1*ULD4+.5*(ULD22*ULD22+ULD23*ULD23)*ULD21)
     & *PRODLD
c
C                           -- FLUOSHER = FLUOSHER - FLULD
c
      AUX     =ABS(PRODLD)
      FM1=FM1-AUX*(ULD21-ULD11)
      FM2=FM2-AUX*(ULD21*ULD22-ULD11*ULD12)
      FM3=FM3-AUX*(ULD21*ULD23-ULD11*ULD13)
      FM4=FM4-AUX*.5*(ULD21*(ULD22*ULD22+ULD23*ULD23)-
     &                ULD11*(ULD12*ULD12+ULD13*ULD13))
C
C  IV/ CALCUL DES QUANTITES AUX POINTS LIMITES POUR LES 2 C.V.N.L.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      AUX        = (PROD2+C2*USG0)/C2
      AX         = .5*(1.+SIGN(1.,AUX))
      UVNL21     = AX*((COEFF*ABS(AUX))**USG0)*UAS12
      UVNL24     = AX*((UVNL21/UAS12)**GAM)*UAS42
      AUX        = USG0*(C2-SQRT(GAM*UVNL24/(UVNL21+(1.-AX))))
      UVNL22     = AX*(UAS22+AUX*xnn)
      UVNL23     = AX*(UAS32+AUX*ynn)
C
      AUX        = (C1*USG0-PROD1)/C1
      AX         = .5*(1.+SIGN(1.,AUX))
      UVNL11     = AX*((COEFF*ABS(AUX))**USG0)*UAS11
      UVNL14     = AX*((UVNL11/UAS11)**GAM)*UAS41
      AUX        = USG0*(C1-SQRT(GAM*UVNL14/(UVNL11+(1.-AX))))
      UVNL12     = AX*(UAS21-AUX*xnn)
      UVNL13     = AX*(UAS31-AUX*ynn)
c                         -- FLUVN = FLUX VRAIMENT NON LINEAIRE
      PROVL1     = UVNL12*xnn+UVNL13*ynn
      FLUVN11    = UVNL11*PROVL1
      FLUVN12    = UVNL11*UVNL12*PROVL1+UVNL14*xnn
      FLUVN13    = UVNL11*UVNL13*PROVL1+UVNL14*ynn
      FLUVN14    = (gg1*UVNL14+.5*(UVNL12*UVNL12+
     &             UVNL13*UVNL13)*UVNL11)*PROVL1
C
      PROVL2     = UVNL22*xnn+UVNL23*ynn
      FLUVN21    = UVNL21*PROVL2
      FLUVN22    = UVNL21*UVNL22*PROVL2+UVNL24*xnn
      FLUVN23    = UVNL21*UVNL23*PROVL2+UVNL24*ynn
      FLUVN24    = (gg1*UVNL24+.5*(UVNL22*UVNL22+
     &             UVNL23*UVNL23)*UVNL21)*PROVL2
C
      SI1        = SIGN(1.,PROD2 -C2)
      SI2        = SIGN(1.,PRODLD-SQRT(GAM*ULD4/ULD21))
      SI3        = SIGN(1.,PRODLD+SQRT(GAM*ULD4/ULD11))
      SI4        = SIGN(1.,PROD1 +C1)
c
C       -- FLUOSHER = FLUOSHER - FLUVN
c
      FM1=FM1-(UAS12*PROD2-FLUVN21)*SI1
     &       +(  FLULD21-FLUVN21)*SI2
     &       -(  FLULD11-FLUVN11)*SI3
     &       +(  UAS11*PROD1-FLUVN11)*SI4
      FM2=FM2-(UAS12*UAS22*PROD2+
     &         UAS42*xnn-FLUVN22)*SI1
     &       +(  FLULD22-FLUVN22)*SI2
     &       -(  FLULD12-FLUVN12)*SI3
     &       +(  UAS11*UAS21*PROD1+
     &           UAS41*xnn-FLUVN12)*SI4
      FM3=FM3-(UAS12*UAS32*PROD2+
     &         UAS42*ynn-FLUVN23)*SI1
     &       +(  FLULD23-FLUVN23)*SI2
     &       -(  FLULD13-FLUVN13)*SI3
     &       +(  UAS11*UAS31*PROD1+
     &           UAS41*ynn-FLUVN13)*SI4
      FM4=FM4-((gg1*UAS42+.5*(UAS22*UAS22+UAS32*UAS32)*
     &        UAS12)*PROD2-FLUVN24)*SI1
     &       +(    FLULD24-FLUVN24)*SI2
     &       -(    FLULD14-FLUVN14)*SI3
     &       +((gg1*UAS41+.5*(UAS21*UAS21+UAS31*UAS31)*
     &        UAS11)*PROD1-FLUVN14)*SI4
C
C 3/ REPORT SUR LES NOEUDS DES SEGMENTS
C======================================
C
      fluro(nsg)=-fm1*rnn
c
      CE(1,NUBO1)=CE(1,NUBO1)-FM1*rnn
      CE(2,NUBO1)=CE(2,NUBO1)-FM2*rnn
      CE(3,NUBO1)=CE(3,NUBO1)-FM3*rnn
      CE(4,NUBO1)=CE(4,NUBO1)-FM4*rnn
C
      CE(1,NUBO2)=CE(1,NUBO2)+FM1*rnn
      CE(2,NUBO2)=CE(2,NUBO2)+FM2*rnn
      CE(3,NUBO2)=CE(3,NUBO2)+FM3*rnn
      CE(4,NUBO2)=CE(4,NUBO2)+FM4*rnn
c
500   continue
c
      return
      end
      subroutine fluroe
c     -----------------------------------------------------------------
c     edgewise convective fluxes computation using roe's
c     approximate riemann solver
c     the nodal gradients are computed using a beta-combination
c     of centered and hermitian (half-upwind) gradients
c     -----------------------------------------------------------------
      include 'nsc2ke.inc'
c     -----------------------------------------------------------------
c     local variables definition
      integer  nsg    , nubo1  , nubo2
      real invas1, invcnt
      real dpm(4),dpor(4),dpex(4),aux1(4),aux2(4),
     1     gradi(4),gradj(4)
c
      usgam1=1.0/gam1
      gsgam1=gam/gam1
      beta2 = beta
      beta3 = 0.5*(1.0 - 2.0*beta)
      if(nordre.eq.2) then
      beta2 = 2.0*beta
      beta3 = 1.0 - 2.0*beta
      endif
      e2 = 1.e-16
c
c     loop on global list of edges
c
      do 500 nsg=1,nseg
c
c        local indexing of the vertices of the current edge
c
         nubo1     = nubo(1,nsg)
         nubo2     = nubo(2,nsg)
c
         aix       = coor(1,nubo2)-coor(1,nubo1)
         aiy       = coor(2,nubo2)-coor(2,nubo1)
c
c        indirect addressing on vertices physical states
c
         uas11     = ua(1,nubo1)
         uas21     = ua(2,nubo1)
         uas31     = ua(3,nubo1)
         uas41     = pres(nubo1)
         uas12     = ua(1,nubo2)
         uas22     = ua(2,nubo2)
         uas32     = ua(3,nubo2)
         uas42     = pres(nubo2)
c
         if(nordre.eq.1) goto 1234
c
         gradI(1)  = beta2*(aix*dx(1,nubo1) + aiy*dy(1,nubo1)) +
     &               beta3*(uas12 - uas11)
c
         gradI(2)  = beta2*(aix*dx(2,nubo1) + aiy*dy(2,nubo1)) +
     &               beta3*(uas22 - uas21)
c
         gradI(3)  = beta2*(aix*dx(3,nubo1) + aiy*dy(3,nubo1)) +
     &               beta3*(uas32 - uas31)
c
         gradI(4)  = beta2*(aix*dx(4,nubo1) + aiy*dy(4,nubo1)) +
     &               beta3*(uas42 - uas41)
c
         gradJ(1)  = beta2*(aix*dx(1,nubo2) + aiy*dy(1,nubo2)) +
     &               beta3*(uas12 - uas11)
c
         gradJ(2)  = beta2*(aix*dx(2,nubo2) + aiy*dy(2,nubo2)) +
     &               beta3*(uas22 - uas21)
c
         gradJ(3)  = beta2*(aix*dx(3,nubo2) + aiy*dy(3,nubo2)) +
     &               beta3*(uas32 - uas31)
c
         gradJ(4)  = beta2*(aix*dx(4,nubo2) + aiy*dy(4,nubo2)) +
     &               beta3*(uas42 - uas41)
c
         if(nordre.eq.2) then
c
         uas11     = uas11 + 0.5*gradI(1)
         uas21     = uas21 + 0.5*gradI(2)
         uas31     = uas31 + 0.5*gradI(3)
         uas41     = uas41 + 0.5*gradI(4)
c
         uas12     = uas12 - 0.5*gradJ(1)
         uas22     = uas22 - 0.5*gradJ(2)
         uas32     = uas32 - 0.5*gradJ(3)
         uas42     = uas42 - 0.5*gradJ(4)
c
         elseif(nordre.eq.3) then
c
         dpm(1)    =-(uas12 - uas11)
         dpex(1)   =-4.0*gradI(1) - dpm(1)
         aux1(1)   = 0.25*(1.0+  SIGN(1.0, dpex(1)*dpm(1)))
         dpor(1)   =-4.0*gradJ(1) - dpm(1)
         aux2(1)   = 0.25*(1.0 + SIGN(1.0, dpor(1)*dpm(1)))
c
         dpm(2)    =-(uas22 - uas21)
         dpex(2)   =-4.0*gradI(2) - dpm(2)
         aux1(2)   = 0.25*(1.0 + SIGN(1.0, dpex(2)*dpm(2)))
         dpor(2)   =-4.0*gradJ(2) - dpm(2)
         aux2(2)   = 0.25*(1.0 + SIGN(1.0, dpor(2)*dpm(2)))
c
         dpm(3)    =-(uas32 - uas31)
         dpex(3)   =-4.0*gradI(3) - dpm(3)
         aux1(3)   = 0.25*(1.0 + SIGN(1.0, dpex(3)*dpm(3)))
         dpor(3)   =-4.0*gradJ(3) - dpm(3)
         aux2(3)   = 0.25*(1.0 + SIGN(1.0, dpor(3)*dpm(3)))
c
         dpm(4)    =-(uas42 - uas41)
         dpex(4)   =-4.0*gradI(4) - dpm(4)
         aux1(4)   = 0.25*(1.0 + SIGN(1.0, dpex(4)*dpm(4)))
         dpor(4)   =-4.0*gradJ(4) - dpm(4)
         aux2(4)   = 0.25*(1.0 + SIGN(1.0, dpor(4)*dpm(4)))
c
         gradI(1)  = aux1(1)*
     &               ((dpex(1)*dpex(1) + e2)*dpm(1)+
     &                (dpm(1)*dpm(1)   + e2)*dpex(1))/
     &               (dpex(1)*dpex(1)  + dpm(1)*dpm(1) + 2.0*e2)
         gradJ(1)  = aux2(1)*
     &               ((dpor(1)*dpor(1) + e2)*dpm(1)+
     &                (dpm(1)*dpm(1)   + e2)*dpor(1))/
     &               (dpor(1)*dpor(1)  + dpm(1)*dpm(1) + 2.0*e2)
c
         gradI(2)  = aux1(2)*
     &               ((dpex(2)*dpex(2) + e2)*dpm(2)+
     &                (dpm(2)*dpm(2)   + e2)*dpex(2))/
     &               (dpex(2)*dpex(2)  + dpm(2)*dpm(2) + 2.0*e2)
         gradJ(2)  = aux2(2)*
     &               ((dpor(2)*dpor(2) + e2)*dpm(2)+
     &                (dpm(2)*dpm(2)   + e2)*dpor(2))/
     &               (dpor(2)*dpor(2)  + dpm(2)*dpm(2) + 2.0*e2)
         gradI(3)  = aux1(3)*
     &               ((dpex(3)*dpex(3) + e2)*dpm(3)+
     &                (dpm(3)*dpm(3)   + e2)*dpex(3))/
     &               (dpex(3)*dpex(3)  + dpm(3)*dpm(3) + 2.0*e2)
         gradJ(3)  = aux2(3)*
     &               ((dpor(3)*dpor(3) + e2)*dpm(3)+
     &                (dpm(3)*dpm(3)   + e2)*dpor(3))/
     &               (dpor(3)*dpor(3)  + dpm(3)*dpm(3) + 2.0*e2)
c
         gradI(4)  = aux1(4)*
     &               ((dpex(4)*dpex(4) + e2)*dpm(4)+
     &                (dpm(4)*dpm(4)   + e2)*dpex(4))/
     &               (dpex(4)*dpex(4)  + dpm(4)*dpm(4) + 2.0*e2)
         gradJ(4)  = aux2(4)*
     &               ((dpor(4)*dpor(4) + e2)*dpm(4)+
     &                (dpm(4)*dpm(4)   + e2)*dpor(4))/
     &               (dpor(4)*dpor(4)  + dpm(4)*dpm(4) + 2.0*e2)
c
         uas11     = uas11 - gradI(1)
         uas21     = uas21 - gradI(2)
         uas31     = uas31 - gradI(3)
         uas41     = uas41 - gradI(4)
c
         uas12     = uas12 + gradJ(1)
         uas22     = uas22 + gradJ(2)
         uas32     = uas32 + gradJ(3)
         uas42     = uas42 + gradJ(4)
c
         endif
c
1234   continue
c
         invas1    = 1.0/uas11
c
         dua1      = sqrt(uas11)
         dua2      = sqrt(uas12)
         dua3      = 1.0/(dua1 + dua2)
c
c        enthalpy
c
         h1        = gsgam1*uas41*invas1 +
     &               0.5*(uas21*uas21 + uas31*uas31)
         h2        = gsgam1*uas42/uas12  +
     &               0.5*(uas22*uas22 + uas32*uas32)
c
c        roe's mean values computation
c
         ucent     = (uas21*dua1 + uas22*dua2)*dua3
         vcent     = (uas31*dua1 + uas32*dua2)*dua3
         hcent     = (h1*dua1 + h2*dua2)*dua3
         uvnc      = 0.5*(ucent*ucent + vcent*vcent)
         ccent     = sqrt(abs(gam1*(hcent-uvnc)))
         invcnt    = 1.0/ccent
c
         pres1     = uas41
         uas41     = usgam1*uas41 + 0.5*uas11*
     &               (uas21*uas21 + uas31*uas31)
         uas42     = usgam1*uas42 + 0.5*uas12*
     &               (uas22*uas22 + uas32*uas32)
         uas21     = uas21*uas11
         uas22     = uas22*uas12
         uas31     = uas31*uas11
         uas32     = uas32*uas12
c
         dua1      = uas12 - uas11
         dua2      = uas22 - uas21
         dua3      = uas32 - uas31
         dua4      = uas42 - uas41
c
c        eigenvalues computation
c
         xnn       = vnocl(1,nsg)
         ynn       = vnocl(2,nsg)
         rnn       = vnocl(3,nsg)
c
         vp1       = rnn*(xnn*ucent + ynn*vcent)
         vp3       = vp1 + ccent*rnn
         vp4       = vp1 - ccent*rnn
         uq41      = uvnc*invcnt + usgam1*ccent
         uq42      = xnn*ucent + ynn*vcent
         rc1       = gam1*invcnt
         rc2       = rc1*invcnt
c
c        computation of the centered part of the flux
c
         flu11     = rnn*(xnn*uas21 + ynn*uas31)
         flu21     = xnn*rnn*pres1 + flu11*uas21*invas1
         flu31     = ynn*rnn*pres1 + flu11*uas31*invas1
         flu41     = h1*flu11
c
c        computation of the diffusive part of the flux
c
         fdc1      = max(vp1,0.0)*
     &               (dua1 + rc2*(-uvnc*dua1 + ucent*dua2 +
     &                            vcent*dua3 - dua4))
         fdc2      = max(vp1,0.0)*
     &               ((xnn*vcent - ynn*ucent)*
     &                dua1 + ynn*dua2 - xnn*dua3)
         fdc3      = max(vp3,0.0)*
     &               (0.5*(-uq42*dua1 + xnn*dua2 + ynn*dua3) +
     &                0.5*rc1*(uvnc*dua1 - ucent*dua2 - vcent*dua3 +
     &                dua4))
         fdc4      = max(vp4,0.0)*
     &               (0.5*( uq42*dua1 - xnn*dua2 - ynn*dua3) +
     &                0.5*rc1*(uvnc*dua1 - ucent*dua2 - vcent*dua3 +
     &                dua4))
c
         duv1      = fdc1 + (fdc3 + fdc4)*invcnt
         duv2      = ucent*fdc1    + ynn*fdc2  +
     &               (ucent*invcnt + xnn)*fdc3 +
     &               (ucent*invcnt - xnn)*fdc4
         duv3      = vcent*fdc1    - xnn*fdc2  +
     &               (vcent*invcnt + ynn)*fdc3 +
     &               (vcent*invcnt - ynn)*fdc4
         duv4      = uvnc*fdc1 + (ynn*ucent - xnn*vcent)* fdc2 +
     &               (uq41+uq42)*fdc3 + (uq41-uq42)*fdc4
c
c        gathering of the elementary fluxes into the global ones
c
c
         fluro(nsg)=flu11+duv1
c
         ce(1,nubo1)    = ce(1,nubo1) + flu11 + duv1
         ce(2,nubo1)    = ce(2,nubo1) + flu21 + duv2
         ce(3,nubo1)    = ce(3,nubo1) + flu31 + duv3
         ce(4,nubo1)    = ce(4,nubo1) + flu41 + duv4
c
         ce(1,nubo2)    = ce(1,nubo2) - flu11 - duv1
         ce(2,nubo2)    = ce(2,nubo2) - flu21 - duv2
         ce(3,nubo2)    = ce(3,nubo2) - flu31 - duv3
         ce(4,nubo2)    = ce(4,nubo2) - flu41 - duv4
c
500   continue
c
      return
      end
       subroutine geowall
       include 'nsc2ke.inc'
c
c maxt  : nombre maximal de tetraedres contenant une normale
c maxp  : nombre maximal de points pour le calcul de la contrainte
c         parietale
c nbt : nombre de tetraedre interceptant la normale en if11
c jtb : numero du tetraedre solution
c nbp : nombre de points pour le calcul de la contrainte en if11
c ipb : numero des points solution
c
      integer deja(nn)
      integer inb,is,is1,is2,jt,ip1(3),ip2(3),nb,isb
      integer k,if11,j
c
      data ip1 /2,3,1/
      data ip2 /3,1,2/
c
      do if11=1,nslog3
         nbt(if11)=0
         do inb=1,maxt
            jtb(inb,if11)=-1
         enddo
      enddo
c
      do 100 jt=1,nt
c
         do 50 k=1,3
c
            is=nu(k,jt)
            is1=nu(ip1(k),jt)
            is2=nu(ip2(k),jt)
c
c si is est un point adherent
c ...........................
c
            if(numc(is).gt.0) then
c
            if11=numc(is)
c
               if(nbt(if11).lt.maxt) then
c
c un nouveau tetraedre est solution
c
                  nbt(if11)=nbt(if11)+1
                  jtb(nbt(if11),if11)=jt
c
               endif
         endif
c
50    continue
100   continue
c
c *** points pour le calcul de la vitesse de frottement *****
c ...........................................................
c
      do 110 is=1,ns
         deja(is)=-1
110   continue
c
      do 120 if11=1,nslog3
         nbp(if11)=0
120   continue
c
c boucle sur les points adherents
c
      do 500 if11=1,nslog3
c
         do 400 nb=1,nbt(if11)
c
            do 300 k=1,3
               isb=nu(k,jtb(nb,if11))
c
               if(numc(isb).eq.0.and.deja(isb).lt.0) then
                   nbp(if11)=nbp(if11)+1
                   ipb(nbp(if11),if11)=isb
               endif
c
               if(nbt(if11).gt.1) deja(isb)=1
c
               if(nbt(if11).eq.nb) deja(isb)=-1
c
300         continue
400      continue
c
      do 450 is=1,ns
         deja(is)=-1
450   continue
500   continue
c
       do is=1,ns
       dist(is)=1.e10
       x1=coor(1,is)
       y1=coor(2,is)
       do j=1,nslog3
       js=log3(j)
       x2=coor(1,js)
       y2=coor(2,js)
       dd=sqrt((x1-x2)**2+(y1-y2)**2)
       if(dd.lt.dist(is)) then
        dist(is)=dd
        nswall(is)=js
       endif
       enddo
       enddo
c
      return
      end
      subroutine gravity
c
cc gravity forces
c bijan mohammadi, INRIA
c
      include 'nsc2ke.inc'
c
      do 100 is=1,ns
      air=airsa(is)
      ce(3,is)=ce(3,is)-air*froud*ua(1,is)
100   continue
c
          return
          end
      subroutine init(kt0, t0)
c     -----------------------------------------------------------------
c     physical solution initialization
c     -----------------------------------------------------------------
      include 'nsc2ke.inc'
c     -----------------------------------------------------------------
c
         kt0=0
         t0=0.0
c
         do 200 is=1,ns
            ua(1,is)    = roin
            ua(2,is)    = roin*uxin
            ua(3,is)    = roin*uyin
            ua(4,is)    = ein
            ua(5,is)    = xkin
            ua(6,is)    = xein
            reyturb(is) = 0.0
            if(iturb.ne.0) reyturb(is) = xkin*xkin*0.09/xein
200      continue
c
c
c     reinitialization strategy
c
         if(ncont.eq.1) then
         open(20,file="INIT_NS",status="unknown")
         do is=1,ns
c
         read(20,*) (ua(i,is),i=1,4)
c
         enddo
          read(20,*) kt0,t0
         close(20)
         endif
         if(ncont1.eq.1) then
         open(21,file="INIT_KE",status="unknown")
         do is=1,ns
         read(21,*) ua(5,is),ua(6,is),reytot,reyturb(is)
         reylam(is)=abs(reytot-reyturb(is))
         enddo
         close(21)
         endif
c
      if (ivis.eq.1.and.ilaw.eq.0) then
         do 9001 i=1,nslog3
            is=log3(i)
            ro=ua(1,is)
            ru=ua(2,is)
            rv=ua(3,is)
            ua(2,is) = 0.0
            ua(3,is) = 0.0
            if(iecc.eq.2) then
              ua(4,is) = ua(1,is)*tbrd2
                          else
              ua(4,is) = ua(4,is)-0.5*(ru**2+rv**2)/ro
              if(ua(4,is).lt.0.) stop 'init ua(4,is)<0'
            endif
            ua(5,is) = 0.0
9001     continue
      endif
c
      return
      end
      SUBROUTINE isovat(ifile,F,COOR,NVAL,VAL,x0,y0,x1,y1)
C
C       TRACE DES LIGNES ISOVALEURS D'UNE FONCTION LINEAIRE
C         SUR UN TRIANGLE
C       IOPFEN=1 RESTRICTION A UNE FENETRE
C       IOPFEN=0 PAS DE RESTRICTION
C       IOPFEN=-1 RESTRICTION A L'EXTERIEUR D'UNE FENETRE
C  *************************************************************
C
      DIMENSION F(3),COOR(2,3),VAL(100)
      DIMENSION IP1(3)
C
      epsi=1.e-5
      IP1(1)=2
      IP1(2)=3
      IP1(3)=1
      FF1=F(1)
      FF2=F(2)
      FF3=F(3)
      FFMA=MAX(ABS(FF1),ABS(FF2))
      FFMA=MAX(ffma,ABS(FF3))
      D12=ABS(FF1-FF2)
      D23=ABS(FF2-FF3)
      IF(D12+D23.LT.AMAX1(epsi,epsi*FFMA)) GOTO 1000
C  PAS DE RESTRICTION
C  ******************
      DO 100 IVAL=1,NVAL
      VAL1=VAL(IVAL)
      ITR=0
      DO 110 K=1,3
      FK=F(K)
      FK1=F(IP1(K))
      FMI=MIN(FK,FK1)
      FMA=MAX(FK,FK1)
      DIF=FMA-FMI
      IF(DIF.LT.epsi) GOTO 110
      EPS=epsi*DIF
      IF(VAL1.LT.FMI-EPS.OR.VAL1.GT.FMA+EPS) GOTO 110
      HH=ABS(FK-VAL1)/DIF
      X=COOR(1,K)+HH*(COOR(1,IP1(K))-COOR(1,K))
      Y=COOR(2,K)+HH*(COOR(2,IP1(K))-COOR(2,K))
      IF(ITR.EQ.0) GOTO 115
      write(ifile,*) x,y
      write(ifile,*)
      GOTO 100
115   ITR=1
      write(ifile,*) x,y
110   CONTINUE
100   CONTINUE
1000   return
      END
      subroutine keps2d
c     -----------------------------------------------------------------
      include 'nsc2ke.inc'
c     -----------------------------------------------------------------
c     local variables definition
      real    dbxx(3)  , dbyy(3), dxt(2), dyt(2),
     &        uph(2,3)
c
      us6               = 1.0/6.0
      us3               = 1.0/3.0
      us43              = 1.0/4.0/3.0
      gampr             = gam/pr
      gamprt            = gam/prt
      cmu               = 0.09
      cep               = 1./1.3
c
c     loop on global list of triangles
c
      do 1000 jt=1,nt
c
c        scatter operation
c        getting the physical solutions for the three vertices
c        of the current triangle
c
         nubo1          = nu(1,jt)
         nubo2          = nu(2,jt)
         nubo3          = nu(3,jt)
c
         uph(1,1)       = ua(5,nubo1)
         uph(2,1)       = ua(6,nubo1)
c
         uph(1,2)       = ua(5,nubo2)
         uph(2,2)       = ua(6,nubo2)
c
         uph(1,3)       = ua(5,nubo3)
         uph(2,3)       = ua(6,nubo3)
c
c        computation of the p1-gradients
c
         x1             = coor(1,nubo1)
         y1             = coor(2,nubo1)
         x2             = coor(1,nubo2)
         y2             = coor(2,nubo2)
         x3             = coor(1,nubo3)
         y3             = coor(2,nubo3)
c
         dbxx(1)        = y2 - y3
         dbxx(2)        = y3 - y1
         dbxx(3)        = y1 - y2
         dbyy(1)        = x3 - x2
         dbyy(2)        = x1 - x3
         dbyy(3)        = x2 - x1
c
            dxt(1)      = uph(1,1)*dbxx(1) +
     &                    uph(1,2)*dbxx(2) +
     &                    uph(1,3)*dbxx(3)
            dyt(1)      = uph(1,1)*dbyy(1) +
     &                    uph(1,2)*dbyy(2) +
     &                    uph(1,3)*dbyy(3)
            dxt(2)      = uph(2,1)*dbxx(1) +
     &                    uph(2,2)*dbxx(2) +
     &                    uph(2,3)*dbxx(3)
            dyt(2)      = uph(2,1)*dbyy(1) +
     &                    uph(2,2)*dbyy(2) +
     &                    uph(2,3)*dbyy(3)
c
       xmlam=us3*(reylam(nubo1)+reylam(nubo2)+reylam(nubo3))
       xmtur=us3*(reyturb(nubo1)+reyturb(nubo2)+reyturb(nubo3))
       xmtot=xmlam+xmtur
c
         aitt           = 0.25*airta(jt)
c
         ce(5,nubo1)    = ce(5,nubo1) - aitt*xmtot*
     &                    (dbxx(1)*dxt(1) + dbyy(1)*dyt(1))
         ce(6,nubo1)    = ce(6,nubo1) - aitt*(xmlam+cep*xmtur)*
     &                    (dbxx(1)*dxt(2) + dbyy(1)*dyt(2))
c
         ce(5,nubo2)    = ce(5,nubo2) - aitt*xmtot*
     &                    (dbxx(2)*dxt(1) + dbyy(2)*dyt(1))
         ce(6,nubo2)    = ce(6,nubo2) - aitt*(xmlam+cep*xmtur)*
     &                    (dbxx(2)*dxt(2) + dbyy(2)*dyt(2))
c
         ce(5,nubo3)    = ce(5,nubo3) - aitt*xmtot*
     &                    (dbxx(3)*dxt(1) + dbyy(3)*dyt(1))
         ce(6,nubo3)    = ce(6,nubo3) - aitt*(xmlam+cep*xmtur)*
     &                    (dbxx(3)*dxt(2) + dbyy(3)*dyt(2))
c
1000  continue
c
       do 500 nsg=1,nseg
         nubo1     = nubo(1,nsg)
         nubo2     = nubo(2,nsg)
c
         uas51     = ua(5,nubo1)
         uas61     = ua(6,nubo1)
         uas52     = ua(5,nubo2)
         uas62     = ua(6,nubo2)
c
         sgn=-sign(.5,fluro(nsg))
c
         ce(5,nubo1)    = ce(5,nubo1) + fluro(nsg)*
     1                   ( (.5+sgn)*uas51+(.5-sgn)*uas52 )
         ce(6,nubo1)    = ce(6,nubo1) + fluro(nsg)*
     1                   ( (.5+sgn)*uas61+(.5-sgn)*uas62 )
         ce(5,nubo2)    = ce(5,nubo2) - fluro(nsg)*
     1                   ( (.5+sgn)*uas51+(.5-sgn)*uas52 )
         ce(6,nubo2)    = ce(6,nubo2) - fluro(nsg)*
     1                   ( (.5+sgn)*uas61+(.5-sgn)*uas62 )
c
500     continue
         if(ilaw.ne.0) call ke_law
         if(ilaw.eq.0) call ke_two
c
      return
      end
      subroutine ke_law
      include 'nsc2ke.inc'
c
      cmu=0.09
      c1=0.1296
      c2=1.83333
      xmact0=0.25
c
      do 200 is=1,ns
      x1=coor(1,is)
      x2=coor(2,is)
      reyturb(is)=0.0
      if(x1.lt.xtmin.or.x1.gt.xtmax.or.
     1   x2.lt.ytmin.or.x2.gt.ytmax) goto 200
      ro=ua(1,is)
      xk=ua(5,is)
      xes=ua(6,is)
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc compressibility corrections
cccccccccccccccccccccccccccccccccccccccccccccccccccc
      xmactu=sqrt( 2.*ro*xk / (gam*pres(is)) )
      xmactu2=xmactu**2
      xtest=xmactu-xmact0
      heavyside=1.
      if(xtest.le.0.) heavyside=0.
      fcomp=heavyside*(1.-exp(-(xtest/0.66)**2))
      xed=xes*(1.+fcomp)
cccccccccccccccccccccccccccccccccccccccccccccccccccc
      reyturb(is)=cmu*ro*xk*xk/xes
      enorm=(dx(3,is)+dy(2,is))**2
c
      sk=enorm*reyturb(is)-ro*xed
      se=c1*enorm*ro*xk-c2*ro*xed**2/xk
c
      ce(5,is)=ce(5,is)+airsa(is)*sk
      ce(6,is)=ce(6,is)+airsa(is)*se
c
200   continue
c
      return
      end
      subroutine ke_two
      include 'nsc2ke.inc'
c
      cmu=0.09
      c1=0.1296
      c2=1.83333
      xkap=0.41
      qs3=4./3.
      ds3=2./3.
      xmact0=0.25
      cl=xkap*cmu**(-3./4.)
      if(xmach.lt.1.) xclmu=0.0142
      if(xmach.ge.1.) xclmu=0.0085
c
      do 300 is=1,ns
      ro=ua(1,is)
      y=dist(is)
      iwall=nswall(is)
      x1=coor(1,is)
      x2=coor(2,is)
      xk=ua(5,is)
      xes=ua(6,is)
      enorm=(dx(3,is)+dy(2,is))**2
      divs=dx(2,is)+dy(3,is)
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc compressibility corrections
cccccccccccccccccccccccccccccccccccccccccccccccccccc
      xmactu=sqrt( 2.*ro*xk / (gam*pres(is)) )
      xmactu2=xmactu**2
      xtest=xmactu-xmact0
      heavyside=1.
      if(xtest.le.0.) heavyside=0.
      fzeman=heavyside*(1.-exp(-(xtest/0.66)**2))
      fcomp=fzeman
      xed=xes*(1.+fcomp)
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      xmu=reylam(iwall)
      yplus=sqrt(xk*ua(1,iwall)*ro)*y/xmu
c
      if(yplus.gt.200..or.y.gt.delta) then
      reyturb(is)=0.0
ccc high reynolds region
      if(x1.lt.xtmin.or.x1.gt.xtmax.or.
     1   x2.lt.ytmin.or.x2.gt.ytmax) goto 300
c
      reyturb(is)=cmu*ro*xk*xk/xes
      sk=enorm*reyturb(is)-ro*xed
      se=c1*enorm*ro*xk-c2*ro*xed*xed/xk
      else
ccc near wall region
      reyturb(is)=0.0
      if(logfr(is).ne.0.or.x1.lt.xtmin.or.x1.gt.xtmax.or.
     1   x2.lt.ytmin.or.x2.gt.ytmax) goto 300
      xlmu=cl*y*(1.-exp(-yplus*xclmu))
      xlep=cl*y*(1.-exp(-yplus/(2.0*cl)))
      reyturb(is)=cmu*xlmu*ro*sqrt(xk)
      un(6,is)=ro*xk**1.5/xlep
      sk=enorm*reyturb(is)-ro*xed
      se=0.0
      ce(6,is)=0.0
      endif
c
      ce(5,is)=ce(5,is)+sk*airsa(is)
      ce(6,is)=ce(6,is)+se*airsa(is)
c
300   continue
c
      return
      end
      subroutine loglaw
      include 'nsc2ke.inc'
c
c classical wall laws technique
c        weak form
c
      rklaw  =0.419
      claw   =5.445
      nitmax =500
      eps1   =1.e-5
      cmu    =0.09
c
      do 1000 isc=1,nslog3
c
        is=log3(isc)
        ro=ua(1,is)
        temp=pres(is)/(gam1*ro)
        rey=reylam(is)
        reytot=reylam(is)+reyturb(is)
c
c on norme le vecteur normal entrant
c
        rnorm=vno(is)
        xn1=-vnox(is)
        xn2=-vnoy(is)
        xt1= xn2
        xt2=-xn1
c
c calcul de la vitesse tangentielle
c
      u=ua(2,is)/ro
      v=ua(3,is)/ro
      utang=xt1*u+xt2*v
c
c methode de newton
c
        uf=sqrt(rey*abs(utang)/delta)
        yplus=ro*delta*uf/rey
c
        if(yplus.lt.5.0) then
           uf=rey*5.0/delta/ro
        endif
c
        if(yplus.gt.11.6) then
c
           nit=0
500     uplus=log(yplus)/rklaw+claw
        uplusp=uplus+1./rklaw
        uf1=uf+(abs(utang)-uf*uplus)/uplusp
           resuf=abs((uf1-uf)/uf)
           uf=uf1
           if(resuf.gt.eps1.and.nit.lt.nitmax) then
             yplus=ro*delta*abs(uf)/rey
             nit=nit+1
             go to 500
           endif
           if(nit.eq.nitmax) then
           print*,'nit petit,is=',is
           stop
           endif
           yplus=ro*delta*uf/rey
c
        endif
c
        if(iturb.ne.0) then
        ce(5,is)=0.0
        ce(6,is)=0.0
        un(5,is)=ro*uf*uf/0.3
        un(6,is)=ro*abs(uf*uf*uf)/(rklaw*delta)
        reyturb(is)=cmu*un(5,is)**2/un(6,is)
        endif
cc        print *,yplus,ua(1,is)*uf**2
c
        utgt=1.0
        if(nbp(isc).gt.0) then
        utgt=0.0
        do k=1,nbp(isc)
        isb=ipb(k,isc)
        utgt=utgt+(xt1*ua(2,isb)+xt2*ua(3,isb))/ua(1,isb)
        enddo
        utgt=utgt/float(nbp(isc))
        utgt=sign(1.,utgt)
        endif
      conpa=utgt*ro*uf*uf
      ce(2,is)=ce(2,is)-(xt1*conpa-xn1*pres(is))*rnorm
      ce(3,is)=ce(3,is)-(xt2*conpa-xn2*pres(is))*rnorm
      ce(4,is)=ce(4,is)-abs(conpa*utang*rnorm)
1000  continue
      return
      end
      subroutine mailla
c     -----------------------------------------------------------------
      include 'nsc2ke.inc'
c     -----------------------------------------------------------------
c     local variables definition
      integer is ,nunit
      integer nslogi, nslog2, nslog3, nslog4, nslog5, nslogt
      real    xmin, xmax
      real    ymin, ymax
      nunit   = 14
      open(nunit,file="mesh.inp",status="old")
      read(nunit, * ) ns,nt
      print*,'number of vertices       : ', ns
      print*,'number of triangles      : ', nt
       do is=1,ns
       read(nunit,*) i,coor(1,is),coor(2,is),logfr(is)
       enddo
       do jt=1,nt
       read(nunit,*) k,nu(1,jt),nu(2,jt),nu(3,jt),i
       enddo
c
      close(nunit)
c
      xmin         =  1.0e+3
      xmax         =-1.0e+3
      ymin         =  1.0e+3
      ymax         =-1.0e+3
      do 702 is=1,ns
         xmin      = min(xmin , coor(1,is))
         xmax      = max(xmax , coor(1,is))
         ymin      = min(ymin , coor(2,is))
         ymax      = max(ymax , coor(2,is))
702   continue
c
      print*,'----------------------------------------------'
      print*,'coordinates extrema   :'
      print*,'xmin = ', xmin
      print*,'xmax = ', xmax
      print*,'ymin = ', ymin
      print*,'ymax = ', ymax
c
c     vertices renumerotation strategy
c
      nslogi       = 0
      nslog2       = 0
      nslog3       = 0
      nslog4       = 0
      nslog5       = 0
      nslog6       = 0
      do 30 is=1,ns
         if (logfr(is) .eq. 0) nslogi  = nslogi+1
c
       if(iaxi.ne.0) then
       coor(2,is)=abs(coor(2,is))
       if(logfr(is).ge.4.and.coor(2,is).lt.1.e-6) logfr(is)=6
       endif
c
         if(logfr(is).eq.2.or.
     1     (logfr(is).eq.3.and.ivis.eq.0)) then
         nslog2  = nslog2+1
         log2(nslog2)=is
         endif
c
         numc(is)=0
         if (logfr(is).eq.3) then
         nslog3 = nslog3+1
         log3(nslog3)=is
         numc(is)=nslog3
         endif
c
         if (logfr(is) .eq. 4) then
         nslog4  = nslog4+1
         log4(nslog4)=is
         endif
         if (logfr(is) .eq. 5) then
         nslog5  = nslog5+1
         log5(nslog5)=is
         endif
         if (logfr(is) .eq. 6) then
         nslog6 =nslog6+1
         log6(nslog6)=is
         endif
30    continue
      nslogt = nslogi+nslog2+nslog3+nslog4+nslog5+nslog6
      print*,'----------------------------------------------'
      print*,'internal mesh vertices   : ',nslogi
      print*,'slipping vertices        : ',nslog2
      print*,'no-slipping vertices     : ',nslog3
      print*,'upstream   vertices      : ',nslog5
      print*,'downstream vertices      : ',nslog4
      print*,'inflow profile vertices  : ',nslog6
      print*,'total for verification   : ',nslogt
c
56    format(2e12.5)
57    format(4i6)
78    format(12i6)
c
      return
      end
      program nsc2ke
c     -----------------------------------------------------------------
c
c     2D and AXI euler and navier-stokes equations solver
c
c     explicit multi-steps time integration process
c     upwind schemes and linear interpolation method for
c     the computation of the convective fluxes using a finite
c     volume formulation.
c     classical central galerkin p1-finite element method
c     for the computation of the diffusive fluxes
c
c     k-epsilon turbulence model with two-layer approach or wall laws
c
c     Bijan Mohammadi-Stephane Lanteri
c     INRIA Domaine de Voluceau 78153 Le Chesnay, France
c
c     release 1.0 Jan. 1994
c     Copyright(C) 1994 Bijan Mohammadi-Stephane Lanteri
c
c     Please send bugs and comments to bijan.mohamadi@inria.fr
c
c     -----------------------------------------------------------------
      include 'nsc2ke.inc'
c     -----------------------------------------------------------------
c     reading the application dependent parameters
c
      call config
c
c     loading the global mesh data structure
c
      call mailla
c
      call geowall
c
c     computation of the control volume and triangle areas
c
      call aires
c
c     computation of the control volume altitudes
c
      call clhaut
c
c     construction of the mesh edges
c
      call seg2d
c
      call axigeo
c
c     initializing the physical solution
c
      call init(kt0,t0)
c
      kt                     = kt0
      t                      = t0
      eps                    = 1.0e-10
c
c     swapping the initial physical solution
c
      do 2000 is=1,ns
      do 2000 ivar=1,nvar
         un(ivar,is) = ua(ivar,is)
2000  continue
c
c     computation of the initial pressures
c
      call calprc
c
c     saving the initial physical solution
c
      call result(dt)
c
      open(30,file="RESIDUAL",status="unknown")
c
c     physical time loop
c
100   kt                     = kt + 1
      call caldtl(dtmin,dt)
c
c     runge-kutta time integration process
c
      do 51 ialpha=1,irk
          call resexp(ialpha)
c
c        computation of the pressures
c
         call calprc
c
c        computation of the residual
c
         if (ialpha .eq. 1) then
               som           = 0.0
            do 115 is=1,ns
               som           = som + ce(1,is)*ce(1,is)
     1                             + ce(2,is)*ce(2,is)
     1                             + ce(3,is)*ce(3,is)
     1                             + ce(4,is)*ce(4,is)
115         continue
            som             = sqrt(som)/(float(ns))
            if ((kt - kt0) .eq. 1) som0 = som
            write(30,*) kt,som/som0
c
         endif
c
51    continue
c
c     swapping of the new physical solution
c
      do 200 is=1,ns
      do 200 ivar=1,nvar
         un(ivar,is)  = ua(ivar,is)
200   continue
c
         if(mod(kt,ifre).eq.0) then
         print *,'  kt  =  ',kt,'  residu  =  ',som/som0
         call result(dt)
         endif
c
c     testing for the end of the physical time loop
c     saving the final physical solution
c
      if (kt  .eq. ktmax) then
         print*,'end of execution : maximal number of time steps : ',kt
         goto 300
      endif
      if (som .lt. resf) then
         print*,'end of execution : minimal residual : ',som
         goto 300
      endif
      if (abs(t - tmax) .lt. eps) then
         print*,'end of execution : maximal physical time: ',t
         goto 300
      endif
      goto 100
300   continue
      close(30,status='delete')
         if(mod(kt,ifre).ne.0) then
         print *,'  kt  =  ',kt,'  residu  =  ',som
         call result(dt)
         endif
      stop
      end
       subroutine resexp(ialpha)
c
c explicit resolution
c
       include 'nsc2ke.inc'
c
         do 200 is=1,ns
         do 200 ivar=1,nvar
            ce(ivar,is) = 0.0
200      continue
c
         do is=1,ns
            ua(2,is)         = ua(2,is)/ua(1,is)
            ua(3,is)         = ua(3,is)/ua(1,is)
            ua(5,is)         = ua(5,is)/ua(1,is)
            ua(6,is)         = ua(6,is)/ua(1,is)
         enddo
         call viscdt
c
         if(iflux.eq.1) call fluroe
         if(iflux.eq.2) call fluosh
         if(iflux.eq.3) call flucin
c
         if(iturb.ne.0) call keps2d
c
         do 777 is=1,ns
            ua(2,is)         = ua(2,is)*ua(1,is)
            ua(3,is)         = ua(3,is)*ua(1,is)
            ua(5,is)         = ua(5,is)*ua(1,is)
            ua(6,is)         = ua(6,is)*ua(1,is)
777      continue
c
         if(abs(froud).gt.1.e-3) call gravity
         if(ilaw.ge.1) call vitfrot
         if(iaxi.eq.1) call sa
c
c        boundary conditions treatment
c
         if(iflux.ne.0) call cdl
c
c        updating the physical solution
c
         do 52 is =1,ns
            usais            = dtl(is)/airsa(is)
         do 52 ivar=1,nvar
            ua(ivar,is)  = un(ivar,is) + alpha(ialpha)*usais*ce(ivar,is)
52       continue
c
         do is=1,ns
            ua(5,is)  = abs(ua(5,is))
            ua(6,is)  = abs(ua(6,is))
         enddo
c
        return
       end
      subroutine result(dt)
c     -----------------------------------------------------------------
c     saving on file the result of a computation
c     -----------------------------------------------------------------
      include 'nsc2ke.inc'
c     -----------------------------------------------------------------
      DIMENSION Ft(3),COORt(2,3),VAL(100)
c
       isave=2
      if(isave.eq.0) goto 1234
c
c     saving on file SOL_NS
c
      open(10,file="SOL_NS",status="unknown")
      do is=1,ns
      write(10,*)(ua(i,is),i=1,4)
      enddo
      write(10,*) kt,t
      close(10,status='delete')
c
c     saving on file SOL_KE for k-epsilon
c
      if(iturb.ne.0) then
      open(11,file="SOL_KE",status="unknown")
      do is=1,ns
      reytot=reylam(is)+reyturb(is)
      write(11,*) ua(5,is),ua(6,is),reytot,reyturb(is)
      enddo
      close(11,status='delete')
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      open(31,file='WALL.DATA',status='unknown')
      xcoef=1./sqrt(reynolds)
      if(ivis.eq.0) nnn = nslog2
      if(ivis.eq.1) nnn = nslog3
      do i=1,nnn
      if(ivis.eq.0) is = log2(i)
      if(ivis.eq.1) is = log3(i)
      ro=ua(1,is)
      reytot=reylam(is)+reyturb(is)
      ux=dx(2,is)
      uy=dy(2,is)
      vx=dx(3,is)
      vy=dy(3,is)
      tx=(dx(4,is)/ro-pres(is)*dx(1,is)/ro**2)/gam1
      ty=(dy(4,is)/ro-pres(is)*dy(1,is)/ro**2)/gam1
      xn=vnox(is)
      yn=vnoy(is)
      cp=-2.0*(pres(is)-pin)
c      cp=pres(is)/pin
      cf=2.0*((ux*xn+uy*yn)*yn-xn*(vx*xn+vy*yn))*reytot
      ch=2.0*gam*((reylam(is)/pr)+(reyturb(is)/prt))
     1            *(tx*xn+ty*yn)
      write(31,*) coor(1,is),cp,cf*xcoef,ch*xcoef
      enddo
      close(31,status='delete')
c
1234  continue
      print*,'----------------------------------------------'
      print 10,t,dt
10      format( '    physical time     = ',f10.4,
     &        /,'    minimal time step = ',e15.8)
c
      nss          = 0
      xmasu        = 0.0
      xmanu        = 1.0e+8
      pmin         = 1.0e+8
      pmax         = 0.0
      rmin         = 1.0e+8
      rmax         = 0.0
      umin         = 1.0e+8
      umax         =-1.0e+8
      vmin         = 1.0e+8
      vmax         =-1.0e+8
      xkmin        = 1.0e+8
      xkmax        =-1.0e+8
      xemin        = 1.0e+8
      xemax        =-1.0e+8
      tempmi        = 1.0e+8
      tempma        =-1.0e+8
      do 100 is=1,ns
         uz        = ua(2,is)/ua(1,is)
         vz        = ua(3,is)/ua(1,is)
         qz        = sqrt(uz*uz + vz*vz)
         pp        = pres(is)
         temp      = pp/(gam1*ua(1,is))
         cc        = sqrt(gam*pp/ua(1,is))
         rmach     = qz/cc
         xmasu     = max(xmasu, rmach)
         xmanu     = min(xmanu, rmach)
         umin      = min(umin , uz)
         umax      = max(umax , uz)
         vmin      = min(vmin , vz)
         vmax      = max(vmax , vz)
         pmin      = min(pmin , pres(is))
         pmax      = max(pmax , pres(is))
         rmin      = min(rmin , ua(1,is))
         rmax      = max(rmax , ua(1,is))
         tempmi      = min(tempmi , temp)
         tempma      = max(tempma , temp)
         xkmin      = min(xkmin , ua(5,is))
         xkmax      = max(xkmax , ua(5,is))
         xemin      = min(xemin , ua(6,is))
         xemax      = max(xemax , ua(6,is))
         if (rmach .gt. 1.) nss = nss+1
100   continue
c
      print 998,nss
      print 1000,xmasu,xmanu,pmax,pmin,rmax,rmin,umax,umin,vmax,vmin
      if(iturb.ne.0) print 1001,xkmax,xkmin,xemax,xemin
998   format('    number of supersonic points = ',i4)
1000  format(/
     &     , '    maximum mach number            = ', e15.8,/
     &     , '    minimum mach number            = ', e15.8,/
     &     , '    maximum pressure               = ', e15.8,/
     &     , '    minimum pressure               = ', e15.8,/
     &     , '    maximum density                = ', e15.8,/
     &     , '    minimum density                = ', e15.8,/
     &     , '    maximum x-velocity             = ', e15.8,/
     &     , '    minimum x-velocity             = ', e15.8,/
     &     , '    maximum y-velocity             = ', e15.8,/
     &     , '    minimum y-velocity             = ', e15.8,/)
1001   format(/
     &     , '    maximum k                      = ', e15.8,/
     &     , '    minimum k                      = ', e15.8,/
     &     , '    maximum epsilon                = ', e15.8,/
     &     , '    minimum epsilon                = ', e15.8,/)
c
        if(isave.ne.2) return
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(.false.)then
        open(21,file='GNU.MESH',status='unknown')
        open(22,file='GNU.PRES',status='unknown')
        open(23,file='GNU.MACH',status='unknown')
        open(24,file='GNU.TURB',status='unknown')
        open(25,file='GNU.TEMP',status='unknown')
        open(26,file='GNU.VECT',status='unknown')
cccc
c mesh for gnuplot
cccc
        do iseg =1,nseg
          nubo1=nubo(1,iseg)
          nubo2=nubo(2,iseg)
          write(21,*) coor(1,nubo1),coor(2,nubo1)
          write(21,*) coor(1,nubo2),coor(2,nubo2)
          write(21,*)
          if((logfr(nubo1).ne.0).and.logfr(nubo2).ne.0) then
            write(22,*) coor(1,nubo1),coor(2,nubo1)
            write(22,*) coor(1,nubo2),coor(2,nubo2)
            write(22,*)
            write(23,*) coor(1,nubo1),coor(2,nubo1)
            write(23,*) coor(1,nubo2),coor(2,nubo2)
            write(23,*)
            write(24,*) coor(1,nubo1),coor(2,nubo1)
            write(24,*) coor(1,nubo2),coor(2,nubo2)
            write(24,*)
            write(25,*) coor(1,nubo1),coor(2,nubo1)
            write(25,*) coor(1,nubo2),coor(2,nubo2)
            write(25,*)
            write(26,*) coor(1,nubo1),coor(2,nubo1)
            write(26,*) coor(1,nubo2),coor(2,nubo2)
            write(26,*)
          endif
        enddo
        niso=30
ccc
c iso-pressure for gnuplot
ccc
        ifile=22
        delro=(pmax-pmin)/niso
        do ii=1,niso+1
          val(ii)=pmin+(ii-1)*delro
        enddo
        do it=1,nt
          do i=1,3
            coort(1,i)=coor(1,nu(i,it))
            coort(2,i)=coor(2,nu(i,it))
            ft(i)=pres(nu(i,it))
          enddo
          call isovat(ifile,ft,coort,niso,val,xx0,yy0,xx1,yy1)
        enddo
ccc
c iso-mach for gnuplot
ccc
        ifile=23
        delro=(xmasu-xmanu)/niso
        do ii=1,niso+1
          val(ii)=xmanu+(ii-1)*delro
        enddo
        do it=1,nt
          do i=1,3
            is=nu(i,it)
            coort(1,i)=coor(1,is)
            coort(2,i)=coor(2,is)
            uz        = ua(2,is)/ua(1,is)
            vz        = ua(3,is)/ua(1,is)
            qz        = sqrt(uz*uz + vz*vz)
            pp        = pres(is)
            cc        = sqrt(gam*pp/ua(1,is))
            rmach     = qz/cc
            ft(i)=rmach
          enddo
          call isovat(ifile,ft,coort,niso,val,xx0,yy0,xx1,yy1)
        enddo
ccc
c iso-temp for gnuplot
ccc
        ifile=25
        delro=(tempma-tempmi)/niso
        do ii=1,niso+1
          val(ii)=tempmi+(ii-1)*delro
        enddo
        do it=1,nt
          do i=1,3
            is=nu(i,it)
            coort(1,i)=coor(1,is)
            coort(2,i)=coor(2,is)
            temp = pres(is)/(gam1*ua(1,is))
            ft(i)=temp
          enddo
          call isovat(ifile,ft,coort,niso,val,xx0,yy0,xx1,yy1)
        enddo
c
cc velocity vector for gnuplot
c
        ifile=26
        deltat=1.e-2
        do is=1,ns
          xx0=coor(1,is)
          yy0=coor(2,is)
          xx1=xx0+ua(2,is)*deltat/ua(1,is)
          yy1=yy0+ua(3,is)*deltat/ua(1,is)
          write(26,*) xx0,yy0
          write(26,*) xx1,yy1
          write(26,*)
        enddo
ccc
c iso-k for gnuplot
ccc
        if(iturb.ne.0) then
          ifile=24
          delro=(xkmax-xkmin)/niso
          do ii=1,niso+1
            val(ii)=xkmin+(ii-1)*delro
          enddo
          do it=1,nt
            do i=1,3
              coort(1,i)=coor(1,nu(i,it))
              coort(2,i)=coor(2,nu(i,it))
              ft(i)=ua(5,nu(i,it))/ua(1,nu(i,it))
            enddo
            call isovat(ifile,ft,coort,niso,val,xx0,yy0,xx1,yy1)
          enddo
        endif
c
        close(21)
        close(22)
        close(23)
        close(24)
        close(25)
        close(26)
      endif

      open(1,file='nsc2ke.2dv',status='unknown')
      open(2,file='nsc2ke.v2d',status='unknown')
      open(3,file='nsc2ke.plt',status='unknown')
      write(3,3)
    3 format('TITLE="nsc2ke by Bijan Mohammadi"')
      write(3,4)
    4 format('VARIABLES="X" "Y" "P" "U" "V"')
      write(1,'(i5)')ns
      write(3,5)ns,nt
    5 format('ZONE T="pressure" N=',I5,', E=',I5,
     &', F=FEPOINT, ET=TRIANGLE')
      write(3,6)
    6 format('DT=(SINGLE SINGLE SINGLE SINGLE SINGLE)')
      do is=1,ns
        write(1,7)coor(1,is),coor(2,is),pres(is)
        write(3,7)coor(1,is),coor(2,is),pres(is),0D0,0D0
    7   format(5E13.5)
      enddo
      write(1,'(i5)')nt
      do it=1,nt
        write(1,8)(nu(i,it),i=1,3)
        write(3,8)(nu(i,it),i=1,3)
    8   format(3I5)
      enddo
      write(3,9)
    9 format('ZONE T="velocities", F=POINT')
      do it=1,nt
        xx=0.
        yy=0.
        uz=0.
        vz=0.
        do i=1,3
          is=nu(i,it)
          xx=xx+coor(1,is)
          yy=yy+coor(2,is)
          uz=uz+ua(2,is)/ua(1,is)
          vz=vz+ua(3,is)/ua(1,is)
        enddo
        xx=xx/3.
        yy=yy/3.
        uz=uz/3.
        vz=vz/3.
        write(2,7)xx,yy,uz,vz
        write(3,7)xx,yy,0.,uz,vz
      enddo
      close(3)
      close(2)
      close(1)

      return
      end
      subroutine sa
c
cc axisymmetric corrections
c
      include 'nsc2ke.inc'
c
       dtier=2./3.
c
      do is=1,ns
c
      airr=airsa(is)
      air=airs(is)
      r=rrs(is)
      reyt=reylam(is)+reyturb(is)
      p=pres(is)
      u=ua(2,is)/ua(1,is)
      v=ua(3,is)/ua(1,is)
      ux=dx(2,is)
      vx=dx(3,is)
      uy=dy(2,is)
      vy=dy(3,is)
      divu=ux+vy+v/r
c
      ce(2,is)=ce(2,is)-air*ivis*dtier*reyt*vx
c
      xcof=p-ivis*dtier*reyt*(vy-v/r)
      tautt=reyt*(2*v/r-dtier*divu)
      ce(3,is)=ce(3,is)+air*(xcof-tautt)
c
      xcof1=u*vx+v*(ux+2*vy)-v**2/r
      ce(4,is)=ce(4,is)-air*ivis*dtier*reyt*xcof1
c
      if(coor(2,is).lt.1.e-6.and.logfr(is).ne.3) then
      ce(3,is)=0.0
      endif
c
      enddo
c
          return
          end
      subroutine seg2d
c     -----------------------------------------------------------------
c     construction of the edges of the triangulation
c     -----------------------------------------------------------------
      include 'nsc2ke.inc'
c     -----------------------------------------------------------------
c     local variables definition
      parameter (nvmax = 15)
      integer is    , jt  , kv    , k
      integer iseg
      integer is1   , is2
      integer nuor  , nuex, nub1
      integer jaret(nn,nvmax,2)   , nor(3) , nex(3)
      real    x1,x2 , x3  , y1, y2, y3
      real    esp   , cap
      real    vx    , vy
      real    dbxx(3)    , dbyy(3)
      save
c
      nex(1)                 = 2
      nex(2)                 = 3
      nex(3)                 = 1
      nor(1)                 = 1
      nor(2)                 = 2
      nor(3)                 = 3
c
      do 2 iseg=1,nnsg
         vnocl(1,iseg)       = 0.0
         vnocl(2,iseg)       = 0.0
         vnocl(3,iseg)       = 0.0
2     continue
      do 10 is=1,ns
         do 5 kv=1,nvmax
            jaret(is,kv,1)   = 0
            jaret(is,kv,2)   = 0
5        continue
10    continue
c
      esp                    = 1.0/6.0
      nseg                   = 0
c
      do 1000 jt=1,nt
         is1                 = nu(1,jt)
         is2                 = nu(2,jt)
         is3                 = nu(3,jt)
         x1                  = coor(1,is1)
         y1                  = coor(2,is1)
         x2                  = coor(1,is2)
         y2                  = coor(2,is2)
         x3                  = coor(1,is3)
         y3                  = coor(2,is3)
         dbxx(1)             = y2 - y3
         dbxx(2)             = y3 - y1
         dbxx(3)             = y1 - y2
         dbyy(1)             = x3 - x2
         dbyy(2)             = x1 - x3
         dbyy(3)             = x2 - x1
c
         do 500 k=1,3
            is1              = nu(nor(k),jt)
            is2              = nu(nex(k),jt)
            do 100 kv=1,nvmax
               if (jaret(is1,kv,1) .eq. 0) goto 110
               if (jaret(is1,kv,1) .eq. is2) then
                  ig0        = jaret(is1,kv,2)
                  nub1       = iabs(nubo(1,ig0))
                  nubo(1,ig0)=-nubo(1,ig0)
                  nub2       = nubo(2,ig0)
                  if (is1 .eq. nub1 .and. is2 .eq. nub2) then
                     esp     = abs(esp)
                  else
                     if (is1 .eq. nub2 .and. is2 .eq. nub1) then
                        esp  =-abs(esp)
                     else
                        print*,'error in seg2d ',is1,is2,nub1,nub2
                     endif
                  endif
                  goto 200
               endif
100         continue
c
110         do 120 kv=1,nvmax
               if (jaret(is2,kv,1) .eq. 0) goto 130
               if (jaret(is2,kv,1) .eq. is1) then
                  ig0        = jaret(is2,kv,2)
                  nub1       = nubo(1,ig0)
                  nubo(1,ig0)=-nubo(1,ig0)
                  nub2       = nubo(2,ig0)
                  if (is1 .eq. nub1 .and. is2 .eq. nub2) then
                     esp     = abs(esp)
                  else
                     if (is1 .eq. nub2 .and. is2 .eq. nub1) then
                        esp  =-abs(esp)
                     else
                        print*,'error in seg2d ',is1,is2,nub1,nub2
                     endif
                  endif
                  goto 200
               endif
120         continue
c
            kv               = kv + 1
            if (kv .gt. nvmax) then
               print*,'increase nvmax : ', nvmax
               stop 'seg2d'
            endif
130         nseg  = nseg + 1
            if (nseg .gt. nnsg) then
               print*,'increase nnsg  : ', nnsg
               stop 'seg2d'
            endif
c
            jaret(is2,kv,1)  = is1
            jaret(is2,kv,2)  = nseg
            nubo(1,nseg)     =-is1
            nubo(2,nseg)     = is2
            ig0              = nseg
            esp              = abs(esp)
200         continue
c
c           computation of the control volume boundary normals
c
            vnocl(1,ig0)     = vnocl(1,ig0) +
     &                         esp*(dbxx(nor(k)) - dbxx(nex(k)))
            vnocl(2,ig0)     = vnocl(2,ig0) +
     &                         esp*(dbyy(nor(k)) - dbyy(nex(k)))
500      continue
1000  continue
c
      do 3000 iseg=1,nseg
         if (abs(vnocl(1,iseg)) .lt. 1.0e-16 .and.
     &       abs(vnocl(2,iseg)) .lt .1.0e-16) then
             vnocl(1,iseg)   = 1.0e-8
             vnocl(2,iseg)   = 1.0e-8
         endif
3000  continue
c
c     normalization of the control volume boundary normals
c     construction of the boundary edges
c
      nfr                    = 0
      do 2000 iseg=1,nseg
         cap                 = sqrt(vnocl(1,iseg)*vnocl(1,iseg) +
     &                              vnocl(2,iseg)*vnocl(2,iseg))
         vnocl(1,iseg)       = vnocl(1,iseg)/cap
         vnocl(2,iseg)       = vnocl(2,iseg)/cap
         vnocl(3,iseg)       = cap
         nub1                = nubo(1,iseg)
         if (nub1 .lt. 0) then
c
c           this is a boundary edge
c
            nfr              = nfr + 1
            nufr(nfr)        = iseg
            nubo(1,iseg)     =-nubo(1,iseg)
         endif
2000  continue
c
      print*,'----------------------------------------------'
      print*,'total number of edges    : ', nseg
      print*,'number of boundary edges : ', nfr
c
c     computation of the boundary vertices normals
c          normales sortantes (tauy, -taux)
c
      do 5000 is=1,ns
         vnox(is)            = 0.0
         vnoy(is)            = 0.0
5000   continue
      do 6000 ifr=1,nfr
         iseg                = nufr(ifr)
         nuor                = nubo(1,iseg)
         nuex                = nubo(2,iseg)
         vx                  = coor(2,nuex) - coor(2,nuor)
         vy                  = coor(1,nuor) - coor(1,nuex)
            vnox(nuor)       = vnox(nuor) + vx * 0.5
            vnoy(nuor)       = vnoy(nuor) + vy * 0.5
            vnox(nuex)       = vnox(nuex) + vx * 0.5
            vnoy(nuex)       = vnoy(nuex) + vy * 0.5
6000  continue
c
      do 5550 is=1,ns
         if(logfr(is).ne.0) then
         vno(is)   = sqrt(vnox(is)**2+vnoy(is)**2)
         vnox(is)  = vnox(is)/vno(is)
         vnoy(is)  = vnoy(is)/vno(is)
                            else
         vno(is)   = 0.0
         vnox(is)  = 0.0
         vnoy(is)  = 0.0
         endif
5550  continue
c
      return
      end
      subroutine viscdt
c     -----------------------------------------------------------------
c     computation of the hermitian nodal gradients
c     computation of the viscous fluxes
c     -----------------------------------------------------------------
      include 'nsc2ke.inc'
c     -----------------------------------------------------------------
c     local variables definition
      real    dbxx(3)  , dbyy(3), dxt(4), dyt(4),
     &        uph(4,3) , um(2)  , ei(3)
c
      us6               = 1.0/6.0
      us3               = 1.0/3.0
      us43              = 1.0/4.0/3.0
      gampr             = gam/pr
      gamprt            = gam/prt
      cmu               = 0.09
      cep               = 1./1.3
c
c     initializing the hermitian nodal gradients
c
      do 2 is=1,ns
         dx(1,is)       = 0.0
         dx(2,is)       = 0.0
         dx(3,is)       = 0.0
         dx(4,is)       = 0.0
         dy(1,is)       = 0.0
         dy(2,is)       = 0.0
         dy(3,is)       = 0.0
         dy(4,is)       = 0.0
2     continue
c
c     loop on global list of triangles
c
      do 1000 jt=1,nt
c
c        scatter operation
c        getting the physical solutions for the three vertices
c        of the current triangle
c
         nubo1          = nu(1,jt)
         nubo2          = nu(2,jt)
         nubo3          = nu(3,jt)
c
         uph(1,1)       = ua(1,nubo1)
         uph(2,1)       = ua(2,nubo1)
         uph(3,1)       = ua(3,nubo1)
         uph(4,1)       = pres(nubo1)
c
         uph(1,2)       = ua(1,nubo2)
         uph(2,2)       = ua(2,nubo2)
         uph(3,2)       = ua(3,nubo2)
         uph(4,2)       = pres(nubo2)
c
         uph(1,3)       = ua(1,nubo3)
         uph(2,3)       = ua(2,nubo3)
         uph(3,3)       = ua(3,nubo3)
         uph(4,3)       = pres(nubo3)
c
c        specific internal energy
c
         ei(1)          = gam4*uph(4,1)/uph(1,1)
         ei(2)          = gam4*uph(4,2)/uph(1,2)
         ei(3)          = gam4*uph(4,3)/uph(1,3)
c
c        computation of the p1-gradients
c
         x1             = coor(1,nubo1)
         y1             = coor(2,nubo1)
         x2             = coor(1,nubo2)
         y2             = coor(2,nubo2)
         x3             = coor(1,nubo3)
         y3             = coor(2,nubo3)
c
         dbxx(1)        = y2 - y3
         dbxx(2)        = y3 - y1
         dbxx(3)        = y1 - y2
         dbyy(1)        = x3 - x2
         dbyy(2)        = x1 - x3
         dbyy(3)        = x2 - x1
c
         do 4 k=1,4
            dxt(k)      = uph(k,1)*dbxx(1) +
     &                    uph(k,2)*dbxx(2) +
     &                    uph(k,3)*dbxx(3)
            dyt(k)      = uph(k,1)*dbyy(1) +
     &                    uph(k,2)*dbyy(2) +
     &                    uph(k,3)*dbyy(3)
4        continue
c
         dx(1,nubo1)    = dx(1,nubo1) + dxt(1)
         dx(1,nubo2)    = dx(1,nubo2) + dxt(1)
         dx(1,nubo3)    = dx(1,nubo3) + dxt(1)
         dx(2,nubo1)    = dx(2,nubo1) + dxt(2)
         dx(2,nubo2)    = dx(2,nubo2) + dxt(2)
         dx(2,nubo3)    = dx(2,nubo3) + dxt(2)
         dx(3,nubo1)    = dx(3,nubo1) + dxt(3)
         dx(3,nubo2)    = dx(3,nubo2) + dxt(3)
         dx(3,nubo3)    = dx(3,nubo3) + dxt(3)
         dx(4,nubo1)    = dx(4,nubo1) + dxt(4)
         dx(4,nubo2)    = dx(4,nubo2) + dxt(4)
         dx(4,nubo3)    = dx(4,nubo3) + dxt(4)
c
         dy(1,nubo1)    = dy(1,nubo1) + dyt(1)
         dy(1,nubo2)    = dy(1,nubo2) + dyt(1)
         dy(1,nubo3)    = dy(1,nubo3) + dyt(1)
         dy(2,nubo1)    = dy(2,nubo1) + dyt(2)
         dy(2,nubo2)    = dy(2,nubo2) + dyt(2)
         dy(2,nubo3)    = dy(2,nubo3) + dyt(2)
         dy(3,nubo1)    = dy(3,nubo1) + dyt(3)
         dy(3,nubo2)    = dy(3,nubo2) + dyt(3)
         dy(3,nubo3)    = dy(3,nubo3) + dyt(3)
         dy(4,nubo1)    = dy(4,nubo1) + dyt(4)
         dy(4,nubo2)    = dy(4,nubo2) + dyt(4)
         dy(4,nubo3)    = dy(4,nubo3) + dyt(4)
c
         if (ivis .eq. 0) goto 1000
c
c        computation of the viscous fluxes
c        mean value of the velocity on the current triangle
c
         um(1)          = us3*(uph(2,1) + uph(2,2) + uph(2,3))
         um(2)          = us3*(uph(3,1) + uph(3,2) + uph(3,3))
c
       xmlam=us3*(reylam(nubo1)+reylam(nubo2)+reylam(nubo3))
       xmtur=us3*(reyturb(nubo1)+reyturb(nubo2)+reyturb(nubo3))
       xmtot=xmlam+xmtur
c
         eix            = dbxx(1)*ei(1) + dbxx(2)*ei(2) +
     &                    dbxx(3)*ei(3)
         eiy            = dbyy(1)*ei(1) + dbyy(2)*ei(2) +
     &                    dbyy(3)*ei(3)
c        deformation tensor components
c
      r3  = (dyt(2) + dxt(3))*xmtot
      r2  = (2.0*us3*(2.0*dxt(2) - dyt(3)))*xmtot
      s3  = (2.0*us3*(2.0*dyt(3) - dxt(2)))*xmtot
      r4  = (um(1)*r2+um(2)*r3)+(gampr*xmlam+gamprt*xmtur)*eix
      s4  = (um(1)*r3+um(2)*s3)+(gampr*xmlam+gamprt*xmtur)*eiy
c
c        gathering of the elementary flux into the global one
c
         aitt           = 0.25*airta(jt)
c
         ce(2,nubo1)    = ce(2,nubo1) -
     &                    aitt*(dbxx(1)*r2 + dbyy(1)*r3)
         ce(3,nubo1)    = ce(3,nubo1) -
     &                    aitt*(dbxx(1)*r3 + dbyy(1)*s3)
         ce(4,nubo1)    = ce(4,nubo1) -
     &                    aitt*(dbxx(1)*r4 + dbyy(1)*s4)
c
         ce(2,nubo2)    = ce(2,nubo2) -
     &                    aitt*(dbxx(2)*r2 + dbyy(2)*r3)
         ce(3,nubo2)    = ce(3,nubo2) -
     &                    aitt*(dbxx(2)*r3 + dbyy(2)*s3)
         ce(4,nubo2)    = ce(4,nubo2) -
     &                    aitt*(dbxx(2)*r4 + dbyy(2)*s4)
c
         ce(2,nubo3)    = ce(2,nubo3) -
     &                    aitt*(dbxx(3)*r2 + dbyy(3)*r3)
         ce(3,nubo3)    = ce(3,nubo3) -
     &                    aitt*(dbxx(3)*r3 + dbyy(3)*s3)
         ce(4,nubo3)    = ce(4,nubo3) -
     &                    aitt*(dbxx(3)*r4 + dbyy(3)*s4)
c
1000  continue
c
c     completing the computation of the nodal gradients
c
      do 110 is=1,ns
         ais            = us6/airs(is)
         dx(1,is)       = dx(1,is)*ais
         dx(2,is)       = dx(2,is)*ais
         dx(3,is)       = dx(3,is)*ais
         dx(4,is)       = dx(4,is)*ais
         dy(1,is)       = dy(1,is)*ais
         dy(2,is)       = dy(2,is)*ais
         dy(3,is)       = dy(3,is)*ais
         dy(4,is)       = dy(4,is)*ais
110   continue
c
      return
      end
      subroutine vitfrot
      include 'nsc2ke.inc'
c
      if(ilaw.lt.3) call loglaw
c
      return
      end

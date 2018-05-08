
      SUBROUTINE RUNKBBPOLY(EAR, NE, PARAM, IFL, PHOTAR, PHOTER, NS,NTH, 
     &                  NG, NENER, FLUX0, FLUX, E, GI, AI, THETAI,
     &                  FLUXE, EAR0, NEX1, NEX2)

      IMPLICIT NONE

      Integer IFL, NE, NS, NTH, NG, NENER
      Real EAR(0:NE), PARAM(*), PHOTAR(NE), PHOTER(NE)
      Real FLUX0(ns,nth,ng,nener), FLUX(nener), E(nener), GI(ng) 
      Real AI(ns), THETAI(nth), fluxe(ne,2), ear0(nener)
      Integer nex1(ne,2), nex2(ne,2)
      
c Modified by DJKB to interpolate table with cubic method.
c Being converted to arbitrary order polynomial.
c Lower case 'c's indicate comments by DJKB.
c Note that the energy grid is currently still linearly interpolated but is in the process of being upgraded to cubic.


C     PROGRAM "kerrbb.f"
C
C     Written by: Li-Xin Li
C              Harvard-Smithsonian Center for Astrophysics
C              60 Garden St., Cambridge, MA 02138
C                      
C     Version:         October, 2004
C
C     Reference:  Li, Zimmerman, Narayan, & McClintock 2004, astro-ph/0411583
C
C     The program computes the black-body emission spectrum from a thin 
C     Keplerian accretion disk around a Kerr black hole. All relativistic 
C     effects are taken into account, including frame-dragging, Doppler boost, 
C     gravitational redshift, and bending of light by the gravity of the black 
C     hole. In particular, self-irradiation of the disk as a result of light 
C     deflection is included. The inner boundary of the disk is fixed at the 
C     marginally stable orbit. Torque at the inner boundary of the disk 
C     is allowed to be nonzero. However, when this program is applied to data
C     reduction, a zero torque (eta = 0) is recommended since we have found
C     that the effect of a nonzero torque on the spectrum can, to a good 
C     approximation, be absorbed into a zero torque model by adjusting the
C     mass accretion rate and the normalization.
C
C     The program reads in the precalculated data file kerrbb.fits and 
C     uses it to fit spectral data by linear interpolation. The first
C     extension (GRIDVALS) contains 4 rows specifying:
C     - a grid of black hole spin (a/M);
C     - a grid of disk inclination angles (in degrees);
C     - a grid of the torque at the inner boundary of the disk (eta, dimensionless);
C     - a grid of spectrum energy (in keV)
C     The second extension (FLUX) contains the flux density at each grid 
C     point specified by the above four parameters, in photons/keV/cm^2/sec.
C     [There are four columns in this file: FLUX1, when both self-irradiation 
C     and limb-darkening are turned off; FLUX2, when self-irradiation is off 
C     but limb-darkening is on; FLUX3, when self-irradiation is on but 
C     limb-darkening is off; FLUX4, when both self-irradiation and 
C     limb-darkening are on.]
C
C     The program requires the following input arguments
C     - "EAR(0:NE)": array of observed energy (frequency) bins (in units of 
C        keV);
C     - "NE": number of observed energy bins (size of the energy array);
C     - "PARAM(*)": values of the model parameters, containing
C       # "eta"(=PARAM(1)): ratio of the disk power produced by a
C         torque at the disk inner boundary to the disk power arising 
C         from accretion. It must be >= 0 and <=1. When eta = 0, the 
C         solution corresponds to that of a standard Keplerian disk with 
C         zero torque at the inner boundary;
C       # "a"(=PARAM(2)): specific angular momentum of the black hole in 
C        units of the black hole mass M (geometrized units G=c=1). a 
C         should be >= -1 and < 1;
C       # "i"(=PARAM(3)): disk's inclination angle (the angle between the 
C        axis of the disk and the line of sight). It is expressed in 
C         degrees. i=0 is for a "face-on" accretion disk. i should be <=  
C         85 degree;
C       # "Mbh"(=PARAM(4)): the mass of the black hole in units of the 
C         solar mass;
C       # "Mdd"(=PARAM(5)): the "effective" mass accretion rate of the disk 
C        in units of 10^18 g/sec. When eta = 0 (zero torque at the inner 
C        boundary), this is just the mass accretion rate of the disk. When 
C         eta is nonzero, the effective mass accretion rate = (1+eta) times 
C         the true mass accretion rate of the disk. The total disk 
C         luminosity is then "epsilon" times "the effective mass accretion 
C         rate" times "c^2", where epsilon is the radiation efficiency of a 
C         standard accretion disk around the Kerr black hole;
C       # "Dbh"(=PARAM(6)): the distance from the observer to the black 
C         hole in units of kpc;
C       # "h"(=PARAM(7)): spectral hardening factor, T_col/T_eff. It should
C         be greater than 1.0, and considered to be 1.5-1.9 for accretion 
C         disks around a stellar-mass black hole. See, e.g., Shimura and 
C         Takahara 1995, ApJ, 445, 780;
C       # "rflag"(=PARAM(8)): a flag to switch on/off the effect of 
C         self-irradiation (never allowed to be free). Self-irradiation is 
C         included when rflag is > 0. Self-irradiation is not included when 
C         rflag is <= 0;
C       # "lflag"(=PARAM(9)): a flag to switch on/off the effect of limb-
C         darkening (never allowed to be free). The disk emission is assumed 
C         to be limb-darkened when lflag is > 0. The disk emission is assumed 
C         to be isotropic when lflag is <= 0.
C     [Formally the program also requires an input parameter "spectrum 
C     number" (IFL)--an integer that specifies which data set the energies 
C     are for, as requested by XSPEC (http://heasarc.gsfc.nasa.gov/docs/
C     xanadu/xspec/manual/node60.html for xspec11; http://heasarc.gsfc.nasa.
C     gov/docs/xanadu/xspec/xspec12, Appendix C in the manual for xspec12). 
C     However, this parameter is irrelevant here.]
C
C     The program outputs an array of the observed photon flux in each bin 
C     of energy: "PHOTAR(NE)", in units of photons/cm^2/second. For
C     example, PHOTAR(i) gives the observed photon number flux with 
C     observed photon energy between EAR(i) and EAR(i+1). Formally, flux 
C     errors are also outputted as PHOTER(NE), although they are irrelevant 
C     here.

      Integer IREAD
      Real eta, a, theta, Mbh, Mdd, Dbh, h, lflag, rflag, lflagsav, 
     $     rflagsav
      Real difx, dify, difz, diftx, difty, diftz, dx, dy, dz, 
     $     dife, difte, lslope
      Real ca, cb, cc, t, s, w
      Real de, earm
      Real flux1, flux2, flux3, flux4, flux5, flux6, flux7, flux8
      Real cflux1, cflux2, cflux3, cflux4, cflux5, cflux6, cflux7, 
     $     cflux8, iflux1, iflux2, iflux12, iflux0
      Integer order, or2f, or2c
      Real xfrac(0:int(PARAM(10))), yfrac(0:int(PARAM(10)))
      Real zfrac(0:int(PARAM(10)))
      Real frac(0:int(PARAM(10)),0:int(PARAM(10)),0:int(PARAM(10)))
      Real enfrac(0:int(PARAM(10)))
      Integer nx, ny, nz, nxb, nyb, nzb, nx1, nx2, ny1, ny2, nz1, nz2
      Integer nex, nexb
      Integer i, j, k, l, n1, n2, n3, n4
      Integer ilun, block, status, hdutyp, icol

      Logical qanyf

      character(256) fgmodf, datdir, contxt
      character(256) filenm
      integer lenact
      external lenact, fgmodf

      save iread, lflagsav, rflagsav
      data iread/0/

c suppress a warning message from the compiler
      i = ifl

c this model does not calculate errors
      DO i = 1, ne
         photer(i) = 0.0
      ENDDO
c      print *,'Im in runkbbdb.f'


      eta = PARAM(1)
      a = PARAM(2)
      theta = PARAM(3)
      Mbh = PARAM(4)
      Mdd = PARAM(5)
      Dbh = PARAM(6)
      h = PARAM(7)
      rflag = PARAM(8)
      lflag = PARAM(9)
      order = PARAM(10)

      ca = Mdd**0.25*Mbh**(-0.5)*h
      cb = Mbh*Mbh/(Dbh*Dbh)*h**(-4.)
      cc = cb*ca**2.
      
      if ((iread .ne. 99) .or. (lflagsav .ne. lflag) .or. 
     $       (rflagsav .ne. rflag)) then
         datdir = fgmodf()
         filenm = datdir(:lenact(datdir))//'kerrbb.fits'

         status = 0

C Open FITS input file

         CALL getlun(ilun)
         CALL ftopen(ilun, filenm, 0, block, status)
         contxt = 'Failed to open '//filenm(:lenact(filenm))
         IF ( status .NE. 0 ) GOTO 999

C Move to first (GRIDVALS) extension

         CALL ftmrhd(ilun, 1, hdutyp, status)
         contxt = 'Failed to move to first extension of '//
     &            filenm(:lenact(filenm))
         IF ( status .NE. 0 ) GOTO 999

C Read four sets of grid values

         CALL ftgcve(ilun, 1, 1, 1, ng, 0.0, gi, qanyf, status)
         CALL ftgcve(ilun, 1, 2, 1, ns, 0.0, ai, qanyf, status)
         CALL ftgcve(ilun, 1, 3, 1, nth, 0.0, thetai, qanyf, status)
         CALL ftgcve(ilun, 1, 4, 1, nener, 0.0, e, qanyf, status)
         contxt = 'Failed to read GRIDVALS data from '
     &            //filenm(:lenact(filenm))
         IF ( status .NE. 0 ) GOTO 999

C Move to second (FLUX) extension
          
         CALL ftmrhd(ilun, 1, hdutyp, status)
         contxt = 'Failed to move to second extension of '//
     &            filenm(:lenact(filenm))
         IF ( status .NE. 0 ) GOTO 999

C Read the appropriate flux data
         
         if ((rflag .le. 0.0).and.(lflag .le. 0.0)) then
            icol = 1
         else  if ((rflag .le. 0.0).and.(lflag .gt. 0.0)) then
            icol = 2
         else  if ((rflag .gt. 0.0).and.(lflag .le. 0.0)) then
            icol = 3
         else
            icol = 4
         endif

         CALL ftgcve(ilun, icol, 1, 1, ng*ns*nth*nener, 0.0, 
     &               flux0, qanyf, status)
         contxt = 'Failed to read FLUX data from '//
     &            filenm(:lenact(filenm))
         IF ( status .NE. 0 ) GOTO 999

       CALL ftclos(ilun, status)
         CALL frelun(ilun)
                
       iread = 99
         lflagsav = lflag
         rflagsav = rflag

 999     CONTINUE
         IF ( status .NE. 0 ) THEN
            CALL xwrite(contxt, 10)
            WRITE(contxt, '(a,i6)') 'RUNKBB: Status = ', status
            CALL xwrite(contxt, 10)
         ENDIF
      
      endif
      
      nx2=ns
      do i = 1,ns
       dx = a-ai(i)
       if (dx .lt. 0) then
          nx2 = i
          GOTO 10
       end if
      enddo

  10  ny2=nth
      do i = 1,nth
       dy = theta-thetai(i)
       if (dy .lt. 0) then
          ny2 = i
          GOTO 20
       end if
      enddo

  20  nz2=ng
      do i = 1,ng
       dz = eta-gi(i)
       if (dz .lt. 0) then
          nz2 = i
          GOTO 30
       end if
      enddo

  30  nx1 = nx2-1
      ny1 = ny2-1
      nz1 = nz2-1

c now nx1 is grid point below value
c now nx2 is grid point above value
c run from nx2-2 to nx2+1 unless caught at ends

      or2f=floor(real(order)/2)
      or2c=ceiling(real(order)/2)
      
c catch end cases
      if (nx2 .lt. (or2c+1)) then
        nx2 = or2c+1
      end if
      if (nx2 .gt. ns-or2f) then
        nx2 = ns-or2f
      end if

      if (ny2 .lt. (or2c+1)) then
        ny2 = or2c+1
      end if
      if (ny2 .gt. nth-or2f ) then
        ny2 = nth-or2f
      end if

      if (nz2 .lt. (or2c+1)) then
        nz2 = or2c+1
      end if
      if (nz2 .gt. ng-or2f ) then
        nz2 = ng-or2f
      end if

c calculate fraction from each point
      call polyfrac(a,ai(nx2-or2c:nx2+or2f),order
     $                       ,xfrac)
      call polyfrac(theta,thetai(ny2-or2c:ny2+or2f
     $                ),order,yfrac)
      call polyfrac(eta,gi(nz2-or2c:nz2+or2f)
     $                 ,order,zfrac)
      
c      print *,gi
c      print *,a,':',ai(nx2-or2c:nx2+or2f)
c      print *,theta,':',thetai(ny2-or2c:ny2+or2f)
c      print *,eta,':',gi(nz2-or2c:nz2+or2f)
c      print *,nx2
c      print *,ny2
c      print *,nz2
      
c      print *,xfrac
c      print *,yfrac
c      print *,zfrac
c      print *,zfrac(0)+zfrac(1)+zfrac(2)+zfrac(3)
      
      do l=0,order
        do j=0,order
          do k=0,order
      frac(l,j,k)=xfrac(l)*yfrac(j)*zfrac(k)
          enddo
        enddo
      enddo

c fractions for linear interpolation
c faster to change cfluxes
c      t = 1 + (dx)/(ai(nx2)-ai(nx1))
c      s = 1 + (dy)/(thetai(ny2)-thetai(ny1))
c      w = 1 + (dz)/(gi(nz2)-gi(nz1))


c for each energy
      do i = 1,nener
       ear0(i) = e(i)*ca

c do sum of fluxes from each point

      flux(i)=0
      do l=0,order
        do j=0,order
          do k=0,order
             flux(i)=flux(i)+frac(l,j,k)*flux0(nx2-or2c+l,
     $       ny2-or2c+j,nz2-or2c+k,i)
          enddo
        enddo
      enddo
      flux(i)=flux(i)*cc

c      if (isnan(flux(i))) then
c        print *,'flux(',i,') is nan'
c      endif

      if (flux(i) .lt. 0) then
        flux(i) = 1e-20
c        print *,'flux(',i,') is negative'
      endif

c       flux1 = flux0(nx1,ny1,nz1,i)
c       flux2 = flux0(nx2,ny1,nz1,i)
c       flux3 = flux0(nx2,ny2,nz1,i)
c       flux4 = flux0(nx1,ny2,nz1,i)
c       flux5 = flux0(nx1,ny1,nz2,i)
c       flux6 = flux0(nx2,ny1,nz2,i)
c       flux7 = flux0(nx2,ny2,nz2,i)
c       flux8 = flux0(nx1,ny2,nz2,i)

c       cflux1 = (1.-t)*(1.-s)*(1.-w)*flux1
c       cflux2 = t*(1.-s)*(1.-w)*flux2
c       cflux3 = t*s*(1.-w)*flux3
c       cflux4 = (1.-t)*s*(1.-w)*flux4
c       cflux5 = (1.-t)*(1.-s)*w*flux5
c       cflux6 = t*(1.-s)*w*flux6
c       cflux7 = t*s*w*flux7
c       cflux8 = (1.-t)*s*w*flux8

c       flux(i) = cc*(cflux1 + cflux2 + cflux3 + cflux4 
c     $           + cflux5 + cflux6 + cflux7 + cflux8)
            
      enddo

c      print *,'kerrbb energy table range: ',ear0(1),'-',ear0(nener)

c     For each energy bin
      do i = 1,NE
c     Look at lower and upper ends of the observed energy bin
         do k=1,2
            earm =ear(i+k-2)

C     When bin energy earm < ear0(1), calculate the specific flux at
C     earm according to law: specific flux proportional to energy^(-2/3)
c     i.e. If the observed energy is below the kerrbb table, extrapolate as a standard disc blackbody

            if (earm .lt. ear0(1)) then
               fluxe(i,k) = flux(1)*(ear0(1)/earm)**(2./3.)
               nex1(i,k) = 1
               nex2(i,k) = 1

C     When bin energy earm > ear0(nener), cutoff the flux density

            else if (earm .gt. ear0(nener)) then
               fluxe(i,k) = 1.e-20
               nex1(i,k) = nener
               nex2(i,k) = nener

C     Use linear interpolation when ear0(1) < earm < ear0(nener)
c     but this is to be changed to cubic interpolation

            else
c     Find the closest table energy
               nex = 1
               dife = abs(earm-ear0(1))
               do j = 1,nener
                  difte = abs(earm-ear0(j))
                  if (difte .lt. dife) then
                     nex = j
                     dife = difte
                  end if
               enddo
c     Find the other end of the tabulated E interval over which to interpolate
               de = earm - ear0(nex)
               nexb = nex + INT(sign(1.,de))
c     ? nex1,2 don't actually need to be stored for each energy?
               nex1(i,k) = min(nex,nexb)
               nex2(i,k) = max(nex,nexb)

c     Linearly interpolate fluxes onto observed energy bin ends (in log space)
               lslope = (log10(flux(nex2(i,k)))-log10(flux(nex1(i,k))))
     $               /(log10(ear0(nex2(i,k)))-log10(ear0(nex1(i,k))))
               fluxe(i,k) = flux(nex1(i,k))*10.**(lslope*(log10(earm)
     $               -log10(ear0(nex1(i,k)))))

c     Now use cubic interpolation BUT NOT YET LOGARITHMICALLY - fine for small dE and slowly changing flux
c               call CUBPOLYFRAC(earm,ear0(nex1(i,k)-1),ear0(nex1(i,k)),
c     $               ear0(nex2(i,k)),ear0(nex2(i,k)+1),enfrac)
c
c               fluxe(i,k) = enfrac(0)*flux(nex1(i,k)-1)+
c     $                      enfrac(1)*flux(nex1(i,k)  )+
c     $                      enfrac(2)*flux(nex2(i,k)  )+
c     $                      enfrac(3)*flux(nex2(i,k)+1)

            endif
c     End of k-loop: have interpolated flux onto each end of the observed energy bin
         enddo

c     Indices of lower and upper ends of observed energy bin in table energy bins
         n1=nex1(i,1)
         n2=nex2(i,1)
         n3=nex1(i,2)
         n4=nex2(i,2)

         iflux1=0.5*(fluxe(i,1)+flux(n2))*(ear0(n2)-ear(i-1))
         iflux2=0.5*(fluxe(i,2)+flux(n3))*(ear(i)-ear0(n3))
         iflux12=0.5*(fluxe(i,1)+fluxe(i,2))*(ear(i)-ear(i-1))

         if (n3 .gt. n2) then
            iflux0=0.0
            do k=n2,n3-1
               iflux0=iflux0+0.5*(flux(k)+flux(k+1))*(ear0(k+1)-ear0(k))
            enddo

            PHOTAR(i) = iflux0 + iflux1 + iflux2

         else if (n3 .eq. n2) then
            PHOTAR(i) = iflux1 + iflux2

         else
            PHOTAR(i) = iflux12

         endif
      
      enddo
      
      RETURN
      END

      SUBROUTINE POLYFRAC(X,X_LIST,ORDER,FRAC)

      IMPLICIT NONE

      Integer order, i, j
      Real x, x_list(0:order), frac(0:order)

      do i=0,order
         frac(i) = 1.
         do j=0,order
            if (j .ne. i) then
               frac(i) = frac(i) * (x-x_list(j))/(x_list(i)-x_list(j))
            endif
         enddo
      enddo
      
      RETURN
      END

C	Synthetic eddy method module for NEK5000
C
C User needs to provide 
C  -dimensions of the eddy box (in usrdat)
c  -A few input parameters 
c     (e.g. # of eddies, max(sigma), ..
C  -sem_input.txt 
c     containing 
c       -1d-pos
c       -k(pos)
c       -epsilon(pos)
c       -umean(pos)
c
c   The latter two are are e.g. 
c   generated with Matlab-Script from DNS-Pipedata
C
C	Lorenz Hufnagel hufnagel@kth.se
C Based on code by Oana Marin
C

      module SEM
         implicit none
         save

         integer             :: nEddy ! number of eddies 

         integer             :: nElInlet ! number of elements at inlet
         ! Warning: This optimization assumes that actually the first
         ! nElInlet mesh-elements are the inlet elements.
         ! In the pipe-case the mesh is created accordingly, 
         ! and nElInlet = nelperface

         real                :: yplus_cutoff 
         ! yplus above which no eddies are generated, to maintain
         ! divergence freedom. 
         ! Typically this should be 10 delta+ from the wall.
         real                :: sigma_max 
         ! upper eddy size bound (due to k&eps)
         real                :: bbox_max
         ! maximum extent of eddies (radial direction here)
         real                :: u0 ! bulk velocity


         ! Deviation of bulk velocity of generated signal
         real                :: bulk_vel_diff ! bulk velocity
         real                :: inlet_area 
         
         character*80 infile
         parameter (infile='sem_input.txt')

         ! extent, Volume of the virtual SEM domain
         real                :: ybmin,ybmax,zbmin,zbmax,xbmin,xbmax 
         real                :: Vb
         real                :: z_inlet ! z value of inlet plane

         real , allocatable  :: umean_inlet(:,:,:)
         real , allocatable  :: sigma_inlet(:,:,:)
         ! eddy size at inlet face
         real , allocatable  :: intensity_inlet(:,:,:) 
         ! cholesky of isotropic diag. reystress tensor
         ! at inlet face

         ! Actual prescribed velocity field 
         ! (relevant elements, are extracted via 'v  ' bc)
         real , allocatable  :: u_sem(:,:,:,:)
         real , allocatable  :: v_sem(:,:,:,:)
         real , allocatable  :: w_sem(:,:,:,:)

         ! Individual eddy energies
         real, allocatable   :: ex(:),ey(:),ez(:),eps(:,:)
         integer, allocatable :: eddy_pt(:)
         
         integer             :: nInputdata
         real   ,allocatable :: pos(:),umean(:),tke(:),dissip(:)
      end module SEM

c-----------------------------------------------------------------------

      subroutine SEMinit()
      use SEM
      use AVG, only: nElperFace
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'  ! L[XYZ]1,LELV
      include 'GEOM_DEF'
      include 'GEOM' ! XM1
      include 'PARALLEL_DEF'
      include 'PARALLEL'

      real vel_interp, tke_interp, dissip_interp, sigmal, radius
      real work
      integer fid, nlines
      integer e, i, j, eg

      logical semstop

      if (nElInlet.ne.nElperFace) then
        if (nid.eq.0) write(*,*) 
     $'ABORT For pipe flow, you probably want nElperFace==nElinlnet'
        call exitti
      endif

c
c     CHECK FOR INPUT FILE
c
      inquire(file=infile,exist=semstop); semstop=.not.semstop
      if (semstop) then

         write(*,*) 'Warning: ',trim(infile),' does not exist:',
     &              ' simulation will stop.'

         stop
      end if
c
c     Read infile
      
      fid = 35
      open(unit=fid,file=infile,form='formatted')
      
      read(fid,*)      ! skip header
      read(fid,*) nlines

      nInputdata = nlines

      allocate(pos(nlines))
      allocate(umean(nlines))
      allocate(tke(nlines))
      allocate(dissip(nlines))

      do i=1, nlines
        read(fid,*) pos(i), umean(i), 
     &              tke(i), dissip(i)
      enddo
      close(fid)


      allocate(u_sem(lx1,ly1,lz1,lelv))
      allocate(v_sem(lx1,ly1,lz1,lelv))
      allocate(w_sem(lx1,ly1,lz1,lelv))


      allocate(ex(neddy))
      allocate(ey(neddy))
      allocate(ez(neddy))
      allocate(eps(3,neddy))
      allocate(eddy_pt(neddy))

      allocate(sigma_inlet(ly1,lz1,nElInlet))
      allocate(intensity_inlet(ly1,lz1,nElInlet))
      allocate(umean_inlet(ly1,lz1,nElInlet))

      inlet_area = 0

      do e=1,nelv
        eg = lglel(e)

        ! Calculate eddy size and intensity at inlet only once
        if (eg.le.nElInlet) then

          do j=1,ly1
          do i=1,lx1
          
          radius = sqrt(xm1(i,j,1,e)**2 + ym1(i,j,1,e)**2)

          call SEMinputData(radius,vel_interp,tke_interp,dissip_interp)

          sigmal = (tke_interp**1.5)/dissip_interp
          sigmal = max(.5*sigmal, 1e-8)  ! avoid numerical instability

          ! Optional: Limit eddy size far away from wall. 
          ! Suggested in Jarrins PhD, but not implemented in Code Saturne
          ! kappa is .41, 0.5 is pipe radius, therefore .41*.5
          ! sigmal = max(.5*min(sigmal,0.41*0.5),1e-8)  

          sigma_inlet(i,j,eg)     = sigmal
          umean_inlet(i,j,eg)     = vel_interp
          if (radius > yplus_cutoff) then
            intensity_inlet(i,j,eg) = 0.0d0
          else
            intensity_inlet(i,j,eg) = sqrt(2./3.*tke_interp)
          endif

!-----Pick face 5 to evaluate surface Jacobian
          inlet_area = inlet_area + area(i,j,5,e)

          enddo
          enddo

        endif
      enddo

      call gop(inlet_area,work,'+  ',1)

      end subroutine SEMinit

c-----------------------------------------------------------------------

!     read parameters SEM
      subroutine SEM_param_in(fid)
        use SEM, only: nEddy, yplus_cutoff, nElInlet, 
     $                 sigma_max, bbox_max, u0
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            !
      include 'PARALLEL_DEF' 
      include 'PARALLEL'        ! ISIZE, WDSIZE, LSIZE,CSIZE

!     argument list
      integer fid               ! file id

!     local variables
      integer ierr

!     namelists
      namelist /SEM_list/ nEddy, nElInlet, yplus_cutoff,
     $        sigma_max, bbox_max, u0
!-----------------------------------------------------------------------
!     default values
      nEddy = 5000
      nElInlet = 256
      yplus_cutoff = 0.4
      sigma_max = 0.25
      bbox_max = 0.05
      u0 = 1.0
!     read the file
      ierr=0
      if (NID.eq.0) then
         read(unit=fid,nml=SEM_list,iostat=ierr)
      endif
      call err_chk(ierr,'Error reading SEM parameters.$')

!     broadcast data
      call bcast(nEddy,ISIZE)
      call bcast(nElInlet,ISIZE)
      call bcast(yplus_cutoff, WDSIZE)
      call bcast(sigma_max, WDSIZE)
      call bcast(bbox_max, WDSIZE)
      call bcast(u0, WDSIZE)

      return
      end  subroutine SEM_param_in
!***********************************************************************
!     write parameters checkpoint
      subroutine SEM_param_out(fid)
        use SEM, only: nEddy, yplus_cutoff, nElInlet, 
     $                 sigma_max, bbox_max, u0
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            !

!     argument list
      integer fid               ! file id

!     local variables
      integer ierr

!     namelists
      namelist /SEM_list/ nEddy, nElInlet, yplus_cutoff,
     $        sigma_max, bbox_max, u0
!-----------------------------------------------------------------------
      ierr=0
      if (NID.eq.0) then
         write(unit=fid,nml=SEM_list,iostat=ierr)
      endif
      call err_chk(ierr,'Error writing SEM parameters.$')

      return
      end subroutine SEM_param_out

C=======================================================================
C SEM Main routine
C=======================================================================
      subroutine synthetic_eddies
        use SEM
      implicit none
      include 'SIZE_DEF'
      include 'SIZE'
      include 'TSTEP_DEF'
      include 'TSTEP' ! ISTEP,IOSTEP
      include 'GEOM_DEF'
      include 'GEOM' ! XM1
      include 'PARALLEL_DEF'
      include 'PARALLEL' ! XM1
    
      real    ff,fx,fy,fz,
     &           rr, rrx,rry,rrz
      real, parameter :: sqrt32 = sqrt(3./2)

      integer neddy_ll
      save    neddy_ll

      real  work
      real wk_e(neddy*3)
      integer clock, e,eg, 
     &      i,j,i_l,ne,nv,iseed
      !     functions
      real dnekclock, sqrtVn
  
c --- generate initial eddy distribution ----
    
      if (istep.eq.0) then
          clock = int(dnekclock())
          iseed=(nid+1)*11 + clock + 1
          call  ZBQLINI(iseed)

          call distribute_eddy(neddy,neddy_ll,eddy_pt)
c     Generate local eddies with locations ex,ey,ez
          do i=1,neddy_ll
              i_l = eddy_pt(i)
              call gen_eddy(ex,ey,ez,eps,i_l)
          enddo   
c     zeroing eddies off this proc
          call zero_nonlocal_eddies(ex,ey,ez,eps,neddy_ll,neddy,eddy_pt)
      else 
          call advect_recycle_eddies(ex,ey,ez,eps,neddy_ll,eddy_pt)
          call zero_nonlocal_eddies(ex,ey,ez,eps,neddy_ll,neddy,eddy_pt)
      endif

c     Gather eddies on each processor
c     possibly replace with glvadd

      call gop(ex,wk_e,'+  ',neddy) ! inplace!
      call gop(ey,wk_e,'+  ',neddy)
      call gop(ez,wk_e,'+  ',neddy)
      call gop(eps,wk_e,'+  ',neddy*3)

c ---- compute velocity contribution of eddies ------

      nv = nx1*ny1*nz1*nelv ! max number of gird points per processor

      call rzero(u_sem,nv)
      call rzero(v_sem,nv)
      call rzero(w_sem,nv)
      bulk_vel_diff = 0

      sqrtVn = sqrt(Vb/real(neddy))

      do e=1,nelv
        eg = lglel(e)
c         if (abs(zm1(1,1,1,e)-z_inlet).lt.1e-13) then
        if (eg.le.nElInlet) then 

        do j=1,ly1
        do i=1,lx1

        ! NOTE that the pipe/mesh is rotated in usrdat2()
        ! hence, the inverted indexing into ym1/zm1

           u_sem(i,j,1,e) = 0
           v_sem(i,j,1,e) = 0
           w_sem(i,j,1,e) = umean_inlet(i,j,eg)

           do ne=1,neddy
            rrx = (xm1(i,j,1,e)-ex(ne))
            rry = (ym1(i,j,1,e)-ey(ne))
            rrz = (zm1(i,j,1,e)-ez(ne))

            rr = sqrt(rrx*rrx + rry*rry + rrz*rrz)

            if (rr.lt.sigma_inlet(i,j,eg)) then
              fx = sqrt32*(1.0-abs(rrx)/sigma_inlet(i,j,eg))
              fy = sqrt32*(1.0-abs(rry)/sigma_inlet(i,j,eg))
              fz = sqrt32*(1.0-abs(rrz)/sigma_inlet(i,j,eg))

              ff=fx*fy*fz*sqrtVn/sqrt(sigma_inlet(i,j,eg))**3

              u_sem(i,j,1,e) = u_sem(i,j,1,e) + 
     $                      intensity_inlet(i,j,eg)*eps(1,ne)*ff
              v_sem(i,j,1,e) = v_sem(i,j,1,e) + 
     $                      intensity_inlet(i,j,eg)*eps(2,ne)*ff
              w_sem(i,j,1,e) = w_sem(i,j,1,e) + 
     $                      intensity_inlet(i,j,eg)*eps(3,ne)*ff
            endif
           enddo         


           bulk_vel_diff = bulk_vel_diff + area(i,j,5,e)*w_sem(i,j,1,e)

        enddo
        enddo
      endif 
      enddo


      call gop(bulk_vel_diff,work,'+  ',1)

      return
      end subroutine synthetic_eddies
c-----------------------------------------------------------------------
c     distribute eddies across all processors
c     neddy - total # eddies
c     nl    - local # eddies
c     pt    - local pointer
      subroutine distribute_eddy(neddy,nl,pt)
      implicit none
      include 'SIZE_DEF'
      include 'SIZE'
      include 'PARALLEL_DEF'
      include 'PARALLEL'

      integer, intent(in) :: neddy
      integer, intent(out) :: nl, pt(neddy)
      integer i, nbatch0,nbatch1,ibnd,istart

      nbatch0 = neddy/np + 1
      nbatch1 = neddy/np
      ibnd = neddy - nbatch1*np - 1
      if(nid.le.ibnd)then
        nl     = nbatch0
        istart = nbatch0*nid + 1
      else
        nl     = nbatch1
        istart = nbatch0*(ibnd+1) + nbatch1*(nid-(ibnd+1))+ 1
      endif

      do i=1,nl
         pt(i) = istart + (i-1)
      enddo

      return
      end subroutine distribute_eddy
c-----------------------------------------------------------------------
c     Generate eddy location randomly in bounding box
      subroutine gen_eddy(ex,ey,ez,eps,n)
      use SEM, only: zbmin, zbmax, ybmax 
      implicit none

      real, parameter :: twoPi = 6.283185307179586476925286766
      double precision, external :: rnd_loc
      include 'SIZE_DEF'
      include 'SIZE'
      include 'TSTEP_DEF'
      include 'TSTEP' ! ISTEP,IOSTEP

      real, intent(inout) :: ex(1),ey(1),ez(1),eps(3,1)
      integer, intent(in) :: n
      real rnd, rho, theta
      integer j,i_l

      i_l = n

c     Generate uniformly distributed random locations 
c     in polar coordinates (!)
c     ybmax contains actual radius for origin-centered pipe
      rho = ybmax*sqrt(rnd_loc(0.0,1.0))  
      theta = rnd_loc(0.,twoPI) 

      ex(i_l) = rho * cos(theta)
      ey(i_l) = rho * sin(theta)
      if(istep.eq.0)then
          ez(i_l) = rnd_loc(zbmin,zbmax)
      else
          ez(i_l) = zbmin
      endif       

      do j=1,3
        rnd = rnd_loc(0.0,1.0)
        if (rnd.gt.0.5) eps(j,i_l)=  1.0
        if (rnd.le.0.5) eps(j,i_l)= -1.0
      enddo

      return
      end subroutine gen_eddy
c-----------------------------------------------------------------------
      real function rnd_loc (lower,upper)
      implicit none

      real, intent(in) :: lower, upper
      real rnd
      double precision, external :: zbqlu01

      rnd = real(zbqlu01(1.0d0))
      rnd_loc = lower + rnd*(upper-lower)
      return
      end function rnd_loc
c-----------------------------------------------------------------------
c     Set energies of nonlocal eddies to zero
      subroutine zero_nonlocal_eddies(ex,ey,ez,eps,neddy_ll,neddy,pt)
      implicit none
c     include 'SIZE'

      integer, intent(in) :: neddy_ll,neddy,pt(1)
      real, intent(inout) :: ex(1),ey(1),ez(1),eps(3,1)
      integer i,k, i1, i2

      i1 = pt(1)-1
      i2 = pt(neddy_ll)+1

      do i=1,i1
         ex(i) = 0.0
         ey(i) = 0.0
         ez(i) = 0.0
         do k=1,3
            eps(k,i)=0.0
         enddo
      enddo
      
      do i=i2,neddy
         ex(i) = 0.0
         ey(i) = 0.0
         ez(i) = 0.0
         do k=1,3
            eps(k,i)=0.0
         enddo
      enddo

      return
      end subroutine zero_nonlocal_eddies

c-----------------------------------------------------------------------
c     Convect eddies and recycle by regenerating the location
      subroutine advect_recycle_eddies(ex,ey,ez,eps,n,pt)
      use SEM, only: zbmax,u0
      implicit none
      include 'SIZE_DEF' ! DT
      include 'SIZE' ! DT
      include 'TSTEP_DEF' ! DT
      include 'TSTEP'

      integer, intent(in) :: n, pt(1)
      real, intent(inout) :: ex(1),ey(1),ez(1),eps(3,1)
     
      integer i,i_l

      do i=1,n
        i_l = pt(i)
        ez(i_l) = ez(i_l) + u0*DT
                            
c     ---- recycle exiting eddies ----
        if (ez(i_l).gt.(zbmax))then
          call gen_eddy(ex,ey,ez,eps,i_l)
        endif
      enddo

      return
      end subroutine advect_recycle_eddies
c-------------------------------------------------------------------
c     Return velocity/turbulence kinetic energy/dissipation rate
c     at given 1D-position, via linear interpolation.

c     In the case of a pipe-geometry that position is the radius of the
c     point of interest

      subroutine SEMinputData(pos_in, vel_out, tke_out, eps_out)
      use SEM, only: nInputdata, pos, umean, tke, dissip
      implicit none

      real, intent(in) :: pos_in
      real, intent(out) :: vel_out, tke_out, eps_out

      integer i
      real a1, a2

      do i=1,nInputdata-1
        if (pos_in.le.pos(i+1)) then
          a1 = (pos(i+1)-pos_in)/(pos(i+1)-pos(i))
          a2 = (pos_in-pos(i))/(pos(i+1)-pos(i))
          vel_out = a1 * umean(i) + a2 * umean(i+1)
          tke_out = a1 * tke(i) + a2 * tke(i+1)
          eps_out = a1 * dissip(i) + a2 * dissip(i+1)
          return
        endif
      enddo

      ! bogus fallback output, to detect interpolation problems

      write(*,*) ' !!! Warning !!! SEM failed at input interpolation'

      vel_out = -1
      tke_out = -1
      eps_out = 0.0

      return
      end subroutine SEMinputData



********  Below follows random-number generation, untouched by Lhuf
********
c-------------------------------------------------------------------
*******************************************************************
********    AUTHORS: Richard Chandler       ***********
********         (richard@stats.ucl.ac.uk)  ***********
********         Paul Northrop              ***********
********         (northrop@stats.ox.ac.uk)  ***********
********    LAST MODIFIED: 26/8/03          ***********
*******************************************************************

      BLOCK DATA ZBQLBD01
*
*       Initializes seed array etc. for random number generator.
*       The values below have themselves been generated using the
*       NAG generator.
*
      COMMON /ZBQL0001/ ZBQLIX,B,C
      DOUBLE PRECISION ZBQLIX(43),B,C
      INTEGER I
      DATA (ZBQLIX(I),I=1,43) /8.001441D7,5.5321801D8,
     +1.69570999D8,2.88589940D8,2.91581871D8,1.03842493D8,
     +7.9952507D7,3.81202335D8,3.11575334D8,4.02878631D8,
     +2.49757109D8,1.15192595D8,2.10629619D8,3.99952890D8,
     +4.12280521D8,1.33873288D8,7.1345525D7,2.23467704D8,
     +2.82934796D8,9.9756750D7,1.68564303D8,2.86817366D8,
     +1.14310713D8,3.47045253D8,9.3762426D7 ,1.09670477D8,
     +3.20029657D8,3.26369301D8,9.441177D6,3.53244738D8,
     +2.44771580D8,1.59804337D8,2.07319904D8,3.37342907D8,
     +3.75423178D8,7.0893571D7 ,4.26059785D8,3.95854390D8,
     +2.0081010D7,5.9250059D7,1.62176640D8,3.20429173D8,
     +2.63576576D8/
      DATA B / 4.294967291D9 /
      DATA C / 0.0D0 /
      END
******************************************************************
******************************************************************
******************************************************************
      SUBROUTINE ZBQLINI(SEED)
******************************************************************
*       To initialize the random number generator - either
*       repeatably or nonrepeatably. Need double precision
*       variables because integer storage can't handle the
*       numbers involved
******************************************************************
*   ARGUMENTS
*   =========
*   SEED    (integer, input). User-input number which generates
*       elements of the array ZBQLIX, which is subsequently used 
*       in the random number generation algorithm. If SEED=0,
*       the array is seeded using the system clock if the 
*       FORTRAN implementation allows it.
******************************************************************
*   PARAMETERS
*   ==========
*   LFLNO   (integer). Number of lowest file handle to try when
*       opening a temporary file to copy the system clock into.
*       Default is 80 to keep out of the way of any existing
*       open files (although the program keeps searching till
*       it finds an available handle). If this causes problems,
*               (which will only happen if handles 80 through 99 are 
*               already in use), decrease the default value.
******************************************************************
      INTEGER LFLNO
      PARAMETER (LFLNO=80)
******************************************************************
*   VARIABLES
*   =========
*   SEED    See above
*   ZBQLIX  Seed array for the random number generator. Defined
*       in ZBQLBD01
*   B,C Used in congruential initialisation of ZBQLIX
*   SS,MM,} System clock secs, mins, hours and days
*   HH,DD }
*   FILNO   File handle used for temporary file
*   INIT    Indicates whether generator has already been initialised
*
      INTEGER SEED,SS,MM,HH,DD,FILNO,I
      INTEGER INIT
      DOUBLE PRECISION ZBQLIX(43),B,C
      DOUBLE PRECISION TMPVAR1,DSS,DMM,DHH,DDD

      COMMON /ZBQL0001/ ZBQLIX,B,C
      SAVE INIT

*
*   Ensure we don't call this more than once in a program
*
      IF (INIT.GE.1) THEN
       IF(INIT.EQ.1) THEN
        WRITE(*,1)
        INIT = 2
       ENDIF
       RETURN
      ELSE
       INIT = 1
      ENDIF
*
*       If SEED = 0, cat the contents of the clock into a file
*       and transform to obtain ZQBLIX(1), then use a congr.
*       algorithm to set remaining elements. Otherwise take
*       specified value of SEED.
*
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*>>>>>>>    NB FOR SYSTEMS WHICH DO NOT SUPPORT THE  >>>>>>>
*>>>>>>>    (NON-STANDARD) 'CALL SYSTEM' COMMAND,    >>>>>>>
*>>>>>>>    THIS WILL NOT WORK, AND THE FIRST CLAUSE >>>>>>>
*>>>>>>>    OF THE FOLLOWING IF BLOCK SHOULD BE  >>>>>>>
*>>>>>>>    COMMENTED OUT.               >>>>>>>
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF (SEED.EQ.0) THEN
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*>>>>>>>    COMMENT OUT FROM HERE IF YOU DON'T HAVE  >>>>>>>
*>>>>>>>    'CALL SYSTEM' CAPABILITY ...         >>>>>>>
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       temp=1.0
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*<<<<<<<<   ... TO HERE (END OF COMMENTING OUT FOR    <<<<<<<
*<<<<<<<<   USERS WITHOUT 'CALL SYSTEM' CAPABILITY    <<<<<<<
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      ELSE
       TMPVAR1 = DMOD(DBLE(SEED),B)
      ENDIF
      ZBQLIX(1) = TMPVAR1
      DO 100 I = 2,43
       TMPVAR1 = ZBQLIX(I-1)*3.0269D4
       TMPVAR1 = DMOD(TMPVAR1,B)       
       ZBQLIX(I) = TMPVAR1
 100  CONTINUE

 1    FORMAT(//5X,'****WARNING**** You have called routine ZBQLINI ',
     +'more than',/5X,'once. I''m ignoring any subsequent calls.',//)
 2    FORMAT(//5X,'**** ERROR **** In routine ZBQLINI, I couldn''t',
     +' find an',/5X,
     +'available file number. To rectify the problem, decrease the ',
     +'value of',/5X,
     +'the parameter LFLNO at the start of this routine (in file ',
     +'randgen.f)',/5X,
     +'and recompile. Any number less than 100 should work.')
      END
******************************************************************
      FUNCTION ZBQLU01(DUMMY)
         implicit none
*
*       Returns a uniform random number between 0 & 1, using
*       a Marsaglia-Zaman type subtract-with-borrow generator.
*       Uses double precision, rather than integer, arithmetic 
*       throughout because MZ's integer constants overflow
*       32-bit integer storage (which goes from -2^31 to 2^31).
*       Ideally, we would explicitly truncate all integer 
*       quantities at each stage to ensure that the double
*       precision representations do not accumulate approximation
*       error; however, on some machines the use of DNINT to
*       accomplish this is *seriously* slow (run-time increased
*       by a factor of about 3). This double precision version 
*       has been tested against an integer implementation that
*       uses long integers (non-standard and, again, slow) -
*       the output was identical up to the 16th decimal place
*       after 10^10 calls, so we're probably OK ...
*
      DOUBLE PRECISION ZBQLU01,DUMMY,B,C,ZBQLIX(43),X,B2,BINV
      INTEGER CURPOS,ID22,ID43

      COMMON /ZBQL0001/ ZBQLIX,B,C
      SAVE /ZBQL0001/
      SAVE CURPOS,ID22,ID43
      DATA CURPOS,ID22,ID43 /1,22,43/

      B2 = B
      BINV = 1.0D0/B
 5    X = ZBQLIX(ID22) - ZBQLIX(ID43) - C
      IF (X.LT.0.0D0) THEN
       X = X + B
       C = 1.0D0
      ELSE
       C = 0.0D0
      ENDIF
      ZBQLIX(ID43) = X
*
*     Update array pointers. Do explicit check for bounds of each to
*     avoid expense of modular arithmetic. If one of them is 0 the others
*     won't be
*
      CURPOS = CURPOS - 1
      ID22 = ID22 - 1
      ID43 = ID43 - 1
      IF (CURPOS.EQ.0) THEN
       CURPOS=43
      ELSEIF (ID22.EQ.0) THEN
       ID22 = 43
      ELSEIF (ID43.EQ.0) THEN
       ID43 = 43
      ENDIF
*
*     The integer arithmetic there can yield X=0, which can cause 
*     problems in subsequent routines (e.g. ZBQLEXP). The problem
*     is simply that X is discrete whereas U is supposed to 
*     be continuous - hence if X is 0, go back and generate another
*     X and return X/B^2 (etc.), which will be uniform on (0,1/B). 
*
      IF (X.LT.BINV) THEN
       B2 = B2*B
       GOTO 5
      ENDIF

      ZBQLU01 = X/B2
      END

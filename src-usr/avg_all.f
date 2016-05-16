C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C
C ---------------------------------------------------------------------- C
C -------------- AVG'ing routines by AZAD   ---------------------------- C
C -------------- quickly copied by Lorenz, added extract Z-Slice ------- C
C ---------------------------------------------------------------------- C
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C
      module AVG
         implicit none
         save

         integer             :: nslices ! number of slices
         real                :: delta_time_avg 
         ! average every time interval
         integer             :: nElperFace ! number of elements per face
         real , allocatable  :: zslices(:)
      end module AVG

C=======================================================================

      subroutine avg_stat_all
        use AVG
        implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'SOLN_DEF'
      include 'SOLN'
      include 'GEOM_DEF'
      include 'GEOM'
      include 'INPUT_DEF'
      include 'INPUT'

      real, external :: glmin, glmax ! defined in math.f
      real, external :: surf_mean ! defined in subs1.f
      external :: comp_vort3 ! defined in navier5.f

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C
C ---------------------------------------------------------------------- C
C -------------- avg_stat_all  ----------------------------------------- C
C This routine computes the 1st, 2nd, 3rd and 4th order statistics in    C
C addition to the terms needed for the Reynolds stress budget.           C
C ---------------------------------------------------------------------- C
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C

      integer nstat, iastep, i, j, k,m, ierr
      integer ntot
      parameter (nstat = 62) ! Number of statistical fields to be saved

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C

      real p0, pmean(lx1,ly1,lz1,lelt)
      common /c_p0/ p0(lx1,ly1,lz1,lelt)

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C
C ---------------------------------------------------------------------- C
C -------------- Define the various quantities on a 3D array------------ C
C ---------------------------------------------------------------------- C
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C

      real stat(lx1*ly1*lz1*lelt, nstat)
      real stat_rot(lx1*ly1*lz1*lelt, nstat)

      real pm1(lx1*ly1*lz1*lelt)
      real wk1(lx1*ly1*lz1)
      real wk2(lx1*ly1*lz1)

      real duidxj(lx1*ly1*lz1, lelt, 1:3*ldim)
      real ur(lx1*ly1*lz1), us(lx1*ly1*lz1), ut(lx1*ly1*lz1)
      real vr(lx1*ly1*lz1), vs(lx1*ly1*lz1), vt(lx1*ly1*lz1)
      real wr(lx1*ly1*lz1), ws(lx1*ly1*lz1), wt(lx1*ly1*lz1)
      real vort(lx1*ly1*lz1*lelv,3),
     $      w3(lx1*ly1*lz1*lelv),w4(lx1*ly1*lz1*lelv)
      real omega_r(lx1*ly1*lz1*lelt), omega_t(lx1*ly1*lz1*lelt)

      real alpha, beta, dtime
      real xlmin, xlmax, domain_x
      real ylmin, ylmax, domain_y
      real zlmin, zlmax, domain_z

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C
C ---------------------------------------------------------------------- C
C -------------- Define the various quantities on a 2D array ----------- C
C ---------------------------------------------------------------------- C
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C

!     real stat_xy(ly1*lely*lz1*lelz, nstat)
      real , allocatable :: stat_xy(:, :), w1(:)

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C

      logical ifverbose
      logical file_exists
      integer icalld
      save    icalld
      data    icalld  /0/

      real atime,timel,times
      common /avgcmnr/ atime,timel

      integer indts, nrec
      save    indts, nrec
      save    times 
      save    domain_x, domain_y, domain_z

      character*80 pippo
      character*80 val1, val2, val3, val4, val5, val6
      character*80 val7, val8, val9, val10
      character*80 inputname1, inputname2, inputname3

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C

      ntot    = nx1*ny1*nz1*nelv
      
      if (icalld.eq.0) then
        icalld = icalld + 1

        atime  = 0.
        timel  = time

        xlmin = glmin(xm1,ntot)
        xlmax = glmax(xm1,ntot)
        domain_x = xlmax - xlmin

        ylmin = glmin(ym1,ntot)          
        ylmax = glmax(ym1,ntot)
        domain_y = ylmax - ylmin

        zlmin = glmin(zm1,ntot)          
        zlmax = glmax(zm1,ntot)
        domain_z = zlmax - zlmin

        call rzero(stat,ntot*nstat)
        ! TODO why not stat_extra and stat_rot?

         nrec  = 0
         times = time

         if(nid.eq.0) then
         
         indts = 0 ! output file counter

         ! Do a stupid search whether there are pre-existing files
         ! (after restarting e.g.)
         if (nslices.gt.0) then
           file_exists = .true.

           do while (file_exists)

           write(pippo,'(F7.3, A, i4.4)') zslices(1), '_',  indts+1
           inputname1 = 'statistics/recordings/stat_z_'//
     $    adjustl(trim(pippo))

           inquire(file=inputname1,exist=file_exists); 

           indts = indts + 1
           end do
           indts = indts - 1
         endif
       endif

      endif
      ifverbose = .FALSE.

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C

      ! fist time avg_stat_all is called this is zero
      dtime = time  - timel

      if (dtime .lt. delta_time_avg) return

      iastep = param(68)
      if  (iastep.eq.0) iastep=param(15)   ! same as iostep
      if  (iastep.eq.0) iastep=500

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C  
     
      atime = atime + dtime

c      call mappr(pm1,pr,wk1,wk2) ! map pressure to mesh 1 (vel. mesh) 
      call copy(p0,pr,ntot)

      pmean = -surf_mean(pr,1,'W  ',ierr)
      call cadd(p0,pmean,ntot)

      call comp_derivat(duidxj,vx,vy,vz,ur,us,ut,vr,vs,vt,wr,ws,wt)
      call comp_vort3(vort,w3,w4,vx,vy,vz)
      call convert_vor(omega_r,omega_t,vort(1,1),vort(1,2),ntot)

      if (atime.ne.0.) then ! closed in l. 1606
        nrec  = nrec + 1
        beta  = dtime/atime
        alpha = 1.-beta
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C 

        call avg1(stat(1,1),vx,alpha,beta,ntot,'velx',ifverbose)     ! <u>
        call avg1(stat(1,2),vy,alpha,beta,ntot,'vely',ifverbose)     ! <v>
        call avg1(stat(1,3),vz,alpha,beta,ntot,'velz',ifverbose)     ! <w>
        call avg1(stat(1,4),p0,alpha,beta,ntot,'pres',ifverbose)     ! <p>

        call avg2(stat(1,5),vx,alpha,beta,ntot,'urms',ifverbose)     ! <uu> (u: instantaneous) 
        call avg2(stat(1,6),vy,alpha,beta,ntot,'vrms',ifverbose)     ! <vv> (v: instantaneous)
        call avg2(stat(1,7),vz,alpha,beta,ntot,'wrms',ifverbose)     ! <ww> (w: instantaneous)
        call avg2(stat(1,8),p0,alpha,beta,ntot,'prms',ifverbose)     ! <pp> (p: instantaneous)

        call avg3(stat(1,9),vx,vy,alpha,beta,ntot,'uvrm',ifverbose)  ! <uv> (u, v: instantaneous)
        call avg3(stat(1,10),vy,vz,alpha,beta,ntot,'vwrm',ifverbose) ! <vw> (v, w: instantaneous)
        call avg3(stat(1,11),vz,vx,alpha,beta,ntot,'wurm',ifverbose) ! <uw> (u, w: instantaneous)

        call avg4(stat(1,12),vx,p0,alpha,beta,ntot,'pu')   ! <pu> (p, u: instantaneous)
        call avg4(stat(1,13),vy,p0,alpha,beta,ntot,'pv')   ! <pv> (p, v: instantaneous)           
        call avg4(stat(1,14),vz,p0,alpha,beta,ntot,'pw')   ! <pw> (p, w: instantaneous)

        call avg5(stat(1,15),p0,duidxj,alpha,beta,ntot,'pux') ! <pdudx> (p, dudx: instantaneous) 
        call avg5(stat(1,16),p0,duidxj,alpha,beta,ntot,'puy') ! <pdudy> (p, dudx: instantaneous) 
        call avg5(stat(1,17),p0,duidxj,alpha,beta,ntot,'puz') ! <pdudz> (p, dudx: instantaneous)

        call avg5(stat(1,18),p0,duidxj,alpha,beta,ntot,'pvx') ! <pdvdx> (p, dvdx: instantaneous) 
        call avg5(stat(1,19),p0,duidxj,alpha,beta,ntot,'pvy') ! <pdvdy> (p, dvdy: instantaneous)
        call avg5(stat(1,20),p0,duidxj,alpha,beta,ntot,'pvz') ! <pdvdz> (p, dudz: instantaneous)   

        call avg5(stat(1,21),p0,duidxj,alpha,beta,ntot,'pwx') ! <pdwdx> (p, dwdx: instantaneous) 
        call avg5(stat(1,22),p0,duidxj,alpha,beta,ntot,'pwy') ! <pdwdy> (p, dwdy: instantaneous)
        call avg5(stat(1,23),p0,duidxj,alpha,beta,ntot,'pwz') ! <pdwdz> (p, dwdz: instantaneous)

        call avg6(stat(1,24),vx,vx,vx,alpha,beta,ntot,'u3')    ! <uuu> (u: instantaneous) 
        call avg6(stat(1,25),vy,vy,vy,alpha,beta,ntot,'v3')    ! <vvv> (v: instantaneous) 
        call avg6(stat(1,26),vz,vz,vz,alpha,beta,ntot,'w3')    ! <www> (w: instantaneous) 
        call avg6(stat(1,27),p0,p0,p0,alpha,beta,ntot,'p3')    ! <ppp> (p: instantaneous) 

        call avg6(stat(1,28),vx,vx,vy,alpha,beta,ntot,'u2v')   ! <uuv> (u, v: instantaneous) 
        call avg6(stat(1,29),vx,vx,vz,alpha,beta,ntot,'u2w')   ! <uuw> (u, w: instantaneous)  
        call avg6(stat(1,30),vy,vy,vx,alpha,beta,ntot,'v2v')   ! <vvu> (v, u: instantaneous)	
        call avg6(stat(1,31),vy,vy,vz,alpha,beta,ntot,'v2w')   ! <vvw> (v, w: instantaneous) 

        call avg6(stat(1,32),vz,vz,vx,alpha,beta,ntot,'w2u')   ! <wwu> (w, u: instantaneous)
        call avg6(stat(1,33),vz,vz,vy,alpha,beta,ntot,'w2v')   ! <wwv> (w, v: instantaneous)
        call avg6(stat(1,34),vx,vy,vz,alpha,beta,ntot,'uvw')   ! <uvw> (u, v, w: instantaneous) 	

        call avg8(stat(1,35),vx,vx,vx,vx,alpha,beta,ntot,'u4')      ! <uuuu> (u: instantaneous)     
        call avg8(stat(1,36),vy,vy,vy,vy,alpha,beta,ntot,'v4')      ! <vvvv> (v: instantaneous)
        call avg8(stat(1,37),vz,vz,vz,vz,alpha,beta,ntot,'w4')      ! <wwww> (w: instantaneous)	 
        call avg8 (stat(1,38),p0,p0,p0,p0,alpha,beta,ntot,'p4')     ! <pppp> (p: instantaneous) 

        call avg8(stat(1,39),vx,vx,vx,vy,alpha,beta,ntot,'u4')      ! <uuuv> (u: instantaneous)     
        call avg8(stat(1,40),vx,vx,vy,vy,alpha,beta,ntot,'v4')      ! <uuvv> (v: instantaneous)
        call avg8(stat(1,41),vx,vy,vy,vy,alpha,beta,ntot,'w4')      ! <uvvv> (w: instantaneous)	 

        call avg7(stat(1,42),duidxj,alpha,beta,ntot,'e11')   ! e11: <du/dx.du/dx + du/dy.du/dy + du/dz.du/dz> (u: instantaneous)
        call avg7(stat(1,43),duidxj,alpha,beta,ntot,'e22')   ! e22: <dv/dx.dv/dx + dv/dy.dv/dy + dv/dz.dv/dz> (v: instantaneous) 
        call avg7(stat(1,44),duidxj,alpha,beta,ntot,'e33')   ! e33: <dw/dx.dw/dx + dw/dy.dw/dy + dw/dz.dw/dz> (w: instantaneous)

        call avg7(stat(1,45),duidxj,alpha,beta,ntot,'e12')   ! e12: <du/dx.dv/dx + du/dy.dv/dy + du/dz.dv/dz> (u, v: instantaneous)       
        call avg7(stat(1,46),duidxj,alpha,beta,ntot,'e13')   ! e13: <du/dx.dw/dx + du/dy.dw/dy + du/dz.dw/dz> (u, w: instantaneous) 
        call avg7(stat(1,47),duidxj,alpha,beta,ntot,'e23')   ! e23: <dv/dx.dw/dx + dv/dy.dw/dy + dv/dz.dw/dz> (v, w: instantaneous) 

        call avg5(stat(1,48),p0,duidxj,alpha,beta,ntot,'nnn') ! <du/dx>
        call avg5(stat(1,49),p0,duidxj,alpha,beta,ntot,'ooo') ! <du/dy>
        call avg5(stat(1,50),p0,duidxj,alpha,beta,ntot,'ppp') ! <du/dz>

        call avg5(stat(1,51),p0,duidxj,alpha,beta,ntot,'qqq') ! <dv/dx>
        call avg5(stat(1,52),p0,duidxj,alpha,beta,ntot,'rrr') ! <dv/dy>
        call avg5(stat(1,53),p0,duidxj,alpha,beta,ntot,'sss') ! <dv/dz>

        call avg5(stat(1,54),p0,duidxj,alpha,beta,ntot,'ttt') ! <dw/dx>
        call avg5(stat(1,55),p0,duidxj,alpha,beta,ntot,'uuu') ! <dw/dy>
        call avg5(stat(1,56),p0,duidxj,alpha,beta,ntot,'vvv') ! <dw/dz>

        call avg1(stat(1,57),omega_r,alpha,beta,ntot,'velx',ifverbose)  ! <omr>
        call avg1(stat(1,58),omega_t,alpha,beta,ntot,'vely',ifverbose)  ! <omt>
        call avg1(stat(1,59),vort(1,3),alpha,beta,ntot,'velz',ifverbose)! <omz>

        call avg2(stat(1,60),omega_r,alpha,beta,ntot,'urms',ifverbose)  ! <omr*omr> 
        call avg2(stat(1,61),omega_t,alpha,beta,ntot,'vrms',ifverbose)  ! <omt*omt> 
        call avg2(stat(1,62),vort(1,3),alpha,beta,ntot,'wrms',ifverbose)!<omz*omz> 


      endif

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C 
C ------------------------------------------------------------------------------- C
C --------- Extract z-slices of the averaged quantities  ------------------------ C
C ------------------------------------------------------------------------------- C
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C 

      if (mod(nrec,iastep).eq.0.and.istep.ge.1) then

        allocate(stat_xy(lx1*lz1*nElperFace, nstat))
        allocate(w1(ly1*lx1*nElperFace))
        
        if(nid.eq.0) indts = indts + 1

        do k=1,nslices

        call extract_z_slice(zslices(k), stat_xy, nElperFace, stat,
     $    nstat,w1) 
        
        if(nid.eq.0) then
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C 
C ------------------------------------------------------------------------------- C
C ------------ Write statistics to file ----------------------------------------- C
C ------------------------------------------------------------------------------- C
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C 
        write(pippo,'(F7.3, A, i4.4)') zslices(k), '_',  indts

        inputname1 = 'statistics/recordings/stat_z_'//
     $    adjustl(trim(pippo))

        write(6,*) inputname1

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C
C ------------------------------------------------------------------------------- C
C ----- Inputname1 -------------------------------------------------------------- C
C ------------------------------------------------------------------------------- C	    
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C

        open(unit=33,form='unformatted',file=inputname1)

        m=lx1*ly1*nElperFace

        write(val1,'(1p,e17.9)') 1./param(2)                 ! Reynolds number	  
        write(val2,'(1p,3e17.9)') domain_x,domain_y,domain_z ! domain size
        write(val3,'(i9)') nElperFace                       ! number of elements 
        write(val4,'(3i9)') nx1-1,ny1-1,nz1-1                ! polynomial order
        write(val5,'(i9)')       nstat                      ! number of saved statistics 
        write(val6,'(1p,e17.9)') times                      ! start time
        write(val7,'(1p,e17.9)') time                       ! end time
        write(val8,'(1p,e17.9)') atime                     ! average time
        write(val9,'(1p,e17.9)') DT                         ! time step
        write(val10,'(i9)')      nrec                       ! number of time records


        write(33) '(Re ='//trim(val1)
     &   //') (Lx, Ly, Lz ='//trim(val2)
     &   //') (nel ='//trim(val3)
     &   //') (Polynomial order ='//trim(val4)
     &   //') (Nstat ='//trim(val5)
     &   //') (start time ='//trim(val6)
     &   //') (end time ='//trim(val7)
     &   //') (average time ='//trim(val8)
     &   //') (time step ='//trim(val9)
     &   //') (nrec ='//trim(val10)
     &   //')'

        write(33) 1/param(2),
     &      domain_x, domain_y, domain_z,
     &      nElperFace,
     &      nx1-1   , ny1-1   , nz1-1,
     &      nstat,
     &      times,
     &      time,
     &      atime,
     &      DT,
     &      nrec


        do i=1,nstat
          write(33) (stat_xy(j,i),j=1,m)  
        enddo

        close(33)

        endif
        enddo

        deallocate(stat_xy)
        deallocate(w1)

        nrec = 0
        times = time
        atime = 0.
      endif

      timel = time


      return
      end subroutine avg_stat_all

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C

c-----------------------------------------------------------------------

!     read parameters AVG
      subroutine AVG_param_in(fid)
        use AVG
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
      namelist /AVG_LIST/ nSlices,nElperFace,delta_time_avg, zslices

!-----------------------------------------------------------------------
!     default values

      nSlices = 50 ! adapt this if necessary
      nelperface = 256
      delta_time_avg = 1e-3
c     zslices = (/0.,1.,2.,5.,10./)
      allocate(zslices(nSlices))
!     read the file
      ierr=0
      if (NID.eq.0) then
         read(unit=fid,nml=AVG_list,iostat=ierr)
      endif
      call err_chk(ierr,'Error reading AVG parameters.$')

!     broadcast data
      call bcast(nSlices,ISIZE)
      call bcast(nelperface,ISIZE)
      call bcast(delta_time_avg,WDSIZE)
      call bcast(zslices, nSlices*WDSIZE)

      return
      end  subroutine AVG_param_in

      subroutine AVG_param_out(fid)
        use AVG
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
      namelist /AVG_LIST/ nSlices,nElperFace,delta_time_avg, zslices
!-----------------------------------------------------------------------
      ierr=0
      if (NID.eq.0) then
         write(unit=fid,nml=AVG_list,iostat=ierr)
      endif
      call err_chk(ierr,'Error writing AVG parameters.$')

      return
      end subroutine AVG_param_out
!***********************************************************************

      subroutine avg4(avg,f,g,alpha,beta,n,name)
c     subroutine avg4(avg,f,g,alpha,beta,n,name,ifverbose)
        implicit none

      real,intent(inout) ::  avg(n)
      real, intent(in) :: alpha, beta, f(n),g(n)
      integer, intent(in) :: n
      character*2,intent(in) :: name
      integer k
c     logical ifverbose

      do k=1,n
         avg(k) = alpha*avg(k) + beta*f(k)*g(k)
      enddo
c
c     if (ifverbose) then
c        avgmax = glmax(avg,n)
c        avgmin = glmin(avg,n)
c        if (nid.eq.0) write(6,1) istep,time,avgmin,avgmax
c    $                           ,alpha,beta,name
c   1    format(i9,1p5e13.5,1x,a4,' av3mnx')
c     endif
c
      return
      end subroutine avg4
c-----------------------------------------------------------------------

      subroutine avg5(avg,f,duidxj,alpha,beta,n,name)
        implicit none
      include 'SIZE_DEF'
      include 'SIZE'

        real,intent(inout) ::  avg(n)
        real, intent(in) :: alpha, beta, f(n),
     $        duidxj(lx1*ly1*lz1,lelt,1:3*ldim)
        integer, intent(in) :: n
        character*3,intent(in) :: name
        integer k

      if (name .eq. 'pux') then
         do k=1,n
            avg(k) = alpha*avg(k) + beta*f(k)*duidxj(k,1,1)
         enddo
      elseif (name .eq. 'pvy') then
         do k=1,n
            avg(k) = alpha*avg(k) + beta*f(k)*duidxj(k,1,2)
         enddo
      elseif (name .eq. 'pwz') then
         do k=1,n
            avg(k) = alpha*avg(k) + beta*f(k)*duidxj(k,1,3)
         enddo
      elseif (name .eq. 'puy') then
         do k=1,n
            avg(k) = alpha*avg(k) + beta*f(k)*duidxj(k,1,4)
         enddo
      elseif (name .eq. 'pvz') then
         do k=1,n
            avg(k) = alpha*avg(k) + beta*f(k)*duidxj(k,1,5)
         enddo
      elseif (name .eq. 'pwx') then
         do k=1,n
            avg(k) = alpha*avg(k) + beta*f(k)*duidxj(k,1,6)
         enddo
      elseif (name .eq. 'puz') then
         do k=1,n
            avg(k) = alpha*avg(k) + beta*f(k)*duidxj(k,1,7)
         enddo
      elseif (name .eq. 'pvx') then
         do k=1,n
            avg(k) = alpha*avg(k) + beta*f(k)*duidxj(k,1,8)
         enddo
      elseif (name .eq. 'pwy') then
         do k=1,n
            avg(k) = alpha*avg(k) + beta*f(k)*duidxj(k,1,9)
         enddo
      elseif (name .eq. 'omz') then  ! <omz>
         do k=1,n
            avg(k) = alpha*avg(k) + beta*(duidxj(k,1,8)-duidxj(k,1,4))
         enddo   
      elseif (name .eq. 'ozz') then  ! <omz*omz>
         do k=1,n
            avg(k) = alpha*avg(k) + beta*
     $     (duidxj(k,1,8)-duidxj(k,1,4))**2
         enddo 
      elseif (name .eq. 'aaa') then  ! <dw/dx*dw/dx>
         do k=1,n
            avg(k) = alpha*avg(k) + beta*duidxj(k,1,6)*duidxj(k,1,6)
         enddo  
      elseif (name .eq. 'bbb') then  ! <dw/dy*dw/dy>
         do k=1,n
            avg(k) = alpha*avg(k) + beta*duidxj(k,1,9)*duidxj(k,1,9)
         enddo 
      elseif (name .eq. 'ccc') then  ! <dw/dx*dw/dy>
         do k=1,n
            avg(k) = alpha*avg(k) + beta*duidxj(k,1,6)*duidxj(k,1,9)
         enddo 
      elseif (name .eq. 'ddd') then  ! <du/dx*du/dx>
         do k=1,n
            avg(k) = alpha*avg(k) + beta*duidxj(k,1,1)*duidxj(k,1,1)
         enddo 
      elseif (name .eq. 'eee') then  ! <du/dy*du/dy>
         do k=1,n
            avg(k) = alpha*avg(k) + beta*duidxj(k,1,4)*duidxj(k,1,4)
         enddo    
      elseif (name .eq. 'fff') then  ! <du/dx*du/dy>
         do k=1,n
            avg(k) = alpha*avg(k) + beta*duidxj(k,1,1)*duidxj(k,1,4)
         enddo  
      elseif (name .eq. 'ggg') then  ! <dv/dx*dv/dx>
         do k=1,n
            avg(k) = alpha*avg(k) + beta*duidxj(k,1,8)*duidxj(k,1,8)
         enddo  
      elseif (name .eq. 'hhh') then  ! <dv/dy*dv/dy>
         do k=1,n
            avg(k) = alpha*avg(k) + beta*duidxj(k,1,2)*duidxj(k,1,2)
         enddo 
      elseif (name .eq. 'iii') then  ! <dv/dx*dv/dy>
         do k=1,n
            avg(k) = alpha*avg(k) + beta*duidxj(k,1,8)*duidxj(k,1,2)
         enddo 
      elseif (name .eq. 'jjj') then  ! <du/dx*dv/dx>
         do k=1,n
            avg(k) = alpha*avg(k) + beta*duidxj(k,1,1)*duidxj(k,1,8)
         enddo 
      elseif (name .eq. 'kkk') then  ! <du/dy*dv/dy>
         do k=1,n
            avg(k) = alpha*avg(k) + beta*duidxj(k,1,4)*duidxj(k,1,2)
         enddo 
      elseif (name .eq. 'lll') then  ! <du/dx*dv/dy>
         do k=1,n
            avg(k) = alpha*avg(k) + beta*duidxj(k,1,1)*duidxj(k,1,2)
         enddo 
      elseif (name .eq. 'mmm') then  ! <du/dy*dv/dx>
         do k=1,n
            avg(k) = alpha*avg(k) + beta*duidxj(k,1,4)*duidxj(k,1,8)
         enddo    
      elseif (name .eq. 'nnn') then  ! <du/dx>
         do k=1,n
            avg(k) = alpha*avg(k) + beta*duidxj(k,1,1)
         enddo  
      elseif (name .eq. 'ooo') then  ! <du/dy>
         do k=1,n
            avg(k) = alpha*avg(k) + beta*duidxj(k,1,4)
         enddo  
      elseif (name .eq. 'ppp') then  ! <du/dz>
         do k=1,n
            avg(k) = alpha*avg(k) + beta*duidxj(k,1,7)
         enddo  
      elseif (name .eq. 'qqq') then  ! <dv/dx>
         do k=1,n
            avg(k) = alpha*avg(k) + beta*duidxj(k,1,8)
         enddo  
      elseif (name .eq. 'rrr') then  ! <dv/dy>
         do k=1,n
            avg(k) = alpha*avg(k) + beta*duidxj(k,1,2)
         enddo  
      elseif (name .eq. 'sss') then  ! <dv/dz>
         do k=1,n
            avg(k) = alpha*avg(k) + beta*duidxj(k,1,5)
         enddo  
      elseif (name .eq. 'ttt') then  ! <dw/dx>
         do k=1,n
            avg(k) = alpha*avg(k) + beta*duidxj(k,1,6)
         enddo 
      elseif (name .eq. 'uuu') then  ! <dw/dy>
         do k=1,n
            avg(k) = alpha*avg(k) + beta*duidxj(k,1,9)
         enddo 
      elseif (name .eq. 'vvv') then  ! <dw/dz>
         do k=1,n
            avg(k) = alpha*avg(k) + beta*duidxj(k,1,3)
         enddo 
      endif
      return
      end subroutine avg5
c----------------------------------------------------------------------

      subroutine avg6(avg,f,g,h,alpha,beta,n,name)
        implicit none

        real,intent(inout) ::  avg(n)
        real, intent(in) :: alpha, beta, f(n), g(n), h(n)
        integer, intent(in) :: n
        character*3,intent(in) :: name
        integer k

        do k=1,n
          avg(k) = alpha*avg(k) + beta*f(k)*g(k)*h(k)
        enddo
      return
      end subroutine avg6
c-----------------------------------------------------------------------

      subroutine avg7(avg,duidxj,alpha,beta,n,name)
        implicit none
      include 'SIZE_DEF'
      include 'SIZE'

        real,intent(inout) ::  avg(n)
        real, intent(in) :: alpha, beta, 
     $      duidxj(lx1*ly1*lz1,lelt,1:3*ldim)
        integer, intent(in) :: n
        character*3,intent(in) :: name
        integer k

      if (name .eq. 'e11') then
         do k=1,n
            avg(k) = alpha*avg(k) + beta*(duidxj(k,1,1)*duidxj(k,1,1) +
     $           duidxj(k,1,4)*duidxj(k,1,4) + 
     $           duidxj(k,1,7)*duidxj(k,1,7))
         enddo
      elseif (name .eq. 'e22') then
         do k=1,n
            avg(k) = alpha*avg(k) + beta*(duidxj(k,1,8)*duidxj(k,1,8) +
     $           duidxj(k,1,2)*duidxj(k,1,2) + 
     $           duidxj(k,1,5)*duidxj(k,1,5))
         enddo
      elseif (name .eq. 'e33') then
         do k=1,n
            avg(k) = alpha*avg(k) + beta*(duidxj(k,1,6)*duidxj(k,1,6) +
     $           duidxj(k,1,9)*duidxj(k,1,9) + 
     $           duidxj(k,1,3)*duidxj(k,1,3))
         enddo
      elseif (name .eq. 'e12') then
         do k=1,n
            avg(k) = alpha*avg(k) + beta*(duidxj(k,1,1)*duidxj(k,1,8) +
     $           duidxj(k,1,4)*duidxj(k,1,2) + 
     $           duidxj(k,1,7)*duidxj(k,1,5))
         enddo
      elseif (name .eq. 'e13') then
         do k=1,n
            avg(k) = alpha*avg(k) + beta*(duidxj(k,1,1)*duidxj(k,1,6) +
     $           duidxj(k,1,4)*duidxj(k,1,9) + 
     $           duidxj(k,1,7)*duidxj(k,1,3))
         enddo
      elseif (name .eq. 'e23') then
         do k=1,n
            avg(k) = alpha*avg(k) + beta*(duidxj(k,1,8)*duidxj(k,1,6) +
     $           duidxj(k,1,2)*duidxj(k,1,9) + 
     $           duidxj(k,1,5)*duidxj(k,1,3))
         enddo
      endif
      return
      end subroutine avg7
c-----------------------------------------------------------------------

      subroutine avg8(avg,f,g,h,s,alpha,beta,n,name)
        implicit none

        real,intent(inout) ::  avg(n)
        real, intent(in) :: alpha, beta, f(n), g(n), h(n), s(n)
        integer, intent(in) :: n
        character*3,intent(in) :: name
        integer k
        do k=1,n
          avg(k) = alpha*avg(k) + beta*f(k)*g(k)*h(k)*s(k)
        enddo
      return
      end subroutine avg8
c-----------------------------------------------------------------------

      subroutine comp_derivat(duidxj,u,v,w,ur,us,ut,vr,vs,vt,wr,ws,wt)
        implicit none
        external :: local_grad3 ! defined in navier1.f
c
      include 'SIZE_DEF'
      include 'SIZE'
      include 'DXYZ_DEF'
      include 'DXYZ'
      include 'GEOM_DEF'
      include 'GEOM'

      real,intent(out) ::  duidxj(lx1*ly1*lz1,lelt,1:3*ldim)
      real, intent(inout) ::  u  (lx1*ly1*lz1,lelt),
     $                     v  (lx1*ly1*lz1,lelt),
     $                     w  (lx1*ly1*lz1,lelt),
     $          ur(lx1*ly1*lz1), us(lx1*ly1*lz1), ut(lx1*ly1*lz1),
     $          vr(lx1*ly1*lz1), vs(lx1*ly1*lz1), vt(lx1*ly1*lz1),
     $          wr(lx1*ly1*lz1), ws(lx1*ly1*lz1), wt(lx1*ly1*lz1)

      integer e,k,n,nxyz

      n    = nx1-1                          ! Polynomial degree
      nxyz = nx1*ny1*nz1
c
      do e=1,nelv
         call local_grad3(ur,us,ut,u,N,e,dxm1,dxtm1)
         call local_grad3(vr,vs,vt,v,N,e,dxm1,dxtm1)
         call local_grad3(wr,ws,wt,w,N,e,dxm1,dxtm1)
c

      do k=1,nxyz
         duidxj(k,e,1) = jacmi(k,e)*(ur(k)*rxm1(k,1,1,e)+
     $        us(k)*sxm1(k,1,1,e)+
     $        ut(k)*txm1(k,1,1,e))
         duidxj(k,e,2) = jacmi(k,e)*(vr(k)*rym1(k,1,1,e)+
     $        vs(k)*sym1(k,1,1,e)+
     $        vt(k)*tym1(k,1,1,e))
         duidxj(k,e,3) = jacmi(k,e)*(wr(k)*rzm1(k,1,1,e)+
     $        ws(k)*szm1(k,1,1,e)+
     $        wt(k)*tzm1(k,1,1,e))
         duidxj(k,e,4) = jacmi(k,e)*(ur(k)*rym1(k,1,1,e)+
     $        us(k)*sym1(k,1,1,e)+
     $        ut(k)*tym1(k,1,1,e))
         duidxj(k,e,5) = jacmi(k,e)*(vr(k)*rzm1(k,1,1,e)+
     $        vs(k)*szm1(k,1,1,e)+
     $        vt(k)*tzm1(k,1,1,e))
         duidxj(k,e,6) = jacmi(k,e)*(wr(k)*rxm1(k,1,1,e)+
     $        ws(k)*sxm1(k,1,1,e)+
     $        wt(k)*txm1(k,1,1,e))
         duidxj(k,e,7) = jacmi(k,e)*(ur(k)*rzm1(k,1,1,e)+
     $        us(k)*szm1(k,1,1,e)+
     $        ut(k)*tzm1(k,1,1,e))
         duidxj(k,e,8) = jacmi(k,e)*(vr(k)*rxm1(k,1,1,e)+
     $        vs(k)*sxm1(k,1,1,e)+
     $        vt(k)*txm1(k,1,1,e))
         duidxj(k,e,9) = jacmi(k,e)*(wr(k)*rym1(k,1,1,e)+
     $        ws(k)*sym1(k,1,1,e)+
     $        wt(k)*tym1(k,1,1,e))
      enddo
      enddo
      return
      end subroutine comp_derivat


      ! Convert vorticity into cylindrical coordinate
      ! adapted from convert_vel I could only find in Azads scratch..

      subroutine convert_vor(or,ot,oxi,oyi,n)
        implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'GEOM_DEF' ! xm1
      include 'GEOM'

      real, intent(out) :: or(n),ot(n)
      real, intent(in) :: oxi(n),oyi(n)
      integer, intent(in) :: n
      integer i
      real  x, y, c, s, prmtrc_t! , prmtrc_r

      do i=1,n  
         x=xm1(i,1,1,1)
         y=ym1(i,1,1,1)
         prmtrc_t=atan2(y,x)
         c = cos(prmtrc_t)
         s = sin(prmtrc_t) 
!        prmtrc_r=sqrt(x*x+y*y)   ! TODO check if neccessary
         or(i) =  c*oxi(i) + s*oyi(i) 
         ot(i) = -s*oxi(i) + c*oyi(i) 
      enddo 
      return
      end subroutine convert_vor
      
C-----------------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C      
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C      
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C 
c-----------------------------------------------------------------------

        ! extract the given z-value (must be mesh aligned)
        ! and contiguous in element indices
      subroutine extract_z_slice(z_val,stat_xy,nelpFac,stat3d,nstat,w1)
        implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'GEOM_DEF'
      include 'GEOM'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
      include 'USERPAR' 

      real angle, circumf, z_unbent

      integer, intent(in) :: nelpFac, nstat
      real, intent(in) :: z_val, stat3d(lx1,ly1,lz1, lelv,nstat)
      real, intent(out) ::  w1(ly1*lz1*nelpFac)
      real, intent(out) :: stat_xy(lx1,ly1,nelpFac, nstat)
      integer e,eg,ex,n, j, k

      call rzero(stat_xy,lx1*ly1*nelpFac*nstat)
      circumf = bent_radius*bent_phi

      do n=1,nstat
      do e=1,nelv
        eg = lglel(e)
        ex = mod(eg-1,nelpFac)+1 ! avoid to access element 0 ..

        z_unbent = zm1(1,1,1,e)

        if (abs(bent_phi).gt.1e-10) then
          if (zm1(1,1,1,e).gt.0) then
            angle = atan2(zm1(1,1,1,e),xm1(1,1,1,e))
            if (angle.le.bent_phi) then
              z_unbent = bent_radius*angle
            else
              z_unbent =
     $  zm1(1,1,1,e)+circumf*(cos(bent_phi)+sin(bent_phi)/cos(bent_phi))
              z_unbent = z_unbent - xm1(1,1,1,e)/cos(bent_phi)
              z_unbent = z_unbent/(cos(bent_phi) + 
     $    sin(bent_phi)*sin(bent_phi)/cos(bent_phi))
            endif
c        elseif (z.le.0.and.x.lt.0) then
        endif
        endif

        ! numerical precission
          if (abs(z_val - z_unbent).lt.1.e-14) then
            do k=1,ly1
              do j=1,lx1
                stat_xy(j,k,ex,n) = stat3d(j,k,1,e,n)
              enddo
            enddo
          endif




      enddo
      call gop(stat_xy(1,1,1,n),w1,'+  ', nelpFac*lx1*ly1)
      enddo

      return
      end subroutine extract_z_slice

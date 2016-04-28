C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C
C ---------------------------------------------------------------------- C
c     dummy module to provide reading interface
c     for parameterfile in 3d case
C ---------------------------------------------------------------------- C
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C
      module AVG_read
         implicit none
         save

         integer             :: nslices ! number of slices
         integer             :: nElperFace ! number of elements per face
         real                :: delta_time_avg
         real , allocatable  :: zslices(:)
      end module AVG_read


!     read parameters AVG
      subroutine AVG_param_in(fid)
        use AVG_read
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
      nSlices = 50
      allocate(zslices(nslices))
      nelperface = 256
c     zslices = (/0.,1.,2.,5.,10./)
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
        use AVG_read
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

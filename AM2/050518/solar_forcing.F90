! [SKC] Custom module for handling solar forcing patterns in the moist model.
! [ZS] Updated 050518

module solar_forcing_mod

use          mpp_mod, only: input_nml_file
use          fms_mod, only: open_namelist_file, check_nml_error, &
                            file_exist, close_file, &
                            error_mesg, FATAL, NOTE, WARNING, mpp_pe, mpp_root_pe, &
                            write_version_number, stdlog, stdout
use       fms_io_mod, only: read_data, field_size
use    constants_mod, only: PI
use horiz_interp_mod, only: horiz_interp_type, horiz_interp_init, &
                            horiz_interp_new, horiz_interp, horiz_interp_del

implicit none
private
!============================================================================================================

character(len=128) :: version = '$Id$'

character(len=128) :: tagname = '$Name$'
character(len=10), parameter :: mod_name='solar_forcing'
integer, parameter :: ipts = 144
integer, parameter :: jpts = 90

!============================================================================================================

public :: solar_forcing

!------- namelist ---------
logical     :: add_solar_forcing = .false.
character(len=32) :: solar_forcing_option = 'gaussian' ! Valid options are 'prescribed' or 'gaussian'; 'prescribed' reads in from a file;
		                                 ! 'gaussian' is computed within the model.
real        :: solar_forcing_mag1 = 0.0
real        :: solar_forcing_center1 = 0.0 ! In degrees (for gaussian forcing only)
real        :: solar_forcing_width1 = 30.0 ! In degrees (for gaussian forcing only)
real        :: solar_forcing_sf1 = 1.0 
real        :: solar_forcing_mag2 = 0.0
real        :: solar_forcing_center2 = 0.0
real        :: solar_forcing_width2 = 30.0
real        :: solar_forcing_sf2 = 1.0
real        :: solar_forcing_mag_prescribed = 1.0

namelist / solar_forcing_nml / add_solar_forcing, solar_forcing_option, &  
                               solar_forcing_mag1, solar_forcing_center1, &
                               solar_forcing_width1, solar_forcing_sf1, &
                               solar_forcing_mag2, solar_forcing_center2, &
                               solar_forcing_width2, solar_forcing_sf2, &
                               solar_forcing_mag_prescribed

!============================================================================================================
contains
!============================================================================================================

subroutine solar_forcing(lonb, latb, lon, lat, forcing)

real, intent(in), dimension(:,:) :: lonb, latb
real, intent(in), dimension(:,:) :: lon, lat
real, intent(out), dimension(:,:) :: forcing

integer :: id_solar_forcing
integer :: io, ierr, nml_unit, stdlog_unit

real, dimension(size(lon,1),size(lat,2))    :: solar_forcing1, solar_forcing2
real  :: c_solar_forcing1, solar_forcing_center1_rad, &
         c_solar_forcing2, solar_forcing_center2_rad

character(len=64) :: file
character(len=256) :: message
integer, dimension(4) :: siz

real, dimension(size(lonb,1)) :: lonb_out
real, dimension(size(latb,2)) :: latb_out
real, allocatable, dimension(:)  :: lonb_in, latb_in
real, allocatable, dimension(:,:) :: forcing_glb
type (horiz_interp_type) :: Interp

call write_version_number(version, tagname)

#ifdef INTERNAL_FILE_NML
   read (input_nml_file, nml=solar_forcing_nml, iostat=io)
   ierr = check_nml_error(io, 'solar_forcing_nml')
#else  
   if ( file_exist('input.nml') ) then
      nml_unit = open_namelist_file()
      ierr = 1
      do while (ierr /= 0
         read (nml_unit, nml=solar_forcing_nml, iostat=io, end=10)
         ierr = check_nml_error(io, 'solar_forcing_nml')
      end do
  10  call close_file(nml_unit)
   endif
#endif

stdlog_unit = stdlog()
if (mpp_pe() == mpp_root_pe()) write(stdlog_unit, solar_forcing_nml)

if(add_solar_forcing) then
   if(solar_forcing_option == 'gaussian') then
      solar_forcing_center1_rad = solar_forcing_center1 * PI / 180.0
      c_solar_forcing1 = solar_forcing_width1 * PI / (2.0 * 180.0 * sqrt(2.0 * log(100.0)))
      solar_forcing1(:,:) =  exp(-(lat(:,:) - solar_forcing_center1_rad)**2 / (2.0 * (c_solar_forcing1**2)))
      solar_forcing_center2_rad = solar_forcing_center2 * PI / 180.0
      c_solar_forcing2 = solar_forcing_width2 * PI / (2.0 * 180.0 * sqrt(2.0 * log(100.0)))
      solar_forcing2(:,:) =  exp(-(lat(:,:) - solar_forcing_center2_rad)**2 / (2.0 * (c_solar_forcing2**2)))
      forcing(:,:) = solar_forcing1 * solar_forcing_mag1 * (1/solar_forcing_sf1) + &
                     solar_forcing2 * solar_forcing_mag2 * (1/solar_forcing_sf2)
   else if(solar_forcing_option == 'prescribed') then
      file = 'INPUT/solar_forcing.nc'
      if(file_exist(trim(file))) then
         if(mpp_pe() == mpp_root_pe()) call error_mesg('solar_forcing_mod','Reading netCDF inputdata: solar_forcing.nc',NOTE)
         call field_size(trim(file), 'solar_forcing', siz)
         if(siz(1) /= ipts .or. siz(2) /= jpts) then
            write(message,*) 'Resolution of solar forcing data does not match model resolution. Solar forcing data: lon_max', &
                        siz(1),', lat_max=', siz(2)
            call error_mesg('solar_forcing_mod', message, FATAL)
         endif
         !if(mpp_pe() == mpp_root_pe()) print *, "size of solar_forcing", siz
         allocate (lonb_in(ipts+1))
         allocate (latb_in(jpts+1))
         allocate (forcing_glb(ipts,jpts))
         call read_data(trim(file), 'lonb', lonb_in, no_domain=.true.)  
         call read_data(trim(file), 'latb', latb_in, no_domain=.true.)  
         call read_data(trim(file), 'solar_forcing', forcing_glb, no_domain=.true.)  
         lonb_in = lonb_in * PI / 180.0
         latb_in = latb_in * PI / 180.0
         lonb_out(:) = lonb(:,1)
         latb_out(:) = latb(1,:)
         !if(mpp_pe() == mpp_root_pe()) print *, "size of lon_in and lat_in", size(lonb_in), size(latb_in)
         !if(mpp_pe() == mpp_root_pe()) print *, "size of lon_out and lat_out", size(lonb_out), size(latb_out)
         !if(mpp_pe() == mpp_root_pe()) print *, "size of forcing_glb", size(forcing_glb,1), size(forcing_glb,2)
         call horiz_interp_init
         call horiz_interp_new(Interp, lonb_in, latb_in, lonb_out, latb_out, interp_method="conservative")
         !if(mpp_pe() == mpp_root_pe()) print *, "size of Interp", Interp%nlon_src, Interp%nlat_src
         call horiz_interp(Interp, forcing_glb, forcing)
         !if(mpp_pe() == mpp_root_pe()) print *, "size of forcing", size(forcing,1), size(forcing,2)
         call horiz_interp_del(Interp)
      else
         call error_mesg('solar_forcing_mod', 'solar forcing file not found for prescribed option; note file must be named solar_forcing.nc', FATAL)
      endif
    
      ! Apply scaling factor from namelist
      forcing = solar_forcing_mag_prescribed * forcing
   else
      call error_mesg('solar_forcing_mod', 'invalid solar_forcing_option provided in namelist; current valid options are gaussian or prescribed', FATAL)
   endif
else
   forcing(:,:) = 0
endif

end subroutine solar_forcing

end module solar_forcing_mod

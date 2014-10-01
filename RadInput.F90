module RadInput

  use shr_kind_mod, only: r8 => shr_kind_r8                                   

  use ppgrid
  use pmgrid
  use error_messages
  use prescribed_aerosols, only: naer_all, &
                                 idxSUL, &
                                 idxSSLT, &
                                 idxDUSTfirst, &
                                 idxOCPHO, &
                                 idxOCPHI, &
                                 idxBCPHO, &
                                 idxBCPHI, &
                                 idxBG, &
                                 idxVOLC

  use ghg_surfvals, only:        co2vmr

  implicit none
  save
!
! Fields input from history file
!
! 3D
!
  real(r8), allocatable :: inp_aermass(:,:,:) ! Aerosol masses per layer
                                              !  Generic name: history vars =
                                              !    MBCPHI_V: BC hydrophilic
                                              !    MBCPHO_V: BC hydrophobic
                                              !    MBG_V   : background (SO4)
                                              !    MDUST1_V: dust size 1
                                              !    MDUST2_V: dust size 2
                                              !    MDUST3_V: dust size 3
                                              !    MDUST4_V: dust size 4
                                              !    MOCPHI_V: OC hydrophilic
                                              !    MOCPHO_V: OC hydrophobic
                                              !    MSSLT_V:  sea salt
                                              !    MSUL_V:   sulfate
                                              !    MVOLC:    volcanic
  real(r8), allocatable :: inp_allaer (:,:,:,:) ! All aerosol masses
  real(r8), allocatable :: inp_cfc11  (:,:,:) ! CFC11 mixing ratio
  real(r8), allocatable :: inp_cfc12  (:,:,:) ! CFC12 mixing ratio
  real(r8), allocatable :: inp_ch4    (:,:,:) ! CH4 mixing ratio
  real(r8), allocatable :: inp_cloud  (:,:,:) ! cloud amount
  real(r8), allocatable :: inp_emis   (:,:,:) ! Cloud longwave emissivity
  real(r8), allocatable :: inp_icldiwp(:,:,:) ! in-cloud ice water path
  real(r8), allocatable :: inp_icldlwp(:,:,:) ! in-cloud total water path
  real(r8), allocatable :: inp_fice(:,:,:)    ! ice fraction
  real(r8), allocatable :: inp_tcldice(:,:,:)    ! ice total grid water box
  real(r8), allocatable :: inp_tcldliq(:,:,:)    ! liq total grid water box
  real(r8), allocatable :: inp_n2o    (:,:,:) ! N2O mixing ratio
  real(r8), allocatable :: inp_o3vmr  (:,:,:) ! Ozone volume mixing ratio
  real(r8), allocatable :: inp_q      (:,:,:) ! Water vapor mixing ratio
  real(r8), allocatable :: inp_rh     (:,:,:) ! Water vapor relative humidity
  real(r8), allocatable :: inp_rel    (:,:,:) ! Liq. droplet eff. radius
  real(r8), allocatable :: inp_rei    (:,:,:) ! Ice effective drop size 
  real(r8), allocatable :: inp_t      (:,:,:) ! Atmospheric temperature
  real(r8), allocatable :: inp_t_cld  (:,:,:) ! Atmospheric temperature for
                                              !   cloud properties
!
! 2D
!
  real(r8), allocatable :: inp_asdir   (:,:)  ! Albedo: shortwave, direct
  real(r8), allocatable :: inp_asdif   (:,:)  ! Albedo: shortwave, diffuse
  real(r8), allocatable :: inp_aldir   (:,:)  ! Albedo: longwave, direct
  real(r8), allocatable :: inp_aldif   (:,:)  ! Albedo: longwave, diffuse
  real(r8), allocatable :: inp_icefrac (:,:)  ! Fractional sea-ice amount
  real(r8), allocatable :: inp_landfrac(:,:)  ! Fractional land amount
  real(r8), allocatable :: inp_landmcos(:,:)  ! Landm coslat field used for rel/rei
  real(r8), allocatable :: inp_ps      (:,:)  ! Surface pressure
  real(r8), allocatable :: inp_snowh   (:,:)  ! Snow height (equivalent water depth)
  real(r8), allocatable :: inp_lwup    (:,:)  ! Upwelling longwave flux
  real(r8), allocatable :: inp_ts      (:,:)  ! Radiative surface temperature (DRF)
  real(r8), allocatable :: inp_taux    (:,:)  ! zonal surface stress (DRF)
  real(r8), allocatable :: inp_tauy    (:,:)  ! meridional surface stress (DRF)
  real(r8), allocatable :: inp_trefht  (:,:)  ! temp at 10m altitude (DRF)
  real(r8), allocatable :: inp_windspeed  (:,:)  ! surface wind speed (DRF)
!
! 1D
!
  real(r8), allocatable :: inp_lat     (:)    ! latitudes
  real(r8), allocatable :: inp_lon     (:)    ! longitudes
!  real, allocatable :: inp_co2vmr      (:)    ! co2vmr DRF

! Indices for aerosol mmmr fields expected by radcswmx
  integer, public, parameter :: aerosol_index(naer_all) = &
      (/ idxSUL, &
         idxSSLT, &
         idxDUSTfirst, &
         idxDUSTfirst+1, &
         idxDUSTfirst+2, &
         idxDUSTfirst+3, &
         idxOCPHO, &
         idxBCPHO, &
         idxOCPHI, &
         idxBCPHI, &
         idxBG, &
         idxVOLC /)

CONTAINS 

  subroutine input_data(path1, path2, fld_option, itime, &
             build_aermmr, build_trace, build_emis, build_re, build_ozone)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Input instantaneous history fields output each radiation time step
!    required to regenerate instantaneous radiation fluxes
!
! Note: Timing data required to reconstitute calday is read using
!    input_times
! 
! Method: 
! Use NetCDF wrapper routines to input a single time slice.
! Data is tranposed onto internal arrays for parallelization in Chunks
! module.
!
! Input fields:
!     ALDIF_r
!     ALDIR_r
!     ASDIF_r
!     ASDIR_r
!     CFC11
!     CFC12
!     CH4
!     CLOUD
!     EMIS
!     ICLDIWP
!     ICLDLWP
!     ICEFRA_r
!     LANDFRAC
!     MBCPHI_V
!     MBCPHO_V
!     MBG_V   
!     MDUST1_V
!     MDUST2_V
!     MDUST3_V
!     MDUST4_V
!     MOCPHI_V
!     MOCPHO_V
!     MSSLT_V
!     MSUL_V
!     MVOLC
!     O3VMR
!     PS
!     Q
!     REI
!     REL
!     SNOWL_r
!     T
!     LWUP_r
!     lat
!     lon 
!     co2vmr  !DRF
! Author: W. Collins
! 
!-----------------------------------------------------------------------

    implicit none
#include <netcdf.inc>
!
! Input arguments
!
    character(len=*), intent(in) :: path1        ! path to data
    character(len=*), intent(in) :: path2        ! path to data
    integer,      intent(in) :: fld_option   ! Option for field switching
                                             !   0 = no switch
                                             !   1 = swap temperatures
                                             !   2 = swap spec. humidity (Q)
                                             !   3 = swap clouds
                                             !   4 = swap albedos
                                             !   5 = perform swaps 1-4

    integer, intent(in) :: itime             ! time slice
    character(len=16), intent(in) :: build_aermmr ! Build AERMMR internally
    logical, intent(in) :: build_trace       ! Build CFCs,CH4,& N2O internally
    logical, intent(in) :: build_emis        ! Build EMIS internally
    logical, intent(in) :: build_re          ! Build RE internally
    logical, intent(in) :: build_ozone       ! Build O3 internally 
!
! Local variables
!
    integer, parameter :: SWAP_TEMP  = 1     ! Flag for temp. swap
    integer, parameter :: SWAP_HUMID = 2     ! Flag for humid. swap
    integer, parameter :: SWAP_CLDS  = 3     ! Flag for temp. swap
    integer, parameter :: SWAP_ALBS  = 4     ! Flag for temp. swap
    integer, parameter :: SWAP_ALL   = 5     ! Flag for swaps 1-4

    integer :: nfid1                         ! NetCDF id for path1
    integer :: nfid2                         ! NetCDF id for path2
    integer :: NFIDX                         ! nfid1 or nfid2
    integer :: start(4)                      ! start indices
    integer :: count(4)                      ! count indices
    integer :: istat                         ! allocate status
    integer :: varid                         ! variable id
    character(len = 10) :: routine = "input_data"

    integer iaer                             ! Aerosol index
    integer ii,jj,kk                         ! DRF indices for recreating IWP
    real(r8) dummy_var                       ! DRF for assigning IWP, LWP values
! Names of aerosol mass fields on history file
    character(len=8), parameter :: aerosol_name(naer_all) =  &
     (/"MSUL_V  "&
      ,"MSSLT_V "&
      ,"MDUST1_V"&
      ,"MDUST2_V"&
      ,"MDUST3_V"&
      ,"MDUST4_V"&
      ,"MOCPHO_V"&
      ,"MBCPHO_V"&
      ,"MOCPHI_V"&
      ,"MBCPHI_V"&
      ,"MBG_V   "&
      ,"MVOLC   "/)


!----------------------------------------------------------------------

    call wrap_open(path1, nf_nowrite, nfid1)
    call wrap_open(path2, nf_nowrite, nfid2)
    
    if (build_aermmr == 'NONE') then
       allocate(inp_aermass(plon,plat,pver), stat = istat) 
       call alloc_err(istat, routine, "inp_aermass",(plon*plat*pver) )

       allocate(inp_allaer(plon,plat,pver,naer_all), stat = istat) 
       call alloc_err(istat, routine, "inp_allaer",(plon*plat*pver*naer_all) )
       
       start = (/1,   1,   1,   itime/) 
       count = (/plon,plat,plev,1    /) 

       do iaer = 1, naer_all
	write(*,*) "iaer = ",iaer
	write(*,*) "aerosol_name(iaer) = ",aerosol_name(iaer)
          call wrap_inq_varid(nfid1, aerosol_name(iaer), varid)
	write(*,*) "nfid1 = ",nfid1
	write(*,*) "varid = ",varid
          call wrap_get_vara_realx(nfid1, varid, start, count, inp_aermass)
          inp_allaer(:,:,:,aerosol_index(iaer)) = inp_aermass
       end do
       
    endif

    if (build_aermmr == 'IPCC') then !!!EVALUATED
       allocate(inp_aermass(plon,plat,pver), stat = istat) 
       call alloc_err(istat, routine, "inp_aermass",(plon*plat*pver) )

       start = (/1,   1,   1,   itime/) 
       count = (/plon,plat,plev,1    /) 

!
! Determine index for sulfates
!
       iaer = minval(minloc(abs(aerosol_index - idxSUL)))
       !call wrap_inq_varid(nfid1, aerosol_name(iaer), varid)
       call wrap_inq_varid(nfid1, "SO4", varid)  !DRF modified 
       call wrap_get_vara_realx(nfid1, varid, start, count, inp_aermass)
       
    endif
    
    if (.not. build_trace) then !!!EVALUATED
       allocate(inp_cfc11(plon,plat,pver), stat = istat) 
       call alloc_err(istat, routine, "inp_cfc11",(plon*plat*pver) )
       start = (/1,   1,   1,   itime/) 
       count = (/plon,plat,pver,1/) 
       call wrap_inq_varid(nfid1, "CFC11", varid)
       call wrap_get_vara_realx(nfid1, varid, start, count, inp_cfc11)
       
       allocate(inp_cfc12(plon,plat,pver), stat = istat) 
       call alloc_err(istat, routine, "inp_cfc12",(plon*plat*pver) )
       start = (/1,   1,   1,   itime/) 
       count = (/plon,plat,pver,1/) 
       call wrap_inq_varid(nfid1, "CFC12", varid)
       call wrap_get_vara_realx(nfid1, varid, start, count, inp_cfc12)
    
       allocate(inp_ch4(plon,plat,pver), stat = istat) 
       call alloc_err(istat, routine, "inp_ch4",(plon*plat*pver) )
       start = (/1,   1,   1,   itime/) 
       count = (/plon,plat,pver,1/) 
       call wrap_inq_varid(nfid1, "CH4", varid)
       call wrap_get_vara_realx(nfid1, varid, start, count, inp_ch4)

       allocate(inp_n2o(plon,plat,pver), stat = istat) 
       call alloc_err(istat, routine, "inp_n2o",(plon*plat*pver) )
       start = (/1,   1,   1,   itime/) 
       count = (/plon,plat,pver,1/) 
       call wrap_inq_varid(nfid1, "N2O", varid)
       call wrap_get_vara_realx(nfid1, varid, start, count, inp_n2o)
    endif

    if (.not. build_emis) then
       allocate(inp_emis(plon,plat,pver), stat = istat) 
       call alloc_err(istat, routine, "inp_emis",(plon*plat*pver) )
       start = (/1,   1,   1,   itime/) 
       count = (/plon,plat,pver,1/) 
       if (fld_option == SWAP_CLDS .or. fld_option == SWAP_ALL) then
          NFIDX = nfid2
       else
          NFIDX = nfid1
       endif
       call wrap_inq_varid(NFIDX, "EMIS", varid)
       call wrap_get_vara_realx(NFIDX, varid, start, count, inp_emis)
    endif
    
    if (.not. build_re) then 
       allocate(inp_rel(plon,plat,pver), stat = istat) 
       call alloc_err(istat, routine, "inp_rel",(plon*plat*pver) )
       start = (/1,   1,   1,   itime/) 
       count = (/plon,plat,pver,1/) 
       if (fld_option == SWAP_CLDS .or. fld_option == SWAP_ALL) then
          NFIDX = nfid2
       else
          NFIDX = nfid1
       endif
       call wrap_inq_varid(NFIDX, "RELL", varid)
       call wrap_get_vara_realx(NFIDX, varid, start, count, inp_rel)
       
       allocate(inp_rei(plon,plat,pver), stat = istat) 
       call alloc_err(istat, routine, "inp_rei",(plon*plat*pver) )
       start = (/1,   1,   1,   itime/) 
       count = (/plon,plat,pver,1/) 
       if (fld_option == SWAP_CLDS .or. fld_option == SWAP_ALL) then
          NFIDX = nfid2
       else
          NFIDX = nfid1
       endif
       call wrap_inq_varid(NFIDX, "REI", varid)
       call wrap_get_vara_realx(NFIDX, varid, start, count, inp_rei)
       
    else !DRF !!!EVALUATED
       allocate(inp_landfrac(plon,plat), stat = istat) 
       call alloc_err(istat, routine, "inp_landfrac",(plon*plat) )
       start = (/1,   1,   itime, -1/) 
       count = (/plon,plat,1,     -1/) 
       call wrap_inq_varid(nfid1, "LANDFRAC", varid)
       call wrap_get_vara_realx(nfid1, varid, start, count, inp_landfrac)
       
       allocate(inp_icefrac(plon,plat), stat = istat) 
       call alloc_err(istat, routine, "inp_icefrac",(plon*plat) )
       start = (/1,   1,   itime, -1/) 
       count = (/plon,plat,1,     -1/) 
       !DRF call wrap_inq_varid(nfid1, "ICEFRA_r", varid)
       call wrap_inq_varid(nfid1, "ICEFRAC", varid)
       call wrap_get_vara_realx(nfid1, varid, start, count, inp_icefrac)

       allocate(inp_snowh(plon,plat), stat = istat) 
       call alloc_err(istat, routine, "inp_snowh",(plon*plat) )
       start = (/1,   1,   itime, -1/) 
       count = (/plon,plat,1,     -1/) 
       !DRF call wrap_inq_varid(nfid1, "SNOWHL_r", varid)
       call wrap_inq_varid(nfid1, "SNOWHLND", varid)
       call wrap_get_vara_realx(nfid1, varid, start, count, inp_snowh)
    endif
!!!EVALUATED
    allocate(inp_taux(plon,plat), stat = istat)
    call alloc_err(istat, routine, "inp_taux",(plon*plat) )
    start = (/1,   1,   itime, -1/) 
    count = (/plon,plat,1    , -1/) 
    call wrap_inq_varid(nfid1, "TAUX", varid)
    call wrap_get_vara_realx(nfid1, varid, start, count, inp_taux)

    allocate(inp_tauy(plon,plat), stat = istat)
    call alloc_err(istat, routine, "inp_tauy",(plon*plat) )
    start = (/1,   1,   itime, -1/) 
    count = (/plon,plat,1    , -1/) 
    call wrap_inq_varid(nfid1, "TAUX", varid)
    call wrap_get_vara_realx(nfid1, varid, start, count, inp_tauy)

    allocate(inp_trefht(plon,plat), stat = istat)
    call alloc_err(istat, routine, "inp_trefht",(plon*plat) )
    start = (/1,   1,   itime, -1/) 
    count = (/plon,plat,1    , -1/) 
    call wrap_inq_varid(nfid1, "TREFHT", varid)
    call wrap_get_vara_realx(nfid1, varid, start, count, inp_trefht)

    allocate(inp_asdir(plon,plat), stat = istat) 
    call alloc_err(istat, routine, "inp_asdir",(plon*plat) )
    start = (/1,   1,   itime, -1/) 
    count = (/plon,plat,1    , -1/) 
    if (fld_option == SWAP_ALBS .or. fld_option == SWAP_ALL) then
       NFIDX = nfid2
    else
       NFIDX = nfid1
    endif
    !DRF call wrap_inq_varid(NFIDX, "ASDIR_r", varid)
    !DRF call wrap_get_vara_realx(NFIDX, varid, start, count, inp_asdir)
    
    allocate(inp_asdif(plon,plat), stat = istat) 
    call alloc_err(istat, routine, "inp_asdif",(plon*plat) )
    start = (/1,   1,   itime, -1/) 
    count = (/plon,plat,1    , -1/) 
    if (fld_option == SWAP_ALBS .or. fld_option == SWAP_ALL) then
       NFIDX = nfid2
    else
       NFIDX = nfid1
    endif
    !DRF call wrap_inq_varid(NFIDX, "ASDIF_r", varid)
    !DRF call wrap_get_vara_realx(NFIDX, varid, start, count, inp_asdif)
    
    allocate(inp_aldir(plon,plat), stat = istat) 
    call alloc_err(istat, routine, "inp_aldir",(plon*plat) )
    start = (/1,   1,   itime, -1/) 
    count = (/plon,plat,1    , -1/) 
    if (fld_option == SWAP_ALBS .or. fld_option == SWAP_ALL) then
       NFIDX = nfid2
    else
       NFIDX = nfid1
    endif
    !DRF call wrap_inq_varid(NFIDX, "ALDIR_r", varid)
    !DRF call wrap_get_vara_realx(NFIDX, varid, start, count, inp_aldir)
    
    allocate(inp_aldif(plon,plat), stat = istat) 
    call alloc_err(istat, routine, "inp_aldif",(plon*plat) )
    start = (/1,   1,   itime, -1/) 
    count = (/plon,plat,1    , -1/) 
    if (fld_option == SWAP_ALBS .or. fld_option == SWAP_ALL) then
       NFIDX = nfid2
    else
       NFIDX = nfid1
    endif
    !DRF call wrap_inq_varid(NFIDX, "ALDIF_r", varid)
    !DRF call wrap_get_vara_realx(NFIDX, varid, start, count, inp_aldif)
    
    allocate(inp_cloud(plon,plat,pver), stat = istat) 
    call alloc_err(istat, routine, "inp_cloud",(plon*plat*pver) )
    start = (/1,   1,   1,   itime/) 
    count = (/plon,plat,pver,1/) 
    if (fld_option == SWAP_CLDS .or. fld_option == SWAP_ALL) then
       NFIDX = nfid2
    else
       NFIDX = nfid1
    endif
    call wrap_inq_varid(NFIDX, "CLOUD", varid)
    call wrap_get_vara_realx(NFIDX, varid, start, count, inp_cloud)
    
    allocate(inp_icldiwp(plon,plat,pver), stat = istat) 
    call alloc_err(istat, routine, "inp_icldiwp",(plon*plat*pver) )
    allocate(inp_fice(plon,plat,pver), stat = istat) 
    call alloc_err(istat, routine, "inp_fice",(plon*plat*pver) )
    start = (/1,   1,   1,   itime/) 
    count = (/plon,plat,pver,1/) 
    if (fld_option == SWAP_CLDS .or. fld_option == SWAP_ALL) then
       NFIDX = nfid2
    else
       NFIDX = nfid1
    endif
    !DRF call wrap_inq_varid(NFIDX, "ICLDIWP", varid)
    call wrap_inq_varid(NFIDX, "FICE", varid)
    call wrap_get_vara_realx(NFIDX, varid, start, count, inp_fice)

    allocate(inp_tcldice(plon,plat,pver), stat = istat) 
    call alloc_err(istat, routine, "inp_tcldice",(plon*plat*pver) )
    call wrap_inq_varid(NFIDX, "CLDICE", varid)
    call wrap_get_vara_realx(NFIDX, varid, start, count, inp_tcldice)
    allocate(inp_tcldliq(plon,plat,pver), stat = istat) 
    call alloc_err(istat, routine, "inp_tcldliq",(plon*plat*pver) )
    call wrap_inq_varid(NFIDX, "CLDLIQ", varid)
    call wrap_get_vara_realx(NFIDX, varid, start, count, inp_tcldliq)
    
    allocate(inp_icldlwp(plon,plat,pver), stat = istat) 
    call alloc_err(istat, routine, "inp_icldlwp",(plon*plat*pver) )
    start = (/1,   1,   1,   itime/) 
    count = (/plon,plat,pver,1/) 
    if (fld_option == SWAP_CLDS .or. fld_option == SWAP_ALL) then
       NFIDX = nfid2
    else
       NFIDX = nfid1
    endif
    call wrap_inq_varid(NFIDX, "ICLDLWP", varid)
    call wrap_get_vara_realx(NFIDX, varid, start, count, inp_icldlwp)
   
    do ii=1,plon
       do jj=1,plat
          do kk=1,pver

             if (inp_icldlwp(ii,jj,kk).gt.0.0_r8) then
               inp_icldiwp(ii,jj,kk) = inp_icldlwp(ii,jj,kk)*inp_tcldice(ii,jj,kk)/(inp_tcldliq(ii,jj,kk)+inp_tcldice(ii,jj,kk))
             else
                inp_icldiwp(ii,jj,kk) = 0.0_r8
             endif
             inp_icldlwp(ii,jj,kk)=inp_icldlwp(ii,jj,kk)+inp_icldiwp(ii,jj,kk)

          end do
       end do
    end do

    if (.not. build_ozone) then
       allocate(inp_o3vmr(plon,plat,pver), stat = istat) 
       call alloc_err(istat, routine, "inp_o3vmr",(plon*plat*pver) )
       start = (/1,   1,   1,   itime/) 
       count = (/plon,plat,pver,1/) 
       !DRF call wrap_inq_varid(nfid1, "O3VMR", varid)
       call wrap_inq_varid(nfid1, "O3", varid)
       call wrap_get_vara_realx(nfid1, varid, start, count, inp_o3vmr)
    endif
    
    allocate(inp_ps(plon,plat), stat = istat) 
    call alloc_err(istat, routine, "inp_ps",(plon*plat) )
    start = (/1,   1,   itime, -1/) 
    count = (/plon,plat,1,     -1/) 
    call wrap_inq_varid(nfid1, "PS", varid)
    call wrap_get_vara_realx(nfid1, varid, start, count, inp_ps)

    !DRF
    allocate(inp_windspeed(plon,plat), stat = istat)
    call alloc_err(istat, routine, "inp_windspeed",(plon*plat) )
    do ii=1,plon
       do jj=1,plat
          inp_windspeed(ii,jj) = sqrt(287.06/0.002*(abs(inp_taux(ii,jj))+abs(inp_tauy(ii,jj)))*inp_trefht(ii,jj)/inp_ps(ii,jj))
       end do
    end do

    !write(*,*) "inp_windspeed = ",inp_windspeed(128,95)
    !write(*,*) "inp_taux= ",inp_taux(30,30)
    !write(*,*) "inp_tauy= ",inp_tauy(30,30)
    !write(*,*) "inp_ps= ",inp_ps(30,30)
    !write(*,*) "inp_trefht= ",inp_trefht(30,30)
    !stop 
    !END DRF
    
    allocate(inp_q(plon,plat,pver), stat = istat) 
    call alloc_err(istat, routine, "inp_q",(plon*plat*pver) )
    start = (/1,   1,   1,   itime/) 
    count = (/plon,plat,pver,1/) 
    if (fld_option == SWAP_HUMID .or. fld_option == SWAP_ALL) then
       NFIDX = nfid2
    else
       NFIDX = nfid1
    endif
    call wrap_inq_varid(NFIDX, "Q", varid)
    call wrap_get_vara_realx(NFIDX, varid, start, count, inp_q)

    allocate(inp_rh(plon,plat,pver), stat = istat) 
    call alloc_err(istat, routine, "inp_rh",(plon*plat*pver) )
    start = (/1,   1,   1,   itime/) 
    count = (/plon,plat,pver,1/) 
    if (fld_option == SWAP_HUMID .or. fld_option == SWAP_ALL) then
       NFIDX = nfid2
    else
       NFIDX = nfid1
    endif
    call wrap_inq_varid(NFIDX, "RELHUM", varid)
    call wrap_get_vara_realx(NFIDX, varid, start, count, inp_rh)



    
    allocate(inp_t(plon,plat,pver), stat = istat) 
    call alloc_err(istat, routine, "inp_t",(plon*plat*pver) )
    start = (/1,   1,   1,   itime/) 
    count = (/plon,plat,pver,1/) 
    if (fld_option == SWAP_TEMP .or. fld_option == SWAP_ALL) then
       NFIDX = nfid2
    else
       NFIDX = nfid1
    endif
    call wrap_inq_varid(NFIDX, "T", varid)
    call wrap_get_vara_realx(NFIDX, varid, start, count, inp_t)
    
    allocate(inp_t_cld(plon,plat,pver), stat = istat) 
    call alloc_err(istat, routine, "inp_t_cld", (plon*plat*pver) )
    start = (/1,   1,   1,   itime/) 
    count = (/plon,plat,pver,1/) 
    if (fld_option == SWAP_CLDS .or. fld_option == SWAP_ALL) then
       NFIDX = nfid2
    else
       NFIDX = nfid1
    endif
    call wrap_inq_varid(NFIDX, "T", varid)
    call wrap_get_vara_realx(NFIDX, varid, start, count, inp_t_cld)
    
    allocate(inp_lwup(plon,plat), stat = istat) 
    call alloc_err(istat, routine, "inp_lwup",(plon*plat) )
    allocate(inp_ts(plon,plat), stat = istat) 
    call alloc_err(istat, routine, "inp_ts",(plon*plat) )
    start = (/1,   1,   itime, -1/)
    count = (/plon,plat,1    , -1/)
    if (fld_option == SWAP_TEMP .or. fld_option == SWAP_ALL) then
       NFIDX = nfid2
    else
       NFIDX = nfid1
    endif
    !DRF call wrap_inq_varid(NFIDX, "LWUP_r", varid)
    !DRF call wrap_get_vara_realx(NFIDX, varid, start, count, inp_lwup)
    call wrap_inq_varid(NFIDX, "TS", varid)
    call wrap_get_vara_realx(NFIDX, varid, start, count, inp_ts)

    allocate(inp_lat(plat), stat = istat) 
    call alloc_err(istat, routine, "inp_lat",(plat) )
    call wrap_inq_varid(nfid1, "lat", varid)
    call wrap_get_var_realx(nfid1, varid, inp_lat)

    allocate(inp_lon(plon), stat = istat) 
    call alloc_err(istat, routine, "inp_lon",(plon) )
    call wrap_inq_varid(nfid1, "lon", varid)
    call wrap_get_var_realx(nfid1, varid, inp_lon)

    start = (/itime, -1, -1, -1/)
    count = (/1    , -1, -1, -1/)
    call wrap_inq_varid(nfid1, "co2vmr", varid)
    call wrap_get_var_realx(nfid1, varid, co2vmr)
   
    call wrap_close(nfid1)
    call wrap_close(nfid2)

    return

  end subroutine input_data

  subroutine input_times(path, nslice, nstep   ,dtime   ,mdbase  ,msbase  ,&
                         mbdate, mbsec, date, datesec)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Input timing data required to reconstitute calday
! 
! Method: 
! Use NetCDF wrapper routines to input base day, date, and corresponding
!     time of day, along with time step information
! Output corresponds to following fields in history files:
!     ndbase
!     nsbase
!     nbdate
!     nbsec
!     nsteph
!     mdt 
!     date
!     datesec
!
! Author: W. Collins
! 
!-----------------------------------------------------------------------
    implicit none
#include <netcdf.inc>

!
! Input arguments
!
    character(len=*), intent(in) :: path     ! path to data
!
! Output arguments
!
    integer, intent(out) ::  nslice      ! number of time slices in file
    integer, pointer     ::  nstep(:)    ! current time step
    integer, pointer     ::  date(:)     ! current date
    integer, pointer     ::  datesec(:)  ! current second in date
    real(r8), intent(out) :: dtime       ! length of time step (seconds)
    integer, intent(out) ::  mdbase      ! base day of run (e.g., 0)
    integer, intent(out) ::  msbase      ! base seconds of base day (e.g., 0)
    integer, intent(out) ::  mbdate      ! base date (yyyymmdd format) of run
    integer, intent(out) ::  mbsec       ! base seconds of base date (e.g., 0)
!
! Local variables
!
    integer :: nfid                      ! NetCDF id
    integer :: istat                     ! allocate status
    integer :: varid                     ! variable id
    integer :: dimid                     ! dimension id
    integer :: mdt                       ! Integer time step
    character (len = (nf_max_name)) :: dimname   ! dimension name
    character (len = 11) :: routine = "input_times"
    
!----------------------------------------------------------------------

    call wrap_open(path, nf_nowrite, nfid)

!
! Get number of time slices
!
    call wrap_inq_dimid(nfid, "time", dimid)
    call wrap_inq_dim(nfid, dimid, dimname, nslice)
!
! Allocate space for time steps, date, datesec
!
    allocate(nstep(nslice), stat = istat)
    call alloc_err(istat, routine, "nstep", nslice)

    allocate(date(nslice), stat = istat)
    call alloc_err(istat, routine, "date", nslice)

    allocate(datesec(nslice), stat = istat)
    call alloc_err(istat, routine, "datesec", nslice)
!
! Get the data
!
    call wrap_inq_varid(nfid, "ndbase", varid)
    call wrap_get_var_int(nfid, varid, mdbase)

    call wrap_inq_varid(nfid, "nsbase", varid)
    call wrap_get_var_int(nfid, varid, msbase)

    call wrap_inq_varid(nfid, "nbdate", varid)
    call wrap_get_var_int(nfid, varid, mbdate)

    call wrap_inq_varid(nfid, "nbsec", varid)
    call wrap_get_var_int(nfid, varid, mbsec)

    call wrap_inq_varid(nfid, "nsteph", varid)
    call wrap_get_var_int(nfid, varid, nstep)

    call wrap_inq_varid(nfid, "mdt", varid)
    call wrap_get_var_int(nfid, varid, mdt)
    dtime = mdt

    call wrap_inq_varid(nfid, "date", varid)
    call wrap_get_var_int(nfid, varid, date)

    call wrap_inq_varid(nfid, "datesec", varid)
    call wrap_get_var_int(nfid, varid, datesec)

    call wrap_close(nfid)
    
    return

  end subroutine input_times

  subroutine input_landmcos(ncdata)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Input landm_coslat field from initial conditions files
! 
! Method: 
! Use NetCDF wrapper routines to input landm_coslat
! Output corresponds to following fields in initial files:
!     landm_coslat
!
! Author: W. Collins
! 
!-----------------------------------------------------------------------
    implicit none
#include <netcdf.inc>

!
! Input arguments
!
    character(len=*), intent(in) :: ncdata     ! path to data
!
! Output arguments
!
!
! Local variables
!
    integer :: nfid                      ! NetCDF id
    integer :: istat                     ! allocate status
    integer :: varid                     ! variable id
    integer :: start(4)                  ! start indices
    integer :: count(4)                  ! count indices

    character (len = 14) :: routine = "input_landmcos"
    
!----------------------------------------------------------------------

    call wrap_open(ncdata, nf_nowrite, nfid)

!
! Allocate space for landm_coslat
!
    allocate(inp_landmcos(plon,plat), stat = istat) 
    call alloc_err(istat, routine, "inp_landmcos",(plon*plat) )
    start = (/1,   1,   1, -1/) 
    count = (/plon,plat,1    , -1/) 
!
! Get the landm_coslat data
!
    call wrap_inq_varid(nfid, "LANDM_COSLAT", varid)
    call wrap_get_vara_realx(nfid, varid, start, count, inp_landmcos)

    call wrap_close(nfid)
    
    return

  end subroutine input_landmcos

  subroutine input_vert_grid(path)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Reads in vertical grid data
! 
! Method: 
! Acquires following fields
!    hyai
!    hybi
!    hyam
!    hybm
!    p0
! 
! Author: W. Collins
! 
!-----------------------------------------------------------------------
    implicit none
#include <netcdf.inc>
#include <comhyb.h>
!
! Input arguments
!
    character(len=*), intent(in) :: path     ! path to data
!
! Local variables
!
    integer :: nfid                      ! NetCDF id
    integer :: varid                     ! variable id
    character (len = 15) :: routine = "input_vert_grid"

!
! Get vertical coordinate info in comhyb common block for constructing
!     pressure fields
!
    call wrap_open(path, nf_nowrite, nfid)
   
    call wrap_inq_varid(nfid, "hyai", varid)
    call wrap_get_var_realx(nfid, varid, hyai)

    call wrap_inq_varid(nfid, "hyam", varid)
    call wrap_get_var_realx(nfid, varid, hyam)

    call wrap_inq_varid(nfid, "hybi", varid)
    call wrap_get_var_realx(nfid, varid, hybi)

    call wrap_inq_varid(nfid, "hybm", varid)
    call wrap_get_var_realx(nfid, varid, hybm)

    call wrap_inq_varid(nfid, "P0", varid)
    call wrap_get_var_realx(nfid, varid, ps0)

    call wrap_close(nfid)

    return

 end subroutine input_vert_grid

  subroutine dump_input_data
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Dump allocated arrays for input data 
! 
! Method: 
! Standard F90 deallocate calls
! 
! Author: W. Collins
! 
!-----------------------------------------------------------------------

!
! Local variables
!
    integer :: istat           ! status variable
    character (len = 15) :: routine = "dump_input_data"

    if (allocated(inp_aermass))   then 
       deallocate(inp_aermass, stat = istat)
       call dealloc_err(istat, routine, "inp_aermass")
    endif
    if (allocated(inp_allaer))   then 
       deallocate(inp_allaer, stat = istat)
       call dealloc_err(istat, routine, "inp_allaer")
    endif
    if (allocated(inp_asdir))    then 
       deallocate(inp_asdir, stat = istat)
       call dealloc_err(istat, routine, "inp_asdir")
    endif
    if (allocated(inp_asdif))    then 
       deallocate(inp_asdif, stat = istat)
       call dealloc_err(istat, routine, "inp_asdif")
    endif
    if (allocated(inp_aldir))    then 
       deallocate(inp_aldir, stat = istat)
       call dealloc_err(istat, routine, "inp_aldir")
    endif
    if (allocated(inp_aldif))    then 
       deallocate(inp_aldif, stat = istat)
       call dealloc_err(istat, routine, "inp_aldif")
    endif
    if (allocated(inp_cfc11))    then 
       deallocate(inp_cfc11, stat = istat)
       call dealloc_err(istat, routine, "inp_cfc11")
    endif
    if (allocated(inp_cfc12))    then 
       deallocate(inp_cfc12, stat = istat)
       call dealloc_err(istat, routine, "inp_cfc12")
    endif
    if (allocated(inp_ch4))      then 
       deallocate(inp_ch4, stat = istat)
       call dealloc_err(istat, routine, "inp_ch4")
    endif
    if (allocated(inp_cloud))    then 
       deallocate(inp_cloud, stat = istat)
       call dealloc_err(istat, routine, "inp_cloud")
    endif
    if (allocated(inp_emis))     then 
       deallocate(inp_emis, stat = istat)
       call dealloc_err(istat, routine, "inp_emis")
    endif
    if (allocated(inp_icldiwp))  then 
       deallocate(inp_icldiwp, stat = istat)
       call dealloc_err(istat, routine, "inp_icldiwp")
    endif
    if (allocated(inp_icldlwp))  then 
       deallocate(inp_icldlwp, stat = istat)
       call dealloc_err(istat, routine, "inp_icldlwp")
    endif
    if (allocated(inp_fice))  then 
       deallocate(inp_fice, stat = istat)
       call dealloc_err(istat, routine, "inp_fice")
    endif
    if (allocated(inp_tcldice))  then 
       deallocate(inp_tcldice, stat = istat)
       call dealloc_err(istat, routine, "inp_tcldice")
    endif
    if (allocated(inp_tcldliq))  then 
       deallocate(inp_tcldliq, stat = istat)
       call dealloc_err(istat, routine, "inp_tcldliq")
    endif
    if (allocated(inp_landfrac)) then 
       deallocate(inp_landfrac, stat = istat)
       call dealloc_err(istat, routine, "inp_landfrac")
    endif
    if (allocated(inp_landmcos)) then 
       deallocate(inp_landmcos, stat = istat)
       call dealloc_err(istat, routine, "inp_landmcos")
    endif
    if (allocated(inp_icefrac)) then 
       deallocate(inp_icefrac, stat = istat)
       call dealloc_err(istat, routine, "inp_icefrac")
    endif
    if (allocated(inp_snowh)) then 
       deallocate(inp_snowh, stat = istat)
       call dealloc_err(istat, routine, "inp_snowh")
    endif
    if (allocated(inp_n2o))      then 
       deallocate(inp_n2o, stat = istat)
       call dealloc_err(istat, routine, "inp_n2o")
    endif
    if (allocated(inp_o3vmr))    then 
       deallocate(inp_o3vmr, stat = istat)
       call dealloc_err(istat, routine, "inp_o3vmr")
    endif
    if (allocated(inp_ps))       then 
       deallocate(inp_ps, stat = istat)
       call dealloc_err(istat, routine, "inp_ps")
    endif
    if (allocated(inp_q))        then 
       deallocate(inp_q, stat = istat)
       call dealloc_err(istat, routine, "inp_q")
    endif
    if (allocated(inp_rh))        then 
       deallocate(inp_rh, stat = istat)
       call dealloc_err(istat, routine, "inp_rh")
    endif
    if (allocated(inp_rel))      then 
       deallocate(inp_rel, stat = istat)
       call dealloc_err(istat, routine, "inp_rel")
    endif
    if (allocated(inp_rei))      then 
       deallocate(inp_rei, stat = istat)
       call dealloc_err(istat, routine, "inp_rei")
    endif
    if (allocated(inp_t))        then 
       deallocate(inp_t, stat = istat)
       call dealloc_err(istat, routine, "inp_t")
    endif
    if (allocated(inp_t_cld))        then 
       deallocate(inp_t_cld, stat = istat)
       call dealloc_err(istat, routine, "inp_t_cld")
    endif
    if (allocated(inp_lwup))       then 
       deallocate(inp_lwup, stat = istat)
       call dealloc_err(istat, routine, "inp_lwup")
    endif
    if (allocated(inp_ts))       then 
       deallocate(inp_ts, stat = istat)
       call dealloc_err(istat, routine, "inp_ts")
    endif
    if (allocated(inp_lat))       then 
       deallocate(inp_lat, stat = istat)
       call dealloc_err(istat, routine, "inp_lat")
    endif
    if (allocated(inp_lon))       then 
       deallocate(inp_lon, stat = istat)
       call dealloc_err(istat, routine, "inp_lon")
    endif

    return

  end subroutine dump_input_data
    
end module RadInput

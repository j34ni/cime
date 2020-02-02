module dshr_strdata_mod

  ! holds data and methods to advance data models
  ! TODO: how to handle  z dimension?
  ! TODO: add scmlon and scmlat functionality

  use ESMF
  use NUOPC
  use shr_kind_mod    , only : in=>shr_kind_in, r8=>shr_kind_r8, cs=>shr_kind_cs, cl=>shr_kind_cl, cxx=>shr_kind_cxx
  use shr_sys_mod     , only : shr_sys_abort, shr_sys_flush
  use shr_const_mod   , only : shr_const_pi
  use shr_mpi_mod     , only : shr_mpi_bcast
  use shr_cal_mod     , only : shr_cal_calendarname, shr_cal_timeSet
  use shr_cal_mod     , only : shr_cal_noleap, shr_cal_gregorian
  use shr_cal_mod     , only : shr_cal_date2ymd, shr_cal_ymd2date
  use shr_orb_mod     , only : shr_orb_decl, shr_orb_cosz, shr_orb_undef_real
  use shr_nl_mod      , only : shr_nl_find_group_name
  use shr_mpi_mod     , only : shr_mpi_bcast
  use shr_const_mod   , only : shr_const_cDay
  use shr_string_mod  , only : shr_string_listGetNum, shr_string_listGetName
  use shr_tinterp_mod , only : shr_tInterp_getAvgCosz
  use shr_tinterp_mod , only : shr_tInterp_getFactors
  use dshr_stream_mod , only : shr_stream_taxis_cycle
  use dshr_stream_mod , only : shr_stream_taxis_extend
  use dshr_stream_mod , only : shr_stream_taxis_limit
  use dshr_stream_mod , only : shr_stream_file_null
  use dshr_stream_mod , only : shr_stream_streamType
  use dshr_stream_mod , only : dshr_stream_findBounds
  use dshr_stream_mod , only : dshr_stream_getFilePath
  use dshr_stream_mod , only : dshr_stream_getPrevFileName
  use dshr_stream_mod , only : dshr_stream_getNextFileName
  use dshr_stream_mod , only : dshr_stream_getFileFieldName
  use dshr_stream_mod , only : dshr_stream_getModelFieldList
  use dshr_stream_mod , only : dshr_stream_getCurrFile
  use dshr_stream_mod , only : dshr_stream_setCurrFile
  use perf_mod       ! timing
  use pio            ! pio

  implicit none
  private

  ! !PUBLIC TYPES:
  public shr_strdata_type

  ! !PUBLIC MEMBER FUNCTIONS:
  public ::  dshr_strdata_init
  public ::  dshr_strdata_readnml
  public ::  dshr_strdata_setOrbs
  public ::  dshr_strdata_advance
  public ::  dshr_strdata_restRead
  public ::  dshr_strdata_restWrite

  private :: dshr_strdata_print
  private :: dshr_strdata_readLBUB
  private :: dshr_strdata_readstrm
  private :: dshr_create_mesh_from_grid
  !private :: dshr_strdata_readstrm_fullfile ! TODO: implememnt this

  ! !PRIVATE:

  ! constants
  integer                              :: debug    = 0  ! local debug flag
  character(len=*) ,parameter, public  :: shr_strdata_nullstr = 'null'
  character(len=*) ,parameter          :: shr_strdata_unset = 'NOT_SET'
  real(R8)         ,parameter, private :: dtlimit_default = 1.5_R8

  ! mapping types
  integer , public, parameter :: nmappers       = 8
  integer , public, parameter :: mapunset       = 0
  integer , public, parameter :: mapbilnr       = 1
  integer , public, parameter :: mapconsf       = 2
  integer , public, parameter :: mapconsd       = 3
  integer , public, parameter :: mappatch       = 4
  integer , public, parameter :: mapfcopy       = 5
  integer , public, parameter :: mapnstod       = 6 ! nearest source to destination
  integer , public, parameter :: mapnstod_consd = 7 ! nearest source to destination followed by conservative dst
  integer , public, parameter :: mapnstod_consf = 8 ! nearest source to destination followed by conservative frac
  character(len=*) , public, parameter :: mapnames(nmappers) = & 
       (/'bilinear   ','consf      ','consd      ',&
         'patch      ','fcopy      ','nstod      ',&
         'nstod_consd','nstod_consf'/)

  ! stream maximum sizes
  integer ,parameter          :: nStrMax = 30
  integer ,parameter          :: nVecMax = 30

  type shr_strdata_type
     ! namelist dependent input (shr_strdata_nml)
     integer                        :: nstreams                    ! actual number of streams (set in 
     character(CL)                  :: dataMode                    ! flags physics options wrt input data
     character(CL)                  :: streams (nStrMax)           ! stream description file names
     character(CL)                  :: tintalgo(nStrMax)           ! time interpolation algorithm
     character(CL)                  :: taxMode (nStrMax)           ! time axis cycling mode
     real(R8)                       :: dtlimit (nStrMax)           ! dt max/min limit
     character(CL)                  :: vectors (nVecMax)           ! define vectors to vector map
     character(CL)                  :: mapalgo (nStrMax)           ! scalar map algorithm
     character(CL)                  :: readmode(nStrMax)           ! file read mode
     type(shr_stream_streamType)    :: stream(nStrMax)             ! stream info
     integer                        :: io_type                     ! pio setting
     integer                        :: io_format                   ! pio setting
     type(iosystem_desc_t), pointer :: pio_subsystem => null()     ! pio setting
     type(io_desc_t)                :: pio_iodesc(nStrMax)         ! pio setting
     type(ESMF_Mesh)                :: mesh_model                  ! model mesh
     type(ESMF_Mesh)                :: mesh_stream(nStrnMax)       ! mesh for each stream
     type(ESMF_FieldBundle)         :: FB_stream_alltimes(nstrMax) ! field bundle for stream n for all time slices for stream
     type(ESMF_FieldBundle)         :: FB_stream_lbound(nstrMax)   ! field bundle for stream n for lb of time period (on stream grid)
     type(ESMF_FieldBundle)         :: FB_stream_ubound(nstrMax)   ! field bundle for stream n for lb of time period (on stream grid)
     type(ESMF_FieldBundle)         :: FB_model_lbound(nstrMax)    ! field bundle for stream n for lb of time period (on model grid)
     type(ESMF_FieldBundle)         :: FB_model_ubound(nstrMax)    ! field bundle for stream n for lb of time period (on model grid)
     type(ESMF_FieldBundle)         :: FB_model(nStrMax)           ! field bundle for stream n at model time and on model grid
     type(ESMF_RouteHandle)         :: RH_stream2model(nstrMax)    ! stream n -> model mesh mapping
     integer                        :: nvectors                    ! number of vectors
     integer                        :: ustrm (nVecMax)             ! vector interpolation stream  n -> model      
     integer                        :: vstrm (nVecMax)             ! vector interpolation strearm n -> model
     type(ESMF_Field)               :: field_coszen(nStrMax)       ! needed for coszen time interp 
     integer                        :: ymdLB(nStrMax)              ! stream n ymd lower bound
     integer                        :: todLB(nStrMax)              ! stream n time of day lower bound
     integer                        :: ymdUB(nStrMax)              ! stream n ymd upper bound
     integer                        :: todUB(nStrMax)              ! stream n time of day upper bound
     integer                        :: ymd                         ! model ymd
     integer                        :: tod                         ! model tod
     real(R8)                       :: dtmin(nStrMax)
     real(R8)                       :: dtmax(nStrMax)
     character(CL)                  :: calendar                    ! model calendar for ymd,tod
     real(R8)                       :: eccen                       ! required by stream  cosz t-interp method
     real(R8)                       :: mvelpp                      ! required by stream  cosz t-interp method
     real(R8)                       :: lambm0                      ! required by stream  cosz t-interp method
     real(R8)                       :: obliqr                      ! required by stream  cosz t-interp method
     integer                        :: modeldt                     ! model dt in seconds
     character(CL)                  :: allocstring = 'strdata_allocated'
  end type shr_strdata_type

  real(R8),parameter,private :: deg2rad = SHR_CONST_PI/180.0_R8

!===============================================================================
contains
!===============================================================================

  subroutine shr_strdata_create(sdat, name, mpicom, compid, &
       yearFirst, yearLast, yearAlign, offset, &
       mesh_filename, stream_filePath, stream_filenames, stream_fieldames, model_fieldnames, 
       stream_readmode, taxMode, dtlimit, tintalgo, mapalgo, model_calendar, rc) 

       ! Set strdata and stream info from fortran interface.
       ! Note: When this is called, previous settings are reset to defaults
       !       and then the values passed are used.

    ! input/output arguments
    type(shr_strdata_type),intent(inout):: sdat                 ! strdata data data-type
    character(*)          ,intent(in)   :: name                 ! name of strdata
    integer(IN)           ,intent(in)   :: mpicom               ! mpi comm
    integer(IN)           ,intent(in)   :: compid               ! component id (needed for pio)
    integer(IN)           ,intent(in)   :: yearFirst            ! first year to use
    integer(IN)           ,intent(in)   :: yearLast             ! last  year to use
    integer(IN)           ,intent(in)   :: yearAlign            ! align yearFirst with this model year
    integer(IN)           ,intent(in)   :: offset               ! offset in seconds of stream data
    character(*)          ,intent(in)   :: mesh_filename        ! mesh filename 
    character(*)          ,intent(in)   :: stream_filepath,     ! path to stream files
    character(*)          ,intent(in)   :: stream_filenames(:)  ! stream filenames in stream_filepath
    character(*)          ,intent(in)   :: stream_fieldnames(:) ! stream field names
    character(*)          ,intent(in)   :: model_fieldnames(:)  ! model field names corresponding to stream_fieldnames
    character(*),optional ,intent(in)   :: stream_readmode      ! file read mode (time slice or all times
    character(*),optional ,intent(in)   :: tintalgo             ! time interpolation algorithm
    character(*),optional ,intent(in)   :: mapalgo              ! scalar map algorithm
    character(*),optional ,intent(in)   :: taxMode
    real(R8)    ,optional ,intent(in)   :: dtlimit
    character(*),optional, intent(in)   :: calendar

    ! local variables
    character(*),parameter :: subName = "(shr_strdata_create) "
    character(*),parameter ::   F00 = "('(shr_strdata_create) ',8a)"
    !-------------------------------------------------------------------------------

    ! The following only sets defaults - but does not read the namelist
    sdat%nstreams = 1
    call shr_strdata_pioinit(sdat, compid)

    call shr_strdata_readnml(sdat, mpicom=mpicom)
    if (present(taxMode)) then
       sdat%taxMode(1) = taxMode
       if (trim(sdat%taxMode(1)) == trim(shr_stream_taxis_extend)) sdat%dtlimit(1) = 1.0e30
    endif
    if (present(dtlimit)) sdat%dtlimit(1) = dtlimit
    if (present(mapalgo)) sdat%mapalgo(1) = mapalgo
    if (present(tintalgo))sdat%tintalgo(1) = tintalgo
    if (present(calendar))sdat%calendar = trim(shr_cal_calendarName(trim(calendar)))
    endif

    !---- Backwards compatibility requires Z be optional

    if (present(domZvarName)) then
       call shr_stream_set(sdat%stream(1),yearFirst,yearLast,yearAlign,offset,taxMode, &
            fldListFile,fldListModel,domFilePath,domFileName, &
            domTvarName,domXvarName,domYvarName,domZvarName, &
            domAreaName,domMaskName, &
            filePath,filename,trim(name))
    else
       call shr_stream_set(sdat%stream(1),yearFirst,yearLast,yearAlign,offset,taxMode, &
            fldListFile,fldListModel,domFilePath,domFileName, &
            domTvarName,domXvarName,domYvarName,'undefined', &
            domAreaName,domMaskName, &
            filePath,filename,trim(name))
    endif

    if (present(nzg)) then
       call shr_strdata_init(sdat, mpicom, compid, gsmap=gsmap, ggrid=ggrid, nxg=nxg, nyg=nyg, nzg=nzg)
    else
       call shr_strdata_init(sdat, mpicom, compid, gsmap=gsmap, ggrid=ggrid, nxg=nxg, nyg=nyg, nzg=1)
    endif

  end subroutine shr_strdata_create_newway

  subroutine dshr_strdata_readnml(sdat, file,  mpicom, logunit, rc)

    ! This is ONLY relevant to data models caps - and not to calls from component models that
    ! leverage stream data functionality

    ! input/output arguments
    type(shr_strdata_type) ,intent(inout) :: sdat    ! strdata data data-type
    character(*)           ,intent(in)    :: file    ! file to read strdata from
    integer                ,intent(in)    :: mpicom  ! mpi comm
    integer                ,intent(in)    :: logunit ! output logunit
    integer, optional      ,intent(out)   :: rc      ! return code

    ! local variables
    integer        :: rCode         ! return code
    integer        :: nUnit         ! fortran i/o unit number
    integer        :: n             ! generic loop index
    integer        :: my_task       ! my task number, 0 is default
    integer        :: master_task   ! master task number, 0 is default
    integer        :: ntasks        ! total number of tasks

    !----- temporary/local namelist vars to read int -----
    character(CL) :: dataMode          ! flags physics options wrt input data
    character(CL) :: domainFile        ! model mesh file
    character(CL) :: streams(nStrMax)  ! stream description file names
    character(CL) :: taxMode(nStrMax)  ! time axis cycling mode
    real(R8)      :: dtlimit(nStrMax)  ! delta time limiter
    character(CL) :: vectors(nVecMax)  ! define vectors to vector map
    character(CL) :: mapalgo(nStrMax)  ! scalar map algorithm
    character(CL) :: tintalgo(nStrMax) ! time interpolation algorithm
    character(CL) :: readmode(nStrMax) ! file read mode (single time slice or full)
    character(CL) :: fileName          ! generic file name
    integer       :: yearFirst         ! first year to use in data stream
    integer       :: yearLast          ! last  year to use in data stream
    integer       :: yearAlign         ! data year that aligns with yearFirst

    !----- define namelist -----
    namelist / shr_strdata_nml / &
           dataMode        &
         , domainFile      &
         , streams         & ! string has txtfile, yearFirst, yearLast, yearAlign
         , taxMode         &
         , dtlimit         &
         , vectors         &
         , mapalgo         &
         , tintalgo        &
         , readmode

    !----- formats -----
    character(*),parameter :: subName = "(shr_strdata_readnml) "
    character(*),parameter ::   F00 = "('(shr_strdata_readnml) ',8a)"
    character(*),parameter ::   F01 = "('(shr_strdata_readnml) ',a,i6,a)"
    character(*),parameter ::   F02 = "('(shr_strdata_readnml) ',a,es13.6)"
    character(*),parameter ::   F03 = "('(shr_strdata_readnml) ',a,l6)"
    character(*),parameter ::   F04 = "('(shr_strdata_readnml) ',a,i2,a,a)"
    character(*),parameter ::   F20 = "('(shr_strdata_readnml) ',a,i6,a)"
    character(*),parameter ::   F90 = "('(shr_strdata_readnml) ',58('-'))"
    !-------------------------------------------------------------------------------

    if (present(rc)) rc = 0

    my_task = 0
    master_task = 0
    ntasks = 1
    if (present(mpicom)) then
       call MPI_COMM_RANK(mpicom,my_task,rCode)
       call MPI_COMM_SIZE(mpicom,ntasks,rCode)
    endif

    !--master--task--
    if (my_task == master_task) then

       !----------------------------------------------------------------------------
       ! set default values for namelist vars
       !----------------------------------------------------------------------------
       dataMode    = 'NULL'
       domainFile  = trim(shr_strdata_nullstr)
       streams(:)  = trim(shr_strdata_nullstr)
       taxMode(:)  = trim(shr_stream_taxis_cycle)
       dtlimit(:)  = dtlimit_default
       vectors(:)  = trim(shr_strdata_nullstr)
       mapalgo(:)  = 'bilinear'
       tintalgo(:) = 'linear'
       readmode(:) = 'single'


       !----------------------------------------------------------------------------
       ! read input namelist
       !----------------------------------------------------------------------------
       if (present(file)) then
          write(logunit,F00) 'reading input namelist file: ',trim(file)
          call shr_sys_flush(logunit)
          open (newunit=nUnit,file=trim(file),status="old",action="read")
          call shr_nl_find_group_name(nUnit, 'shr_strdata_nml', status=rCode)
          if (rCode == 0) then
             read (nUnit, nml=shr_strdata_nml, iostat=rCode)
             if (rCode /= 0) then
                write(logunit,F01) 'ERROR: reading input namelist shr_strdata_input from file, &
                     &'//trim(file)//' iostat=',rCode
                call shr_sys_abort(subName//": namelist read error "//trim(file))
             end if
          end if
          close(nUnit)
       endif

       !----------------------------------------------------------------------------
       ! copy temporary/local namelist vars into data structure
       !----------------------------------------------------------------------------
       sdat%nstreams    = 0
       do n=1,nStrMax
          call dshr_stream_default(sdat%stream(n))
       enddo
       sdat%dataMode    = dataMode
       sdat%domainFile  = domainFile
       sdat%streams(:)  = streams(:)
       sdat%taxMode(:)  = taxMode(:)
       sdat%dtlimit(:)  = dtlimit(:)
       sdat%vectors(:)  = vectors(:)
       sdat%mapalgo(:)  = mapalgo(:)
       sdat%tintalgo(:) = tintalgo(:)
       sdat%readmode(:) = readmode(:)
       do n=1,nStrMax
          if (trim(streams(n)) /= trim(shr_strdata_nullstr)) sdat%nstreams = max(sdat%nstreams,n)
          if (trim(sdat%taxMode(n)) == trim(shr_stream_taxis_extend)) sdat%dtlimit(n) = 1.0e30
       end do
       sdat%nvectors = 0
       do n=1,nVecMax
          if (trim(vectors(n)) /= trim(shr_strdata_nullstr)) sdat%nvectors = n
       end do

       do n = 1,sdat%nstreams
          if (trim(sdat%streams(n)) /= shr_strdata_nullstr) then
             ! extract fileName (stream description text file), yearAlign, yearFirst, yearLast from sdat%streams(n)
             call dshr_stream_parseInput(sdat%streams(n), fileName, yearAlign, yearFirst, yearLast)

             ! initialize stream datatype, read description text file
             call dshr_stream_init(sdat%stream(n), fileName, yearFirst, yearLast, yearAlign, trim(sdat%taxMode(n)))
          end if
       enddo

       !   call dshr_strdata_print(sdat,trim(file)//' NML_ONLY')

    endif   ! master_task
    !--master--task--

    if (present(mpicom)) then
       call shr_mpi_bcast(sdat%dataMode  ,mpicom,'dataMode')
       call shr_mpi_bcast(sdat%domainFile,mpicom,'domainFile')
       call shr_mpi_bcast(sdat%calendar  ,mpicom,'calendar')
       call shr_mpi_bcast(sdat%nstreams  ,mpicom,'nstreams')
       call shr_mpi_bcast(sdat%nvectors  ,mpicom,'nvectors')
       call shr_mpi_bcast(sdat%streams   ,mpicom,'streams')
       call shr_mpi_bcast(sdat%taxMode   ,mpicom,'taxMode')
       call shr_mpi_bcast(sdat%dtlimit   ,mpicom,'dtlimit')
       call shr_mpi_bcast(sdat%vectors   ,mpicom,'vectors')
       call shr_mpi_bcast(sdat%mapalgo   ,mpicom,'mapalgo')
       call shr_mpi_bcast(sdat%tintalgo  ,mpicom,'tintalgo')
       call shr_mpi_bcast(sdat%readmode  ,mpicom,'readmode')
    endif

    sdat%ymdLB = -1
    sdat%todLB = -1
    sdat%ymdUB = -1
    sdat%todUB = -1
    sdat%dtmin = 1.0e30
    sdat%dtmax = 0.0
    sdat%nxg   = 0
    sdat%nyg   = 0
    sdat%nzg   = 0
    sdat%eccen  = SHR_ORB_UNDEF_REAL
    sdat%mvelpp = SHR_ORB_UNDEF_REAL
    sdat%lambm0 = SHR_ORB_UNDEF_REAL
    sdat%obliqr = SHR_ORB_UNDEF_REAL
    sdat%modeldt = 0
    sdat%calendar = shr_cal_noleap

  end subroutine dshr_strdata_readnml

  !===============================================================================

  subroutine dshr_strdata_init(gcomp, mesh_file, sdat, mesh_model, &
       reset_domain_mask, compid, mpicom, my_task, logunit, rc) 

    ! Note: for the variable strdata_domain_fracname_from_stream:
    ! This variable is applicable for data models that read the model domain from the
    ! domain of the first stream; it is an error to provide this variable in other
    ! situations. If present, then we read the data model's domain fraction from the
    ! first stream file, and this variable provides the name of the frac field on this
    ! file. If absent, then (if we are taking the model domain from the domain of the
    ! first stream) we do not read a frac field, and instead we set frac to 1 wherever
    ! mask is 1, and set frac to 0 wherever mask is 0.
    ! BUG(wjs, 2018-05-01, ESMCI/cime#2515) Ideally we'd like to get this frac variable
    ! name in a more general and robust way; see comments in that issue for more details.

    use shr_pio_mod, only: shr_pio_getiosys, shr_pio_getiotype, shr_pio_getioformat

    ! input/output variables
    type(ESMF_GridComp)    , intent(inout) :: gcomp 
    type(shr_strdata_type) , intent(inout) :: sdat
    type(ESMF_Mesh)        , intent(in)    :: mesh
    character(CL), optional, intent(in)    :: domain_filename 
    integer                , intent(in)    :: compid
    integer                , intent(in)    :: mpicom
    integer                , intent(in)    :: my_task
    integer                , intent(in)    :: logunit
    integer                , intent(out)   :: rc 

    ! local variables
    integer                     :: n,m,k    ! generic index
    character(CL)               :: filePath ! generic file path
    character(CL)               :: fileName ! generic file name
    character(CS)               :: timeName ! domain file: time variable name
    character(CS)               :: lonName  ! domain file: lon  variable name
    character(CS)               :: latName  ! domain file: lat  variable name
    character(CS)               :: hgtName  ! domain file: hgt  variable name
    character(CS)               :: maskName ! domain file: mask variable name
    character(CS)               :: areaName ! domain file: area variable name
    integer                     :: nfiles
    character(CL)               :: localFn ! name of acquired file
    character(CXX)              :: fldList ! list of fields
    character(CS)               :: uname   ! u vector field name
    character(CS)               :: vname   ! v vector field name
    type(ESMF_Field)            :: lfield_src
    type(ESMF_Field)            :: lfield_dst
    logical                     :: scmMode      ! single column mode
    real(R8)                    :: scmLat       ! single column lat
    real(R8)                    :: scmLon       ! single column lon
    integer          ,pointer   :: dof(:)
    integer          ,parameter :: master_task = 0
    character(*)     ,parameter :: F00 = "('(shr_strdata_init_mapping) ',8a)"
    character(len=*) ,parameter :: subname = "(shr_strdata_init_mapping) "
    character(*)     ,parameter :: F00 = "('(shr_strdata_init) ',8a)"
    integer         , parameter :: master_task = 0
    character(len=*), parameter :: subname = "(shr_strdata_init) "
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! ------------------------------------
    ! Initialize sdat PIO 
    ! ------------------------------------

    call NUOPC_CompAttributeGet(gcomp, name='MCTID', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) compid

    sdat%pio_subsystem => shr_pio_getiosys(compid)
    sdat%io_type       =  shr_pio_getiotype(compid)
    sdat%io_format     =  shr_pio_getioformat(compid)

    ! ------------------------------------
    ! Determine if will use single column
    ! ------------------------------------

    ! TODO: implement single column functionality

    call NUOPC_CompAttributeGet(gcomp, name='scmlon', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) scmlon
    call NUOPC_CompAttributeGet(gcomp, name='scmlat', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) scmlat
    call NUOPC_CompAttributeGet(gcomp, name='single_column', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) scmMode

    ! ------------------------------------
    ! Determine the model mesh
    ! ------------------------------------

    if (trim(mesh_filename) == 'create_mesh') then
       if (.not. present(domain_filename)) call shr_sys_abort('domain filename must be present')
       call dshr_create_mesh_from_grid(trim(domain_filename), sdat%mesh, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if (my_task == master_task) then
          write(logunit,*) trim(subname)// " creating mesh from " // trim(domain_filename)
       end if
    else
       sdat%mesh = ESMF_MeshCreate(trim(mesh_filename), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if (my_task == master_task) then
          write(logunit,*) trim(subname)// " obtaining mesh from " // trim(mesh_filename)
       end if
    end if

    ! ------------------------------------
    ! Update sdat streams info - assumes that readnml has already been called for the shr_strdata namelist
    ! ------------------------------------

    if (my_task == master_task) then

       do n=1,nStrMax
          ! Determine the number of streams
          if (trim(sdat%streams(n)) /= trim(shr_strdata_nullstr)) then
             sdat%nstreams = max(sdat%nstreams,n)
          end if
          ! Determine number of streams - first check if a filename is defined in the stream
          call dshr_stream_getNFiles(sdat%stream(n), nfiles)
          if (nfiles > 0) then
             sdat%nstreams = max(sdat%nstreams,n)
          end if
          if (trim(sdat%taxMode(n)) == trim(shr_stream_taxis_extend)) then
             sdat%dtlimit(n) = 1.0e30
          end if
       end do
       ! Determine the number of vectors to map
       sdat%nvectors = 0
       do n=1,nVecMax
          if (trim(sdat%vectors(n)) /= trim(shr_strdata_nullstr)) then
             sdat%nvectors = n
          end if
       end do
    endif
    call shr_mpi_bcast(sdat%nstreams  ,mpicom,'nstreams')
    call shr_mpi_bcast(sdat%nvectors  ,mpicom,'nvectors')
    call shr_mpi_bcast(sdat%dtlimit   ,mpicom,'dtlimit')

    ! ------------------------------------
    ! Loop through the streams
    ! ------------------------------------

    do n = 1,sdat%nstreams
       
       ! ------------------------------------
       ! Create the stream mesh
       ! ------------------------------------

       if (my_task == master_task) then
          call dshr_stream_getDomainInfo(sdat%stream(n), filePath, fileName, timeName)
          stream_grid_file = trim(filePath)//adjustl(fileName)
       endif

       ! ASSUME THAT WILL READ IN MESH FOR EACH STREAM
       call shr_mpi_bcast(stream_grid_file, mpicom)
       grid_stream = ESMF_GridCreate(stream_grid_file, fileformat=ESMF_FILEFORMAT_SCRIP, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       mesh_stream =  ESMF_MeshCreate(grid_stream, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! ------------------------------------
       ! Initialize the stream pio decomp - use the stream mesh distgrid 
       ! ------------------------------------

       call ESMF_DistGridGet(distgrid, dimCount=dimCount, tileCount=tileCount, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       allocate(minIndexPTile(dimCount, tileCount), maxIndexPTile(dimCount, tileCount))
       call ESMF_DistGridGet(distgrid, minIndexPTile=minIndexPTile, maxIndexPTile=maxIndexPTile, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       lnx = maxval(maxIndexPTile)
       lny = 1

       call ESMF_DistGridGet(distgrid, localDE=0, elementCount=ns, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       allocate(dof(ns))
       call ESMF_DistGridGet(distgrid, localDE=0, seqIndexList=dof, rc=rc)
       write(tmpstr,*) subname,' dof = ',ns, size(dof), dof(1), dof(ns)  !,minval(dof),maxval(dof)
       call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)
       call pio_initdecomp(sdat%pio_subsystem, pio_double, (/lnx, lny/), dof, sdat%pio_iodesc(n))
       deallocate(dof)

       ! ------------------------------------
       ! broadcast the stream calendar
       ! ------------------------------------

       call shr_mpi_bcast(sdat%stream(n)%calendar, mpicom)

       ! ------------------------------------
       ! Create sdat field bundles 
       ! ------------------------------------

       if (my_task == master_task) then
          call dshr_stream_getModelFieldList(sdat%stream(n),fldList)
       endif

       ! Create empty field bundles
       sdat%FB_stream_lbound(n) = ESMF_FieldBundleCreate(rc=rc) ! stream mesh at lower time bound 
       sdat%FB_stream_ubound(n) = ESMF_FieldBundleCreate(rc=rc) ! stream mesh at upper time bound
       sdat%FB_model_lbound(n)  = ESMF_FieldBundleCreate(rc=rc) ! spatial interpolation to model mesh
       sdat%FB_model_lbound(n)  = ESMF_FieldBundleCreate(rc=rc) ! spatial interpolation to model mesh
       sdat%FB_model(n)         = ESMF_FieldBundleCreate(rc=rc) ! time interpolation on model mesh

       ! get number of field names in fldList (colon delimited string)
       nflds = shr_string_listGetNum(fldList)
       ! loop over field names in fldList
       do nfld = 1,nflds
          ! get fldname in colon delimited string
          call shr_string_listGetName(fldlist,nfld,fldname)

          ! create fields on stream mesh at two time levels
          lfield = ESMF_FieldCreate(mesh_stream, ESMF_TYPEKIND_R8, name=trim(fldname), meshloc=ESMF_MESHLOC_ELEMENT, r=rc)
          call ESMF_FieldBundleAdd(sdat%FB_stream_lbound, (/lfield/), rc=rc)
          call ESMF_FieldBundleAdd(sdat%FB_stream_ubound, (/lfield/), rc=rc)
          call ESMF_FieldDestroy(lfield, rc)

          ! create fields on model mesh at two time levels and for the final interpolated time
          lfield = ESMF_FieldCreate(mesh_model, ESMF_TYPEKIND_R8, name=trim(fldname), meshloc=ESMF_MESHLOC_ELEMENT, r=rc)
          call ESMF_FieldBundleAdd(sdat%FB_model_lbound(n), (/lfield/), rc=rc)
          call ESMF_FieldBundleAdd(sdat%FB_model_ubound(n), (/lfield/), rc=rc)
          call ESMF_FieldBundleAdd(sdat%FB_model(n)       , (/lfield/), rc=rc)
          call ESMF_FieldDestroy(lfield, rc)
       end do

       ! Create a field for coszen time interpolation for this stream if needed
       if (trim(sdat%tintalgo(n)) == 'coszen') then
          sdat%field_coszen = ESMF_FieldCreate(mesh_stream, ESMF_TYPEKIND_R8, name='tavCosz', meshloc=ESMF_MESHLOC_ELEMENT, r=rc)
       endif

    end do

    ! ------------------------------------
    ! Create the route handles
    ! ------------------------------------

    ! create the source and destination fields needed for the route handles
    ! these fields will be used to create the route handles
    ! since all fields in a stream share the same mesh and there is only a unique model mesh
    ! can do this outside of a stream loop by just using the first stream index

    call FB_getFieldN(sdat%FB_stream_lbound(1), 1, lfield_src, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_getFieldN(sdat%FB_FB_model(1), 1, lfield_dst, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Loop over all active streams
    do n = 1,sdat%nstreams

       ! Create bilinear route handle for stream(n) -> model mapping
       if (sdat%mapalgo(n) == 'bilinear') then
          call ESMF_FieldRegridStore(lfield_src, lfield_dst, routehandle=is_local%wrap%RH(mapfcopy,n), &
               regridmethod=ESMF_REGRIDMETHOD_BILINEAR, polemethod=ESMF_POLEMETHOD_ALLAVG, &
               extrapMethod=ESMF_EXTRAPMETHOD_NEAREST_STOD, &
              !srcTermProcessing=0, ignoreDegenerate=.true., unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, rc=rc)
               srcTermProcessing=0, ignoreDegenerate=.true., rc=rc)
       end if

    end do ! end of loop over streams

    ! ------------------------------------
    ! check vectors and compute ustrm,vstrm
    ! ------------------------------------

    do m = 1,sdat%nvectors
       if (.not. shr_string_listIsValid(sdat%vectors(m))) then
          write(logunit,*) trim(subname),' vec fldlist invalid m=',m,trim(sdat%vectors(m))
          call shr_sys_abort(subname//': vec fldlist invalid:'//trim(sdat%vectors(m)))
       endif
       if (shr_string_listGetNum(sdat%vectors(m)) /= 2) then
          write(logunit,*) trim(subname),' vec fldlist ne 2 m=',m,trim(sdat%vectors(m))
          call shr_sys_abort(subname//': vec fldlist ne 2:'//trim(sdat%vectors(m)))
       endif
       nu = 0
       nv = 0
       call shr_string_listGetName(sdat%vectors(m),1,uname)
       call shr_string_listGetName(sdat%vectors(m),2,vname)
       do n = 1,sdat%nstreams
          if (FB_fldchk(FB_stream_lbound(n), trim(uname), rc=rc)) nu = n
          if (FB_fldchk(FB_stream_lbound(n), trim(vname), rc=rc)) nv = n
       enddo
       if (nu == 0  .or. nv == 0) then
          write(logunit,*) trim(subname),' vec flds not found  m=',m,trim(sdat%vectors(m))
          call shr_sys_abort(subname//': vec flds not found:'//trim(sdat%vectors(m)))
       else if (nu /= nv) then
          ! TODO: a grid comparison was made for the vector fields before - is this necessary?
          write(logunit,*) trim(subname),' vec fld doms not same m=',m,trim(sdat%vectors(m))
          call shr_sys_abort(subname//': vec fld doms not same:'//trim(sdat%vectors(m)))
       else
          sdat%ustrm(m)= nu
          sdat%vstrm(m)= nv
       end if
    enddo

  end subroutine dshr_strdata_init

  !===============================================================================

  subroutine dshr_strdata_advance(sdat, ymd, tod, mpicom, , logunit, istr, timers)

    ! ------------------------------------------------------- !
    ! tcraig, Oct 11 2010.  Mismatching calendars: 4 cases    !
    ! ------------------------------------------------------- !
    ! ymdmod and todmod are the ymd and tod to time           !
    ! interpolate to.  Generally, these are just the model    !
    ! date and time.  Also, always use the stream calendar    !
    ! for time interpolation for reasons described below.     !
    ! When there is a calendar mismatch, support Feb 29 in a  !
    ! special way as needed to get reasonable values.         !
    ! Note that when Feb 29 needs to be treated specially,    !
    ! a discontinuity will be introduced.  The size of that   !
    ! discontinuity will depend on the time series input data.!
    ! ------------------------------------------------------- !
    ! (0) The stream calendar and model calendar are          !
    ! identical.  Proceed in the standard way.                !
    ! ------------------------------------------------------- !
    ! (1) If the stream is a no leap calendar and the model   !
    ! is gregorian, then time interpolate on the noleap       !
    ! calendar.  Then if the model date is Feb 29, compute    !
    ! stream data for Feb 28 by setting ymdmod and todmod to  !
    ! Feb 28.  This results in duplicate stream data on       !
    ! Feb 28 and Feb 29 and a discontinuity at the start of   !
    ! Feb 29.                                                 !
    ! This could be potentially updated by using the gregorian!
    ! calendar for time interpolation when the input data is  !
    ! relatively infrequent (say greater than daily) with the !
    ! following concerns.
    !   - The forcing will not be reproduced identically on   !
    !     the same day with climatological inputs data        !
    !   - Input data with variable input frequency might      !
    !     behave funny
    !   - An arbitrary discontinuity will be introduced in    !
    !     the time interpolation method based upon the        !
    !     logic chosen to transition from reproducing Feb 28  !
    !     on Feb 29 and interpolating to Feb 29.              !
    !   - The time gradient of data will change by adding a   !
    !     day arbitrarily.
    ! ------------------------------------------------------- !
    ! (2) If the stream is a gregorian calendar and the model !
    ! is a noleap calendar, then just time interpolate on the !
    ! gregorian calendar.  The causes Feb 29 stream data      !
    ! to be skipped and will lead to a discontinuity at the   !
    ! start of March 1.                                       !
    ! ------------------------------------------------------- !
    ! (3) If the calendars mismatch and neither of the three  !
    ! cases above are recognized, then abort.                 !
    ! ------------------------------------------------------- !
    
    ! input/output variables
    type(shr_strdata_type) ,intent(inout)       :: sdat
    integer                ,intent(in)          :: ymd    ! current model date
    integer                ,intent(in)          :: tod    ! current model date
    integer                ,intent(in)          :: mpicom
    integer                ,intent(in)          :: logunit
    character(len=*)       ,intent(in),optional :: istr
    logical                ,intent(in),optional :: timers

    ! local variables
    integer                    :: n,m,i,kf             ! generic index
    integer                    :: my_task,npes
    integer    ,parameter      :: master_task = 0
    logical                    :: mssrmlf
    logical,allocatable        :: newData(:)
    integer                    :: ierr
    integer                    :: nu,nv
    integer                    :: lsize,lsizeR,lsizeF
    integer    ,allocatable    :: ymdmod(:)            ! modified model dates to handle Feb 29
    integer                    :: todmod               ! modified model dates to handle Feb 29
    character(len=32)          :: lstr
    logical                    :: ltimers
    real(R8)                   :: flb,fub              ! factor for lb and ub
    real(R8),pointer           :: xlon(:),ylon(:)
    real(R8),pointer           :: lon_model(:)         ! lon radians
    real(R8),pointer           :: lat_model(:)         ! lat radians
    real(R8),pointer           :: cosz(:)              ! cosz
    real(R8),pointer           :: tavCosz(:)           ! cosz, time avg over [LB,UB]
    real(R8),pointer           :: dataptr(:)           ! pointer into field bundle
    real(R8),pointer           :: dataptr_lbound(:)    ! pointer into field bundle 
    real(R8),parameter         :: solZenMin = 0.001_R8 ! minimum solar zenith angle
    type(ESMF_Time)            :: timeLB, timeUB       ! lb and ub times
    type(ESMF_TimeInterval)    :: timeint              ! delta time
    integer                    :: dday                 ! delta days
    real(R8)                   :: dtime                ! delta time
    integer                    :: uvar,vvar
    character(CS)              :: uname                ! u vector field name
    character(CS)              :: vname                ! v vector field name
    integer                    :: year,month,day       ! date year month day
    integer         ,parameter :: tadj = 2
    character(len=*),parameter :: timname = "_strd_adv"
    character(*)    ,parameter :: subname = "(shr_strdata_advance) "
    !-------------------------------------------------------------------------------

    if (sdat%nstreams < 1) return

    lstr = ''
    if (present(istr)) then
       lstr = trim(istr)
    endif

    ltimers = .true.
    if (present(timers)) then
       ltimers = timers
    endif

    if (.not.ltimers) call t_adj_detailf(tadj)

    call t_barrierf(trim(lstr)//trim(timname)//'_total_BARRIER',mpicom)
    call t_startf(trim(lstr)//trim(timname)//'_total')

    ! use vm here and gcomp
    call MPI_COMM_SIZE(mpicom,npes,ierr)
    call MPI_COMM_RANK(mpicom,my_task,ierr)

    mssrmlf = .false.

    sdat%ymd = ymd
    sdat%tod = tod

    if (sdat%nstreams > 0) then
       allocate(newData(sdat%nstreams))
       allocate(ymdmod(sdat%nstreams))

       do n = 1,sdat%nstreams

          ! ---------------------------------------------------------
          ! Consistency checks 
          ! ---------------------------------------------------------

          ! case(0)
          ymdmod(n) = ymd
          todmod    = tod
          if (trim(sdat%calendar) /= trim(sdat%stream(n)%calendar)) then
             if ( (trim(sdat%calendar) == trim(shr_cal_gregorian)) .and. &
                  (trim(sdat%stream(n)%calendar) == trim(shr_cal_noleap))) then

                ! case (1), set feb 29 = feb 28
                call shr_cal_date2ymd (ymd,year,month,day)
                if (month == 2 .and. day == 29) then
                   call shr_cal_ymd2date(year,2,28,ymdmod(n))
                endif

             else if ((trim(sdat%calendar) == trim(shr_cal_noleap)) .and. &
                      (trim(sdat%stream(n)%calendar) == trim(shr_cal_gregorian))) then
                ! case (2), feb 29 input data will be skipped automatically

             else
                ! case (3), abort
                write(logunit,*) trim(subname),' ERROR: mismatch calendar ', &
                     trim(sdat%calendar),':',trim(sdat%stream(n)%calendar)
                call shr_sys_abort(trim(subname)//' ERROR: mismatch calendar ')
             endif
          endif

          ! ---------------------------------------------------------
          ! Determine if new data is read in - if so
          ! Copy FB_stream_ubound to FB_stream_lbound
          ! Read in new FB_stream_ubound data
          ! ---------------------------------------------------------

          call t_barrierf(trim(lstr)//trim(timname)//'_readLBUB_BARRIER',mpicom)
          call t_startf(trim(lstr)//trim(timname)//'_readLBUB')

          call dshr_strdata_readLBUB(sdat%stream(n),             &
               sdat%pio_subsystem, sdat%io_type, sdat%pio_iodesc(n),   &
               ymdmod(n), todmod, mpicom,                              &
               sdat%FB_stream_lbound(n), sdat%ymdLB(n), sdat%todLB(n), &
               sdat%FB_stream_ubound(n), sdat%ymdUB(n), sdat%todUB(n), &
               sdat%avRFile(n),  trim(sdat%readmode(n)), newData(n),   &
               istr=trim(lstr)//'_readLBUB')

          if (debug > 0) then
             write(logunit,*) trim(subname),' newData flag = ',n,newData(n)
             write(logunit,*) trim(subname),' LB ymd,tod = ',n,sdat%ymdLB(n),sdat%todLB(n)
             write(logunit,*) trim(subname),' UB ymd,tod = ',n,sdat%ymdUB(n),sdat%todUB(n)
          endif

          ! ---------------------------------------------------------
          ! Reset time bounds if newdata read in
          ! ---------------------------------------------------------

          ! If new data was read in - then set sdat%ymdLB, sdat%ymdUB, sdat%dtmin and sdat%dtmax
          if (newData(n)) then
             if (debug > 0) then
                ! diagnostic
                write(logunit,*) trim(subname),' newData RLB = ',n,minval(sdat%FB_stream_lbound(n)%rAttr), &
                     maxval(sdat%FB_stream_lbound(n)%rAttr),sum(sdat%FB_stream_lbound(n)%rAttr)
                write(logunit,*) trim(subname),' newData RUB = ',n,minval(sdat%FB_stream_ubound(n)%rAttr), &
                     maxval(sdat%FB_stream_ubound(n)%rAttr),sum(sdat%FB_stream_ubound(n)%rAttr)
             endif

             call shr_cal_date2ymd(sdat%ymdLB(n),year,month,day)
             call shr_cal_timeSet(timeLB,sdat%ymdLB(n),0,sdat%stream(n)%calendar)
             call shr_cal_timeSet(timeUB,sdat%ymdUB(n),0,sdat%stream(n)%calendar)
             timeint = timeUB-timeLB
             call ESMF_TimeIntervalGet(timeint,StartTimeIn=timeLB,d=dday)
             dtime = abs(real(dday,R8) + real(sdat%todUB(n)-sdat%todLB(n),R8)/shr_const_cDay)

             sdat%dtmin(n) = min(sdat%dtmin(n),dtime)
             sdat%dtmax(n) = max(sdat%dtmax(n),dtime)
             if ((sdat%dtmax(n)/sdat%dtmin(n)) > sdat%dtlimit(n)) then
                write(logunit,*) trim(subname),' ERROR: for stream ',n
                write(logunit,*) trim(subName),' ERROR: dt limit1 ',sdat%dtmax(n),sdat%dtmin(n),sdat%dtlimit(n)
                write(logunit,*) trim(subName),' ERROR: dt limit2 ',sdat%ymdLB(n),sdat%todLB(n), &
                     sdat%ymdUB(n),sdat%todUB(n)
                call shr_sys_abort(trim(subName)//' ERROR dt limit for stream')
             endif
          endif
          call t_stopf(trim(lstr)//trim(timname)//'_readLBUB')

          ! ---------------------------------------------------------
          ! Do spatial interpolation if newdata read in
          ! ---------------------------------------------------------

          ! If new data was read in, interpolate the lower and upper bound data to the model grid
          if (newData(n)) then
             if (sdat%domaps(n)) then
                call t_startf(trim(lstr)//trim(timname)//'_map')
                call FB_fieldRegrid(sdat%FB_stream_lbound, sdat%FB_model_lbound, sdat%rh_stream2model(n), rc=rc) 
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
                call FB_fieldRegrid(sdat%FB_stream_ubound, sdat%FB_model_ubound, sdat%rh_stream2model(n), rc=rc) 
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
                call t_stopf(trim(lstr)//trim(timname)//'_map')
             end if
             if (debug > 0) then
                call FB_diagnose(sdat%FB_model_lb, subname//':FB_model_lb',rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             endif
          endif
       enddo

       ! ---------------------------------------------------------
       ! remap with vectors if needed
       ! ---------------------------------------------------------

       do m = 1,sdat%nvectors
          nu = sdat%ustrm(m)
          nv = sdat%vstrm(m)
          if ((sdat%domaps(nu) .or. sdat%domaps(nv)) .and. (newdata(nu) .or. newdata(nv))) then

             call t_startf(trim(lstr)//trim(timname)//'_vect')
             call shr_string_listGetName(sdat%vectors(m), 1, uname)
             call shr_string_listGetName(sdat%vectors(m), 2, vname)

             ! get source mesh and coordinates
             call ESMF_MeshGet(sdat%mesh_stream(n), spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             allocate(ownedElemCoords_src(spatialDim*numOwnedElements))
             call ESMF_MeshGet(sdat%mesh_stream(n), ownedElemCoords=ownedElemCoords_src)
             if (chkerr(rc,__LINE__,u_FILE_u)) return

             ! get destination mesh and coordinates
             call ESMF_MeshGet(sdat%mesh_model(n), spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             allocate(ownedElemCoords_src(spatialDim*numOwnedElements))
             call ESMF_MeshGet(sdat%mesh_model(n), ownedElemCoords=ownedElemCoords_src)
             if (chkerr(rc,__LINE__,u_FILE_u)) return

             ! create source field and destination fields
             field2d_src = ESMF_FieldCreate(lmesh_src, ESMF_TYPEKIND_R8, name='src2d', &
                  ungriddedLbound=(/1/), ungriddedUbound=(/2/), gridToFieldMap=(/2/), &
                  meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             field2d_dst = ESMF_FieldCreate(lmesh_dst, ESMF_TYPEKIND_R8, name='dst2d', &
                  ungriddedLbound=(/1/), ungriddedUbound=(/2/), gridToFieldMap=(/2/), &
                  meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return

             ! get pointers to source and destination data that will be filled in with rotation to cart3d
             call ESMF_FieldGet(field2d_src, farrayPtr=data2d_src, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call ESMF_FieldGet(field2d_dst, farrayPtr=data2d_dst, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             
             ! map LB: rotate source data, regrid, then rotate back 
             call FB_getFldPtr(sdat%FB_stream_lbound, trim(uname), data_u_src, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call FB_getFldPtr(sdat%FB_stream_lbound, trim(vname), data_v_src, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call FB_getFldPtr(sdat%FB_stream_lbound, trim(uname), data_u_dst, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call FB_getFldPtr(sdat%FB_stream_lbound, trim(vname), data_v_dst, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             do n = 1,size(data_u_src)
                lon = ownedElemCoords_src(2*n-1)
                lat = ownedElemCoords_src(2*n)
                sinlon = sin(lon*deg2rad)
                coslon = cos(lon*deg2rad)
                sinlat = sin(lat*deg2rad)
                coslat = cos(lat*deg2rad)
                data2d_src(1,n) = coslon * data_u_src(n) - sinlon * data_v_src(n)
                data2d_src(2,n) = sinlon * data_u_src(n) + coslon ( data_v_src(n)
             enddo
             call ESMF_FieldRegrid(field2d_src, field2d_dst, RouteHandle, &
                  termorderflag=ESMF_TERMORDER_SRCSEQ, checkflag=checkflag, zeroregion=ESMF_REGION_TOTAL, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             do n = 1,size(data_u_dst)
                lon = ownedElemCoords_dst(2*n-1)
                lat = ownedElemCoords_dst(2*n)
                sinlon = sin(lon*deg2rad)
                coslon = cos(lon*deg2rad)
                sinlat = sin(lat*deg2rad)
                coslat = cos(lat*deg2rad)
                ux = data2d_dst(1,n)
                uy = data2d_dst(2,n)
                uz = data2d_dst(3,n)
                data_u_dst(n) =  coslon * data_u_dst(n) + sinlon * data_v_dst(n)
                data_v_dst(n) = -sinlon * data_u_dst(n) + coslon * data_v_dst(n)
             enddo

             ! map UB: rotate source data, regrid, then rotate back 
             call FB_getFldPtr(sdat%FB_stream_ubound, trim(uname), data_u_src, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call FB_getFldPtr(sdat%FB_stream_ubound, trim(vname), data_v_src, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call FB_getFldPtr(sdat%FB_stream_ubound, trim(uname), data_u_dst, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call FB_getFldPtr(sdat%FB_stream_ubound, trim(vname), data_v_dst, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             do n = 1,size(data_u_src)
                lon = ownedElemCoords_src(2*n-1)
                lat = ownedElemCoords_src(2*n)
                sinlon = sin(lon*deg2rad)
                coslon = cos(lon*deg2rad)
                sinlat = sin(lat*deg2rad)
                coslat = cos(lat*deg2rad)
                data2d_src(1,n) = coslon * data_u_src(n) - sinlon * data_v_src(n)
                data2d_src(2,n) = sinlon * data_u_src(n) + coslon ( data_v_src(n)
             enddo
             call ESMF_FieldRegrid(field2d_src, field2d_dst, RouteHandle, &
                  termorderflag=ESMF_TERMORDER_SRCSEQ, checkflag=checkflag, zeroregion=ESMF_REGION_TOTAL, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             do n = 1,size(data_u_dst)
                lon = ownedElemCoords_dst(2*n-1)
                lat = ownedElemCoords_dst(2*n)
                sinlon = sin(lon*deg2rad)
                coslon = cos(lon*deg2rad)
                sinlat = sin(lat*deg2rad)
                coslat = cos(lat*deg2rad)
                ux = data2d_dst(1,n)
                uy = data2d_dst(2,n)
                uz = data2d_dst(3,n)
                data_u_dst(n) =  coslon * data_u_dst(n) + sinlon * data_v_dst(n)
                data_v_dst(n) = -sinlon * data_u_dst(n) + coslon * data_v_dst(n)
             enddo

             call t_stopf(trim(lstr)//trim(timname)//'_vect')
          endif
       enddo

       ! ---------------------------------------------------------
       ! Do time interpolation to create FB_model
       ! ---------------------------------------------------------

       ! Get field namelist
       call ESMF_FieldBundleGet(FB_model, fieldCount=fieldCount, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       allocate(lfieldNameList(fieldCount))
       call ESMF_FieldBundleGet(FB_model, fieldNameList=lfieldNameList, rc=rc)

       ! get model lat/lon data
       call ESMF_MeshGet(mesh_model, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       allocate(ownedElemCoords(spatialDim*numOwnedElements))
       call ESMF_MeshGet(mesh_model, ownedElemCoords=ownedElemCoords)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       allocate(lon_model(numOwnedElements))
       allocate(lat_model(numOwnedElements))
       do n = 1, numOwnedElements
          lon_model(n) = ownedElemCoords(2*n-1)
          lat_model(n) = ownedElemCoords(2*n)
       end do
       deallocate(ownedElementCoords)

       ! allocate cosz and avgcosz arrays 
       allocate(tavCosz(numOwnedElements),cosz(numOwnedElements))

       do n = 1,sdat%nstreams

          if (trim(sdat%tintalgo(n)) == 'coszen') then

             ! ------------------------------------------
             ! time interpolation method is coszen
             ! ------------------------------------------

             call t_startf(trim(lstr)//trim(timname)//'_coszen')

             ! make sure orb info has been set
             if (sdat%eccen == SHR_ORB_UNDEF_REAL) then
                call shr_sys_abort(subname//' ERROR in orb params for coszen tinterp')
             else if (sdat%modeldt < 1) then
                call shr_sys_abort(subname//' ERROR: model dt < 1 for coszen tinterp')
             endif

             ! get cosz factor
             call t_startf(trim(lstr)//trim(timname)//'_coszenC')
             cosz(:) = 0.0_r8
             call dshr_tInterp_getCosz(cosz, lon_model, lat_model, ymdmod(n), todmod, &
                  sdat%eccen, sdat%mvelpp, sdat%lambm0, sdat%obliqr, sdat%stream(n)%calendar)
             call t_stopf(trim(lstr)//trim(timname)//'_coszenC')

             ! get avg cosz factor
             if (newdata(n)) then
                ! compute a new avg cosz
                call t_startf(trim(lstr)//trim(timname)//'_coszenN')
                call shr_tInterp_getAvgCosz(tavCosz, lon_model, lat_model,  &
                     sdat%ymdLB(n), sdat%todLB(n),  sdat%ymdUB(n), sdat%todUB(n),  &
                     sdat%eccen, sdat%mvelpp, sdat%lambm0, sdat%obliqr, sdat%modeldt, &
                     sdat%stream(n)%calendar)
                sdat%avgCoszen(:) = tavgCosz(:)
                call t_stopf(trim(lstr)//trim(timname)//'_coszenN')
             else
                ! reuse existing avg cosz
                tavgCosz(:) = sdat%avgCoszen(:)
             endif

             ! compute time interperpolate value - LB data normalized with this factor: cosz/tavCosz
             do n = 1,fieldcount
                call FB_getfldptr(sdat%FB_model       , fieldNameList(n), dataptr   , rc=rc)
                call FB_getfldptr(sdat%FB_model_lbound, fieldNameList(n), dataptr_lb, rc=rc)
                do i = 1,lsize
                   if (cosz(i) > solZenMin) then
                      dataptr(i) = dataptr_lb(i)*cosz(i)/sdat%tavCosz(i)
                   else
                      dataptr(i) = 0._r8
                   endif
                end do
             end do

             ! deallocate memory
             call t_stopf(trim(lstr)//trim(timname)//'_coszen')

          elseif (trim(sdat%tintalgo(n)) /= trim(shr_strdata_nullstr)) then

             ! ------------------------------------------
             ! time interpolation method is not coszen
             ! ------------------------------------------

             call t_startf(trim(lstr)//trim(timname)//'_tint')
             call shr_tInterp_getFactors(&
                  sdat%ymdlb(n), sdat%todlb(n), sdat%ymdub(n), sdat%todub(n), &
                  ymdmod(n),todmod,flb,fub, &
                  calendar=sdat%stream(n)%calendar,algo=trim(sdat%tintalgo(n)))
             if (debug > 0) then
                write(logunit,*) trim(subname),' interp = ',n,flb,fub
             endif
             do n = 1,fieldcount
                call FB_getfldptr(sdat%FB_model       , fieldNameList(n), dataptr   , rc=rc)
                call FB_getfldptr(sdat%FB_model_lbound, fieldNameList(n), dataptr_lb, rc=rc)
                call FB_getfldptr(sdat%FB_model_ubound, fieldNameList(n), dataptr_ub, rc=rc)
                dataptr(:) = dataptr_lb(:) * flb + dataptr_ub(:) * fub
             end do
             call t_stopf(trim(lstr)//trim(timname)//'_tint')

          else

             ! ------------------------------------------
             ! zero out stream data for this field
             ! ------------------------------------------

             call t_startf(trim(lstr)//trim(timname)//'_zero')
             do n = 1,fieldcount
                call FB_getfldptr(sdat%FB_model       , fieldNameList(n), dataptr   , rc=rc)
                dataptr(:) = 0._r8
             end do
             call t_stopf(trim(lstr)//trim(timname)//'_zero')

          endif
          deallocate(lfieldNameList(fieldCount))

          if (debug > 0) then
             ! TODO: call field bundle diagnose here
          endif

       enddo
       deallocate(newData)
       deallocate(ymdmod)

    endif    ! nstreams > 0

    deallocate(tavCosz,cosz,lonr,latr)

    call t_stopf(trim(lstr)//trim(timname)//'_total')
    if (.not.ltimers) call t_adj_detailf(-tadj)

  end subroutine dshr_strdata_advance

  !===============================================================================

  subroutine dshr_strdata_restWrite(filename,sdat,mpicom,str1,str2)

    character(len=*)      ,intent(in)    :: filename
    type(shr_strdata_type),intent(inout) :: sdat
    integer               ,intent(in)    :: mpicom
    character(len=*)      ,intent(in)    :: str1
    character(len=*)      ,intent(in)    :: str2

    !--- local ----
    integer     :: my_task,ier

    !----- formats -----
    character(len=*),parameter :: subname = "(shr_strdata_restWrite) "
    !-------------------------------------------------------------------------------

    call MPI_COMM_RANK(mpicom,my_task,ier)
    if (my_task == 0) then
       call dshr_stream_restWrite(sdat%stream,trim(filename),trim(str1),trim(str2),sdat%nstreams)
    endif

  end subroutine dshr_strdata_restWrite

  !===============================================================================

  subroutine dshr_strdata_restRead(filename,sdat,mpicom)

    character(len=*)      ,intent(in)    :: filename
    type(shr_strdata_type),intent(inout) :: sdat
    integer               ,intent(in)    :: mpicom

    !--- local ----
    integer     :: my_task,ier

    !----- formats -----
    character(len=*),parameter :: subname = "(shr_strdata_restRead) "
    !-------------------------------------------------------------------------------

    call MPI_COMM_RANK(mpicom,my_task,ier)
    if (my_task == 0) then
       call dshr_stream_restRead(sdat%stream,trim(filename),sdat%nstreams)
    endif

  end subroutine dshr_strdata_restRead

  !===============================================================================

  subroutine dshr_strdata_setOrbs(sdat,eccen,mvelpp,lambm0,obliqr,modeldt)

    ! input/output variables
    type(shr_strdata_type) ,intent(inout) :: sdat
    real(R8)               ,intent(in)    :: eccen
    real(R8)               ,intent(in)    :: mvelpp
    real(R8)               ,intent(in)    :: lambm0
    real(R8)               ,intent(in)    :: obliqr
    integer                ,intent(in)    :: modeldt
    !-------------------------------------------------------------------------------

    sdat%eccen   = eccen
    sdat%mvelpp  = mvelpp
    sdat%lambm0  = lambm0
    sdat%obliqr  = obliqr
    sdat%modeldt = modeldt

  end subroutine dshr_strdata_setOrbs

  !===============================================================================

  subroutine dshr_strdata_print(sdat, name, logunit)

    !  Print strdata common to all data models

    ! !input/output parameters
    type(shr_strdata_type)    ,intent(in) :: sdat  ! strdata data data-type
    character(len=*),optional ,intent(in) :: name  ! just a name for tracking
    integer                   ,intent(in) :: logunit

    ! local variables
    integer                :: n
    character(CL)          :: lname
    character(*),parameter :: F00 = "('(shr_strdata_print) ',8a)"
    character(*),parameter :: F01 = "('(shr_strdata_print) ',a,i6,a)"
    character(*),parameter :: F02 = "('(shr_strdata_print) ',a,es13.6)"
    character(*),parameter :: F03 = "('(shr_strdata_print) ',a,l6)"
    character(*),parameter :: F04 = "('(shr_strdata_print) ',a,i2,a,a)"
    character(*),parameter :: F05 = "('(shr_strdata_print) ',a,i2,a,i6)"
    character(*),parameter :: F06 = "('(shr_strdata_print) ',a,i2,a,l2)"
    character(*),parameter :: F07 = "('(shr_strdata_print) ',a,i2,a,es13.6)"
    character(*),parameter :: F20 = "('(shr_strdata_print) ',a,i6,a)"
    character(*),parameter :: F90 = "('(shr_strdata_print) ',58('-'))"
    character(*),parameter :: subName = "(shr_strdata_print) "
    !-------------------------------------------------------------------------------

    lname = 'unknown'
    if (present(name)) then
       lname = trim(name)
    endif

    !----------------------------------------------------------------------------
    ! document datatype settings
    !----------------------------------------------------------------------------
    write(logunit,F90)
    write(logunit,F00) "name        = ",trim(lname)
    write(logunit,F00) "dataMode    = ",trim(sdat%dataMode)
    write(logunit,F00) "domainFile  = ",trim(sdat%domainFile)
    write(logunit,F01) "nxg         = ",sdat%nxg
    write(logunit,F01) "nyg         = ",sdat%nyg
    write(logunit,F01) "nzg         = ",sdat%nzg
    write(logunit,F00) "calendar    = ",trim(sdat%calendar)
    write(logunit,F01) "io_type     = ",sdat%io_type
    write(logunit,F02) "eccen       = ",sdat%eccen
    write(logunit,F02) "mvelpp      = ",sdat%mvelpp
    write(logunit,F02) "lambm0      = ",sdat%lambm0
    write(logunit,F02) "obliqr      = ",sdat%obliqr
    write(logunit,F01) "nstreams    = ",sdat%nstreams
    write(logunit,F01) "pio_iotype  = ",sdat%io_type

    do n=1, sdat%nstreams
       write(logunit,F04) "  streams (",n,") = ",trim(sdat%streams(n))
       write(logunit,F04) "  taxMode (",n,") = ",trim(sdat%taxMode(n))
       write(logunit,F07) "  dtlimit (",n,") = ",sdat%dtlimit(n)
       write(logunit,F05) "  strnxg  (",n,") = ",sdat%strnxg(n)
       write(logunit,F05) "  strnyg  (",n,") = ",sdat%strnyg(n)
       write(logunit,F05) "  strnzg  (",n,") = ",sdat%strnzg(n)
       write(logunit,F06) "  domaps  (",n,") = ",sdat%domaps(n)
       write(logunit,F04) "  mapalgo (",n,") = ",trim(sdat%mapalgo(n))
       write(logunit,F04) "  tintalgo(",n,") = ",trim(sdat%tintalgo(n))
       write(logunit,F04) "  readmode(",n,") = ",trim(sdat%readmode(n))
       write(logunit,F01) " "
    end do
    write(logunit,F01) "nvectors    = ",sdat%nvectors
    do n=1, sdat%nvectors
       write(logunit,F04) "  vectors (",n,") = ",trim(sdat%vectors(n))
    end do
    write(logunit,F90)
    call shr_sys_flush(logunit)

  end subroutine dshr_strdata_print

  !===============================================================================

  subroutine dshr_strdata_readLBUB(stream, &
       pio_subsystem, pio_iotype, pio_iodesc, mDate, mSec, mpicom,  &
       avLB, mDateLB, mSecLB, avUB, mDateUB, mSecUB, &
       avFile, readMode, newData, rmOldFile, istr)


    !----- arguments -----
    type(shr_stream_streamType)   ,intent(inout) :: stream
    type(iosystem_desc_t), target ,intent(inout) :: pio_subsystem
    integer                       ,intent(in)    :: pio_iotype
    type(io_desc_t)               ,intent(inout) :: pio_iodesc
    integer                       ,intent(in)    :: mDate  ,mSec
    integer                       ,intent(in)    :: mpicom
    type(ESMF_FieldBundle)        ,intent(inout) :: FB_stream_lbound
    integer                       ,intent(inout) :: mDateLB,mSecLB
    type(ESMF_FieldBundle)        ,intent(inout) :: FB_stream_ubound
    integer                       ,intent(inout) :: mDateUB,mSecUB
    type(ESMF_FieldBundle)        ,intent(inout) :: FB_stream_all(:)
    character(len=*)              ,intent(in)    :: readMode
    logical                       ,intent(out)   :: newData
    logical          ,optional    ,intent(in)    :: rmOldFile
    character(len=*) ,optional    ,intent(in)    :: istr

    ! local variables
    integer           :: my_task, master_task
    integer           :: ierr       ! error code
    integer           :: rCode      ! return code
    logical           :: localCopy,fileexists
    integer           :: ivals(6)   ! bcast buffer
    integer           :: oDateLB,oSecLB,dDateLB
    integer           :: oDateUB,oSecUB,dDateUB
    real(R8)          :: rDateM,rDateLB,rDateUB  ! model,LB,UB dates with fractional days
    integer           :: n_lb, n_ub
    character(CL)     :: fn_lb,fn_ub,fn_next,fn_prev
    character(CL)     :: path
    character(len=32) :: lstr
    real(R8)          :: spd

    character(*), parameter :: subname = '(dshr_strdata_readLBUB) '
    character(*), parameter :: F00   = "('(dshr_strdata_readLBUB) ',8a)"
    character(*), parameter :: F01   = "('(dshr_strdata_readLBUB) ',a,5i8)"

    !-------------------------------------------------------------------------------
    ! PURPOSE:  Read LB and UB stream data
    !----------------------------------------------------------------------------

    lstr = 'dshr_strdata_readLBUB'
    if (present(istr)) then
       lstr = trim(istr)
    endif

    call t_startf(trim(lstr)//'_setup')
    call MPI_COMM_RANK(mpicom,my_task,ierr)
    master_task = 0
    spd = shr_const_cday

    newData = .false.
    n_lb = -1
    n_ub = -1
    fn_lb = 'undefinedlb'
    fn_ub = 'undefinedub'

    oDateLB = mDateLB
    oSecLB  = mSecLB
    oDateUB = mDateUB
    oSecUB  = mSecUB

    rDateM  = real(mDate  ,R8) + real(mSec  ,R8)/spd
    rDateLB = real(mDateLB,R8) + real(mSecLB,R8)/spd
    rDateUB = real(mDateUB,R8) + real(mSecUB,R8)/spd
    call t_stopf(trim(lstr)//'_setup')

    if (rDateM < rDateLB .or. rDateM > rDateUB) then
       call t_startf(trim(lstr)//'_fbound')
       if (my_task == master_task) then
          call dshr_stream_findBounds(stream,mDate,mSec,                 &
               ivals(1),dDateLB,ivals(2),ivals(5),fn_lb, &
               ivals(3),dDateUB,ivals(4),ivals(6),fn_ub  )
          call dshr_stream_getFilePath(stream,path)
       endif
       call t_stopf(trim(lstr)//'_fbound')
       call t_startf(trim(lstr)//'_bcast')
       call shr_mpi_bcast(stream%calendar,mpicom)
       call shr_mpi_bcast(ivals,mpicom)
       mDateLB = ivals(1)
       mSecLB  = ivals(2)
       mDateUB = ivals(3)
       mSecUB  = ivals(4)
       n_lb    = ivals(5)
       n_ub    = ivals(6)
       call t_stopf(trim(lstr)//'_bcast')
    endif

    if (mDateLB /= oDateLB .or. mSecLB /= oSecLB) then
       newdata = .true.

       if (mDateLB == oDateUB .and. mSecLB == oSecUB) then

          ! copy FB_stream_ubound to FB_stream_lbound
          call t_startf(trim(lstr)//'_LB_copy')
          call ESMF_FieldBundleGet(FB_stream_ubound, fieldCount=fieldCount, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          allocate(lfieldNameList(fieldCount))
          call ESMF_FieldBundleGet(FB_stream_ubound, itemNameList=lfieldNameList, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
.          do n = 1,fieldCount
             fldname = trim(lfieldnamelist(n))
             call ESMF_FieldBundleGet(FB_stream_ubound, fieldName=fldname, field=lfield_ub, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call ESMF_FieldBundleGet(FB_stream_lbound, fieldName=fldname, field=lfield_ub, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call ESMF_FieldGet(lfield_ub, farrayPtr=dataptr_ub, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call ESMF_FieldGet(lfield_lb, farrayPtr=dataptr_lb, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             dataptr_lb(:) = dataptr_ub(:)
          end do
          call t_stopf(trim(lstr)//'_LB_copy')

       else

          select case(readMode)
          case ('single')
             call dshr_strdata_readstrm(stream, FB_stream, &
                  pio_subsystem, pio_iotype, pio_iodesc, mpicom, path, &
                  fn_lb, n_lb, istr=trim(lstr)//'_UB', boundstr = 'ub')
          case ('full_file')
             ! TODO: implement full file read
          case default
             write(logunit,F00) "ERROR: Unsupported readmode : ", trim(readMode)
             call shr_sys_abort(subName//"ERROR: Unsupported readmode: "//trim(readMode))
          end select
       endif
    endif

    if (mDateUB /= oDateUB .or. mSecUB /= oSecUB) then
       newdata = .true.

       select case(readMode)
       case ('single')
          call dshr_strdata_readstrm(stream, pio_subsystem, pio_iotype, pio_iodesc, gsMap, avUB, mpicom, &
               path, fn_ub, n_ub,istr=trim(lstr)//'_UB', boundstr = 'ub')
       case ('full_file')
          call dshr_strdata_readstrm_fullfile(stream, pio_subsystem, pio_iotype, &
               gsMap, avUB, avFile, mpicom, &
               path, fn_ub, n_ub,istr=trim(lstr)//'_UB', boundstr = 'ub')
       case default
          write(logunit,F00) "ERROR: Unsupported readmode : ", trim(readMode)
          call shr_sys_abort(subName//"ERROR: Unsupported readmode: "//trim(readMode))
       end select

    endif

    call t_startf(trim(lstr)//'_filemgt')
    !--- determine previous & next data files in list of files ---
    if (my_task == master_task .and. newdata) then
       call dshr_stream_getFilePath(stream,path)
    endif
    call t_stopf(trim(lstr)//'_filemgt')

  end subroutine dshr_strdata_readLBUB

  !===============================================================================

  subroutine dshr_strdata_readstrm(stream, FB_stream, &
       pio_subsystem, pio_iotype, pio_iodesc, mpicom, path, fn, nt, istr, boundstr)

    ! input/output variables
    type(shr_stream_streamType) ,intent(inout)         :: stream
    type(ESMF_FieldBundle)      ,intent(inout)         :: FB_stream
    type(iosystem_desc_t)       ,intent(inout), target :: pio_subsystem
    integer                     ,intent(in)            :: pio_iotype
    type(io_desc_t)             ,intent(inout)         :: pio_iodesc
    integer                     ,intent(in)            :: mpicom
    character(len=*)            ,intent(in)            :: path
    character(len=*)            ,intent(in)            :: fn
    integer                     ,intent(in)            :: nt
    character(len=*),optional   ,intent(in)            :: istr
    character(len=*),optional   ,intent(in)            :: boundstr

    ! local variables
    type(ESMF_Field)              :: lfield
    character(CS)                 :: fldname
    character(CL)                 :: fileName
    character(CL)                 :: currfile
    type(file_desc_t)             :: pioid
    type(var_desc_t)              :: varid
    integer(kind=pio_offset_kind) :: frame
    integer                       :: k, n
    integer                       :: my_task
    integer                       :: master_task
    integer                       :: ierr
    logical                       :: fileexists
    integer                       :: rCode      ! return code
    character(len=32)             :: lstr
    logical                       :: fileopen
    character(ESMF_MAXSTR), allocatable :: lfieldNameList(:)
    character(*), parameter :: subname = '(dshr_strdata_readstrm) '
    character(*), parameter :: F00   = "('(dshr_strdata_readstrm) ',8a)"
    character(*), parameter :: F02   = "('(dshr_strdata_readstrm) ',2a,i8)"
    !-------------------------------------------------------------------------------

    lstr = 'shr_strdata_readstrm'
    if (present(istr)) then
       lstr = trim(istr)
    endif

    bstr = ''
    if (present(boundstr)) then
       bstr = trim(boundstr)
    endif

    call MPI_COMM_RANK(mpicom,my_task,ierr)
    master_task = 0

    ! Set up file to read from
    call t_barrierf(trim(lstr)//'_BARRIER',mpicom)
    call t_startf(trim(lstr)//'_setup')
    if (my_task == master_task) then
       fileName = trim(path)//fn
       inquire(file=trim(fileName),exist=fileExists)
       if (.not. fileExists) then
          write(logunit,F00) "ERROR: file does not exist: ", trim(fileName)
          call shr_sys_abort(subName//"ERROR: file does not exist: "//trim(fileName))
       end if
    endif
    call shr_mpi_bcast(filename,mpicom,'filename')
    call t_stopf(trim(lstr)//'_setup')

    ! Get current file and determine if it is open
    call dshr_stream_getCurrFile(stream, fileopen=fileopen, currfile=currfile, currpioid=pioid)
    if (fileopen .and. currfile==filename) then
       ! don't reopen file, all good
    else
       ! otherwise close the old file if open and open new file
       if (fileopen) then
          if (my_task == master_task) then
             write(logunit,F00) 'close  : ',trim(currfile)
             call shr_sys_flush(logunit)
          endif
          call pio_closefile(pioid)
       endif
       if (my_task == master_task) then
          write(logunit,F00) 'open   : ',trim(filename)
          call shr_sys_flush(logunit)
       endif
       rcode = pio_openfile(pio_subsystem, pioid, pio_iotype, trim(filename), pio_nowrite)
       call dshr_stream_setCurrFile(stream, fileopen=.true., currfile=trim(filename), currpioid=pioid)
    endif
    if (my_task == master_task) then
       write(logunit,*) 'file '// trim(filename),nt
    endif

    ! Always use pio to read in stream data
    call t_startf(trim(lstr)//'_readpio')
    if (my_task == master_task) then
       write(logunit,F02) 'file ' // trim(bstr) //': ',trim(filename), nt
    endif
    call pio_seterrorhandling(pioid,PIO_INTERNAL_ERROR)

    ! Loop over all stream fields in FB_stream
    call ESMF_FieldBundleGet(FB_stream, fieldCount=fieldCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(lfieldNameList(fieldCount))
    call ESMF_FieldBundleGet(FB_stream, itemNameList=lfieldNameList, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    do n = 1,fieldCount
       fldname = trim(lfieldnamelist(n))
       ! get varid of field n
       rcode = pio_inq_varid(pioid, fldname, varid)
       ! set frame to time index
       frame = nt
       call pio_setframe(pioid,varid,frame)
       ! set pointer rdata to field data
       call ESMF_FieldBundleGet(FB_stream, fieldName=fldname, field=lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(field, farrayPtr=dataptr, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       ! read dataptr - which sets FB_stream value for stream field fldname
       call pio_read_darray(pioid, varid, pio_iodesc, dataptr, rcode)
    enddo
    call t_stopf(trim(lstr)//'_readpio')

  end subroutine dshr_strdata_readstrm

!===============================================================================

  subroutine dshr_create_mesh_from_grid(filename, mesh, rc)

    use netcdf , only : nf90_open, nf90_nowrite, nf90_noerr, nf90_close, nf90_strerror
    use netcdf , only : nf90_inq_dimid, nf90_inq_varid, nf90_get_var
    use netcdf , only : nf90_inquire_dimension, nf90_inquire_variable

    ! input/output variables
    character(len=*), intent(in)  :: filename
    type(ESMF_Mesh) , intent(out) :: mesh
    integer         , intent(out) :: rc

    ! local variables
    integer               :: ncid, ierr
    integer               :: dimid_ni, dimid_nj, dimid_nv
    integer               :: ni, nj, nv
    integer               :: varid_xv, varid_yv
    integer               :: maxIndex(2)
    real(r8)              :: mincornerCoord(2)
    real(r8)              :: maxcornerCoord(2)
    type(ESMF_Grid)       :: lgrid
    real(r8), allocatable :: xv(:,:,:), yv(:,:,:)
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! open file
    ierr = nf90_open(filename, NF90_NOWRITE, ncid)
    call nc_check_err(ierr, 'nf90_open', trim(filename))

    ! get dimension ids
    ierr = nf90_inq_dimid(ncid, 'ni', dimid_ni)
    call nc_check_err(ierr, 'nf90_inq_dimid for ni', trim(filename))
    ierr = nf90_inq_dimid(ncid, 'nj', dimid_nj)
    call nc_check_err(ierr, 'nf90_inq_dimid for nj', trim(filename))
    ierr = nf90_inq_dimid(ncid, 'nv', dimid_nv)
    call nc_check_err(ierr, 'nf90_inq_dimid for nv', trim(filename))

    ! get dimension values
    ierr = nf90_inquire_dimension(ncid, dimid_ni, len=ni)
    call nc_check_err(ierr, 'nf90_inq_dimension for ni', trim(filename))
    ierr = nf90_inquire_dimension(ncid, dimid_nj, len=nj)
    call nc_check_err(ierr, 'nf90_inq_dimension for nj', trim(filename))
    ierr = nf90_inquire_dimension(ncid, dimid_nv, len=nv)
    call nc_check_err(ierr, 'nf90_inq_dimension for nv', trim(filename))

    ! get variable ids
    ierr = nf90_inq_varid(ncid, 'xv', varid_xv)
    call nc_check_err(ierr, 'nf90_inq_varid for xv', trim(filename))
    ierr = nf90_inq_varid(ncid, 'yv', varid_yv)
    call nc_check_err(ierr, 'nf90_inq_varid for yv', trim(filename))

    ! allocate memory for variables and get variable values
    allocate(xv(nv,ni,nj), yv(nv,ni,nj))
    ierr = nf90_get_var(ncid, varid_xv, xv)
    call nc_check_err(ierr, 'nf90_get_var for xv', trim(filename))
    ierr = nf90_get_var(ncid, varid_yv, yv)
    call nc_check_err(ierr, 'nf90_get_var for yv', trim(filename))

    ! close file
    ierr = nf90_close(ncid)
    call nc_check_err(ierr, 'nf90_close', trim(filename))

    ! create the grid
    maxIndex(1)       = ni          ! number of lons
    maxIndex(2)       = nj          ! number of lats
    mincornerCoord(1) = xv(1,1,1)   ! min lon
    mincornerCoord(2) = yv(1,1,1)   ! min lat
    maxcornerCoord(1) = xv(3,ni,nj) ! max lon
    maxcornerCoord(2) = yv(3,ni,nj) ! max lat
    deallocate(xv,yv)
    lgrid = ESMF_GridCreateNoPeriDimUfrm (maxindex=maxindex, &
         mincornercoord=mincornercoord, maxcornercoord= maxcornercoord, &
         staggerloclist=(/ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! create the mesh from the grid
    mesh =  ESMF_MeshCreate(lgrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  contains

    subroutine nc_check_err(ierror, description, filename)
      integer     , intent(in) :: ierror
      character(*), intent(in) :: description
      character(*), intent(in) :: filename

      if (ierror /= nf90_noerr) then
         write (*,'(6a)') 'ERROR ', trim(description),'. NetCDF file : "', trim(filename),&
              '". Error message:', trim(nf90_strerror(ierror))
         call shr_sys_abort()
      endif
    end subroutine nc_check_err

  end subroutine dshr_create_mesh_from_grid

end module dshr_strdata_mod

program prepare_gadget
  
  ! This routine reads in the initial particles generated from 2LPT simulations (Manera et al. 2012, Scoccimarro et al. 2012)
  ! and then rewrites then in the correct form for the FoF algorithm to read in. Although this can be used it is not particularly
  ! general and is included mainly as an example routine of how to put the data in the correct format to be read in. 
  !
  ! The exact format that the initial conditions need to be in was chosen as the most general and simplest form and this utility 
  ! is provided such that the main FoF algorithm does not need to be modified. Hence the usage of this routine is strongly advised. 
  ! However, should you decide to read the data straight into the FoF algorithm as it is, not using this routine to reformat it first, 
  ! the lines 94-101 in fof_data.F90 are the lines you will want to look at. 

  character(len=5) :: ninput_str
  character(len=1000) :: inputfile, outputfile

  type io_gheader
    integer :: npart(0:5)  !should be unsigned integer but not supported by fortran90
    real*8  :: mass(0:5)
    real*8  :: time
    real*8  :: redshift
    integer :: flag_sfr
    integer :: flag_feedback
    integer :: npartTotal(0:5) !should be unsigned integer but not supported by fortran90
    integer :: flag_cooling
    integer :: num_files
    real*8  :: BoxSize   !in Gadget code units
    real*8  :: Omega0
    real*8  :: OmegaL0
    real*8  :: HubbleParam
    integer :: flag_stellarage
    integer :: flag_metals
    integer :: hashtabsize
    character :: fill*84   !to fill the header
  end type io_gheader

  character :: red_str1, red_str2
  character(len=5) :: file_str
  character(len=1000) :: filename
  integer*4 :: i, j, k, ninput_int
  integer*4, allocatable, dimension(:) :: id
  real*4, allocatable, dimension(:,:) :: r, v
  type(io_gheader) :: gheader

  ! Reads the command line arguments
  if (command_argument_count().eq.2) then
    call get_command_argument(1,inputfile)
    call get_command_argument(2,ninput_str)
  else 
    write(*,*) 'Usage: prepare_gadget, [Input_Filename, Number_of_Input_Files]'
    STOP
  end if

  read(ninput_str,'(i4.4)') ninput_int

  do i = 0, ninput_int-1
    write(file_str,'(i5)') i
    filename=trim(adjustl(inputfile))//'.'//trim(adjustl(file_str))
    open(unit=11,file=filename,form='unformatted')
    print*, 'Opened: '//trim(adjustl(filename))
    read(11) (gheader%npart(j),j=0,5),(gheader%mass(j),j=0,5),gheader%time,gheader%redshift, &
       & gheader%flag_sfr,gheader%flag_feedback,(gheader%npartTotal(j),j=0,5),gheader%flag_cooling,& 
       & gheader%num_files,gheader%BoxSize,gheader%Omega0,gheader%OmegaL0,gheader%HubbleParam, &
       & gheader%flag_stellarage,gheader%flag_metals,gheader%hashtabsize,gheader%fill
    print*, gheader%redshift
    allocate(r(3,gheader%npart(1)))
    allocate(v(3,gheader%npart(1)))
    allocate(id(gheader%npart(1)))
    read(11)((r(j,k),j=1,3),k=1,gheader%npart(1))
    read(11)((v(j,k),j=1,3),k=1,gheader%npart(1))
    read(11)(id(j),j=1,gheader%npart(1)) 
    write(red_str1,'(i1)') int(gheader%redshift)
    write(red_str2,'(i1)') nint(10.0*(gheader%redshift-int(gheader%redshift)))
    close(11)
    
    filename=trim(adjustl(inputfile))//'_z'//trim(adjustl(red_str1))//'p' & 
            &   //trim(adjustl(red_str2))//'_reform.'//trim(adjustl(file_str))
    open(10,file=filename,form='formatted')
    do j = 1, gheader%npart(1)
      write(10,'(I15, 6F15.6)') id(j), r(1,j), r(2,j), r(3,j), v(1,j), v(2,j), v(3,j)
    end do
    deallocate(id, r, v)
    close(10)
  end do

end program prepare_gadget

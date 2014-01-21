program prepare_lightcone
  
  ! This routine reads in the unformatted binary output from PICOLA's lightcone mode and converts it to ASCII

  character(len=5) :: ninput_str
  character(len=1000) :: inputfile, outputfile

  character(len=5) :: file_str
  character(len=1000) :: foutname, finname
  logical*1 :: errorflag
  integer*4 :: i, j, EOF, nchunk, ninput_int
  real*4, allocatable, dimension(:,:) :: r, v

  ! Reads the command line arguments
  if (command_argument_count().eq.2) then
    call get_command_argument(1,inputfile)
    call get_command_argument(2,ninput_str)
  else 
    write(*,*) 'Usage: prepare_lightcone, [Input_Filename, Number_of_Input_Files]'
    STOP
  end if

  read(ninput_str,'(i4.4)') ninput_int

  do i = 0, ninput_int-1

    write(file_str,'(i5)') i

    foutname=trim(adjustl(inputfile))//'_reform.'//trim(adjustl(file_str))
    open(10,file=foutname,form='formatted')

    finname=trim(adjustl(inputfile))//'.'//trim(adjustl(file_str))
    open(unit=11,file=finname,form='unformatted')
    print*, 'Opened: '//trim(adjustl(finname))

    errorflag = .false.
    do 
      read(11, IOSTAT=EOF) nchunk
      print*, nchunk
      if (EOF .gt. 0) then
        print*, 'Read error in file: ', finname
        print*, 'Exiting program' 
        errorflag = .true.
      else if (EOF .lt. 0) then
        exit
      else
        allocate(r(3,nchunk))
        allocate(v(3,nchunk))
        read(11)((r(1,j),r(2,j),r(3,j),v(1,j),v(2,j),v(3,j)),j=1,nchunk)
        do j = 1, nchunk
          write(10,'(6F13.6)') r(1,j), r(2,j), r(3,j), v(1,j), v(2,j), v(3,j)
        enddo
        deallocate(r, v)
      endif
    enddo
 
    close(10)
    if (errorflag) then 
      exit
    endif

  enddo

end program prepare_lightcone

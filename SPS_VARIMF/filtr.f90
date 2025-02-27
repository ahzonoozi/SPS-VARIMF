SUBROUTINE Filtr()

  ! Routine to read and process all filters
  USE params; USE sub_func
  IMPLICIT NONE


  ! Constants and variables
  INTEGER :: stat=1,i,j,jj,i_min,i_max
  REAL(KIND(1.d0)), DIMENSION(50000) :: filtr_lm=0.,filtr_se=0.

  CHARACTER(3) :: filtr_type

  !----------------------------------------------------------------!

  ! Loop over all filters
  DO i=1,n_filtrs
      ! Generate the filter file name     
     WRITE(filtr_type,'(I3.3)') i

      ! Open the filter file
     OPEN(110,FILE='../FILTERS/filter_'//TRIM(filtr_type)//'.txt',&
        STATUS = 'OLD', IOSTAT = stat, ACTION = 'READ')
         
         
    ! Check if the file opened successfully
    IF (stat /= 0) THEN
      WRITE(*, *) 'ERROR: Unable to open filter file ', i
      STOP
    ENDIF

    ! Read the header line (assuming one header line exists)
    READ(110, *)

     
    ! Reset arrays for filter data
    jj=0
    filtr_lm = 0.0
    filtr_se = 0.0

    ! Read the filter data
    DO j=1,50000
        READ(110, *, IOSTAT = stat) filtr_lm(j), filtr_se(j)
        IF (stat /= 0) EXIT  ! Exit loop if end of file or error
      jj = j

      ! Ensure the filter sensitivity is non-negative
      IF (filtr_se(j) < 0.0) filtr_se(j) = 0.0
    ENDDO
    
    
    ! Validate the number of entries read
    IF (jj == 50000 .OR. jj == 1) THEN
      WRITE(*, *) 'ERROR: Problem reading filter ', i
      STOP
    ENDIF


    ! Interpolate the filter onto the wavelength array
    i_min = locate(lambda, filtr_lm(1))
    IF (i_min <= 1) i_min = 1
    i_max = locate(lambda, filtr_lm(jj))

    IF (i_min .NE. i_max) THEN
      Slm_filtr(i_min:i_max, i) = &
          linterparr(filtr_lm(1:jj), filtr_se(1:jj), lambda(i_min:i_max))
    ENDIF


    ! Ensure filter values are non-negative
    Slm_filtr(:, i) = MAX(Slm_filtr(:, i), 0.0)

     
    CLOSE(110)
  ENDDO


END SUBROUTINE Filtr



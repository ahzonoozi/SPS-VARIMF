FUNCTION Locate(values, target)

  ! Routine to locate the index of the interval in values that contains target

  USE params
  IMPLICIT NONE

  REAL(KIND(1.d0)), DIMENSION(:), INTENT(IN) :: values
  REAL(KIND(1.d0)), INTENT(IN) :: target
  INTEGER :: Locate
  INTEGER :: sizeValues, lower, middle, upper
  LOGICAL :: isAscending

  ! Get the size of the input array
  sizeValues = SIZE(values)
  isAscending = (values(sizeValues) >= values(1))

  ! Initialize indices for binary search
  lower = 0
  upper = sizeValues + 1

  ! Perform binary search
  DO
     IF (upper - lower <= 1) EXIT
     middle = (upper + lower) / 2
     IF (isAscending .EQV. (target >= values(middle))) THEN
        lower = middle
     ELSE
        upper = middle
     END IF
  END DO

  ! Determine the result index
  IF (target == values(1)) THEN
     Locate = 1
  ELSE IF (target == values(sizeValues)) THEN
     Locate = sizeValues - 1
  ELSE
     Locate = lower
  END IF

END FUNCTION Locate


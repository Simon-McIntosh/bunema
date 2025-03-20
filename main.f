      program main
      implicit none
      
      ! Declare variables
      real*8 :: psi(100), sia(100)
      integer :: nwb, nhb, nwnh, i

      ! Initialize variables
      nwb = 10
      nhb = 10
      nwnh = 100

      ! Fill psi and sia with random values to avoid zero division
      call random_seed()
      do i = 1, 100
        call random_number(psi(i))
        call random_number(sia(i))
      end do

      ! Call the subroutine
      call buneto(psi, nwb, nhb, sia, nwnh)

      ! End the program
      stop
      end program main

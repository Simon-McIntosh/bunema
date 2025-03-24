      program main
      implicit none
      
      ! Declare variables
      integer, parameter :: nwnh0 = 8385
      real*8 :: psi(nwnh0), sia(nwnh0)
      real*8 :: s,shift,dr,dz
      integer :: m,n
      common/bunemn/m,n,s,shift,dr,dz
      integer :: nwb, nhb, nwnh
      integer :: i, nbunema = 55

      ! Initialize variables
      open (unit=nbunema,file='bunema_in.dat',form='formatted')    
      read (nbunema,2022) m,n,nwb,nhb,nwnh
      write(*,*) 'nwnh =', nwnh
      read (nbunema,2020) s,shift,dr,dz
      read (nbunema,2020) (psi(i),i=1,nwnh)
      read (nbunema,2022) nwnh
      write(*,*) 'nwnh =', nwnh
      read (nbunema,2020) (sia(i),i=1,nwnh)
      read (nbunema,2022) nwnh
      write(*,*) 'nwnh =', nwnh
      close(nbunema)



      ! Fill psi and sia with random values to avoid zero division
      !call random_seed()
      !do i = 1, 100
      !  call random_number(psi(i))
      !  call random_number(sia(i))
      !end do

      ! Call the subroutine
      call buneto(psi, nwb, nhb, sia, nwnh)


2020  format (5e16.9)
2022  format (5i5)

      ! End the program
      stop
      end program main

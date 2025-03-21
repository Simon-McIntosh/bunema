      program main
      implicit none
      
      ! Declare variables
      integer, parameter :: nwnh = 8385
      real*8 :: psi(nwnh), sia(nwnh)
      real*8 :: s,shift,dr,dz
      integer :: m,n
      common/bunemn/m,n,s,shift,dr,dz
      integer :: nwb, nhb, i
      integer :: nbunema = 55

      ! Initialize variables
      open (unit=nbunema,file='bunema_in.dat',form='formatted')    
      read (nbunema,2022) m,n,nwb,nhb
      read (nbunema,2020) s,shift,dr,dz
      read (nbunema,2020) (psi(i),i=1,nwnh)
      close(nbunema)

      ! Call the subroutine
      call buneto(psi, nwb, nhb, sia, nwnh)

      ! Write the output to a file
      open (unit=nbunema,file='bunema_out.dat',form='formatted')
      write (nbunema,2022) m,n,nwb,nhb
      write (nbunema,2020) s,shift,dr,dz
      write (nbunema,2020) (psi(i),i=1,nwnh)
      write (nbunema,2020) (sia(i),i=1,nwnh)
      close(nbunema)


2020  format (4e16.9)
2022  format (5i5)

      ! End the program
      stop
      end program main

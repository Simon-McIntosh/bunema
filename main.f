      program bunema
      include 'double.inc'
c**********************************************************************
c**                                                                  **
c**     MAIN PROGRAM:  MHD FITTING CODE                              **
c**                                                                  **
c**     PROGRAM DESCRIPTION:                                         **
c**          This is the main driver for the MHD fitting code that   **
c**          reads input data, calls the buneto subroutine to        **
c**          process data with Buneman's solver, and outputs results **
c**                                                                  **
c**     REFERENCES:                                                  **
c**          (1)                                                     **
c**          (2)                                                     **
c**                                                                  **
c**     RECORD OF MODIFICATION:                                      **
c**          01 04/23..........first created                         **
c**                                                                  **
c**                                                                  **
c**********************************************************************
      parameter (NWB=64, NHB=64)
      parameter (NWNH=NWB*NHB)
      
      dimension psi(NWNH), sia(NWNH)
      dimension refdata(NWNH)
      
      common/bunemn/m,n,s,shift,dr,dz
c----------------------------------------------------------------------
c Initialize parameters                                              --
c----------------------------------------------------------------------
      m = NWB
      n = NHB
      s = 1.0d0
      shift = 0.5d0
      dr = 0.05d0
      dz = 0.05d0
      
c----------------------------------------------------------------------
c Read input data                                                    --
c----------------------------------------------------------------------
      open(10, file='bunema_in.dat', status='old')
      read(10,2020) (psi(i), i=1,NWNH)
2020  format(5e16.9)
      close(10)
      
c----------------------------------------------------------------------
c Process with Buneman solver                                        --
c----------------------------------------------------------------------
      call buneto(psi,NWB,NHB,sia,NWNH)

c----------------------------------------------------------------------
c Write output                                                       --
c----------------------------------------------------------------------
      open(11, file='bunema_out.dat', status='unknown')
      write(11,2020) (psi(i), i=1,NWNH)
      close(11)
      
c----------------------------------------------------------------------
c Compare with reference data (if available)                         --
c----------------------------------------------------------------------

      open(12, file='bunema_ref.dat', status='old', err=100)
      read(12,2020) (refdata(i), i=1,NWNH)
      close(12)
      
      maxdiff = 0.0d0
      do 30 i = 1,NWNH
         diff = abs(psi(i) - refdata(i))
         if (diff .gt. maxdiff) maxdiff = diff
30    continue
      
      write(*,*) 'Maximum difference from reference: ', maxdiff
      
100   continue
      stop
      end

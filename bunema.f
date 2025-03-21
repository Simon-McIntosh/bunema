      subroutine buneto(psi,nwb,nhb,sia,nwnh)
	  include 'double.inc'
c**********************************************************************
c**                                                                  **
c**     MAIN PROGRAM:  MHD FITTING CODE                              **
c**                                                                  **
c**                                                                  **
c**     SUBPROGRAM DESCRIPTION:                                      **
c**          buneto sets up the appropriate arrays for the           **
c**          Buneman's solver.                                       **
c**                                                                  **
c**     CALLING ARGUMENTS:                                           **
c**                                                                  **
c**     REFERENCES:                                                  **
c**          (1)                                                     **
c**          (2)                                                     **
c**                                                                  **
c**     RECORD OF MODIFICATION:                                      **
c**          05 03/85..........first created                         **
c**                                                                  **
c**                                                                  **
c**********************************************************************
      dimension   psi(nwnh), sia(nwnh)
      common/bunemn/m,n,s,shift,dr,dz
c----------------------------------------------------------------------
c copy psi into sia rowwise                                          --
c----------------------------------------------------------------------
      do 2 i = 1,nwb
      ii = (i-1)*nhb + 1
      do 2 j = 1, nhb
      sia(i+j*nwb-nwb) = psi(ii-1+j)
    2 continue
      ia = nwb+nwb
      ju = n*nwb
c-----------------------------
c set up for rzpois
c-----------------------------
      do 3 i = ia,ju,nwb
      sia(i-m+1) = sia(i-m+1)+(.5+.25/(1.+shift/dr))*sia(i-m)/s
      sia(i-1) = sia(i-1)+(.5-.25/(m-1+shift/dr))*sia(i)/s
    3 continue
      call rzpois(sia,nwnh)
      nwhbb = nwb*nhb
c------------------------------------------------------------------------------
c copy sia back into psi columnwise                                          --
c------------------------------------------------------------------------------
      do 6 i = 2,n
      ii = (i-1)*nwb + 1
      do 6 j = 2,m
      psi(i+j*nhb-nhb) = sia(ii-1+j)
    6 continue
      return
      end
      subroutine rzpois(q,nwnh)
c**********************************************************************
c**                                                                  **
c**     MAIN PROGRAM:  MHD FITTING CODE                              **
c**                                                                  **
c**                                                                  **
c**     SUBPROGRAM DESCRIPTION:                                      **
c**          rzpois solves for the poloidal flux using the           **
c**          Buneman's method.                                       **
c**                                                                  **
c**     CALLING ARGUMENTS:                                           **
c**                                                                  **
c**     REFERENCES:                                                  **
c**          (1)                                                     **
c**          (2)                                                     **
c**                                                                  **
c**     RECORD OF MODIFICATION:                                      **
c**          05 03/85..........first created                         **
c**                                                                  **
c**                                                                  **
c**********************************************************************
	  implicit real *8 (a-h,o-z)

      common/bunemn/m,n,s,shift,dr,dz
      dimension   g(300),p(300),c(300),d(300),temp(300)

      dimension  q(nwnh)
      include 'double_bunema.inc'

      do 50 i = 1,300
      g(i) = 0.
      p(i) = 0.
      d(i) = 0.
   50 continue
c      if (flag .ne. 0.) go to 200
      do 60 i = 1,300
      temp(i) = 0.
      c(i) = 0.
   60 continue
      shftdr = shift/dr
      do 40 i = 2,m
      temp(i) = 1. - .5/(i+shftdr-1.)
   40 continue
      ju = (n-1)*(m+1)
      n222 = n/2
      c(n222) = 0.
      lo = n/2
   1  l = lo/2
      c(l) = sqrt(2.+c(lo))
      lo = l
   2  c(n-l) = -c(l)
      l = l+2*lo
      if ((2*l/n)*(2*lo-3))4,3,1
   3  c(l) = (c(l+lo)+c(l-lo))/c(lo)
      go to 2
   4  do 5 l = 2,n
   5  c(l-1) = 1./(2.+s*(2.-c(l-1)))
      flag = 1.
  200 lo = n/2
      ko = 2
      id = 1
  15  li = 2*lo
      k4 = 2*ko-li/n
      jd = (m+1)*n/li
      jh = (m+1)*(n/(2*li))
      jt = jd+jh
      ji = 2*jd
      jo = jd*ko
      do 11 j = jo,ju,ji
      j2 = j+2
      iu = j+m
      go to (20,24,26,28),k4
  28  do 29 i = j2,iu
      pi = q(i)-q(i+jt)-q(i-jt)
      q(i) = q(i)-q(i+jh)-q(i-jh)+q(i+jd)+q(i-jd)
  29  p(i-j) = pi+q(i)
      go to 10
  26  do 27 i = j2,iu
      p(i-j) = 2.*q(i)
  27  q(i) = q(i+jd)+q(i-jd)
      go to 10
  24  do 25 i = j2,iu
      p(i-j) = 2.*q(i)+q(i+jd)+q(i-jd)
  25  q(i) = q(i)-q(i+jh)-q(i-jh)
      go to 10
  20  do 23 i = j2,iu
      p(i-j) = 2.*q(i)+q(i+jd)+q(i-jd)
  23  q(i) = 0.
  10  do 22 l = lo,n,li
      a = c(l)
      as = a*s
      do 18 i = 2,m
      p(i) = as*p(i)
      d(i) = a*temp(i)
  18  g(i) = 2*a-d(i)
      g(2) = 0.
      d(m) = 0.
  19  ii = 2*id
      io = ii+1
      do 21 i = io,m,ii
      a = 1./(1.-d(i)*g(i+id)-g(i)*d(i-id))
      p(i) = a*(p(i)+d(i)*p(i+id)+g(i)*p(i-id))
      d(i) = d(i)*d(i+id)*a
  21  g(i) = g(i)*g(i-id)*a
      id = ii
      if (id-m/2.lt.0)go to 19
  16  id = ii/2
      io = id+1
      do 17 i = io,m,ii
  17  p(i) = p(i)+d(i)*p(i+id)+g(i)*p(i-id)
      ii = id
      if (id.gt.1)go to 16
  22  continue
      do 11 i = j2,iu
  11  q(i) = q(i)+p(i-j)
      go to (13,12),ko
  12  lo = lo/2
      if (lo.eq.1)ko = 1
      go to 15
  13  lo = 2*lo
      if (lo.lt.n)go to 15
      return
      end


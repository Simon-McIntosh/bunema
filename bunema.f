      subroutine buneto(psi, nwb, nhb, sia, nwnh)
            implicit none
            double precision, dimension(:), intent(in) :: psi
            double precision, dimension(:), intent(out) :: sia
            integer, intent(in) :: nwb, nhb, nwnh
            include 'double.inc'
            integer :: ia, ju
            interface
                  subroutine rzpois(q, nwnh)
                        implicit none
                        double precision, dimension(:), intent(inout) :: q
                        integer, intent(in) :: nwnh
                  end subroutine rzpois
            end interface
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
            double precision, dimension(nwnh) :: psi, sia
            common /bunemn/ m, n, s, shift, dr, dz
            integer :: m, n
            double precision :: s, shift, dr, dz

c Copy psi into sia row-wise
            call copy_psi_to_sia(psi, sia, nwb, nhb)

            ia = nwb + nwb
            ju = n * nwb

c Set up for rzpois
            call setup_rzpois(sia, ia, ju, nwb, m, shift, dr, s)

            call rzpois(sia, nwnh)

c Copy sia back into psi column-wise
            call copy_sia_to_psi(sia, psi, n, nwb, nhb, m)

            return
            end subroutine buneto

            subroutine rzpois(q, nwnh)
            implicit none
            double precision, dimension(:), intent(inout) :: q
            integer, intent(in) :: nwnh
            integer :: ju, n222, lo, ko, id, li, jd, jh, jt, ji, jo, j2, iu, i, j, k4
            double precision :: shftdr
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
            double precision, dimension(300) :: g, p, c, d, temp
            double precision, dimension(nwnh) :: q
            include 'double_bunema.inc'

c Initialize arrays
            call initialize_arrays(g, p, 300)
            call initialize_arrays(temp, c, 300)

            shftdr = shift / dr

c Compute temp array values
            call compute_temp(temp, m, shftdr)

            ju = (n - 1) * (m + 1)
            n222 = n / 2
            c(n222) = 0.d0
            lo = n / 2

c Compute c array values
            call compute_c(c, lo, n, s)

            lo = n / 2
            ko = 2
            id = 1

c Main computation loop
            do while (lo >= 1)
                  li = 2 * lo
                  k4 = 2 * ko - li / n
                  jd = (m + 1) * n / li
                  jh = (m + 1) * (n / (2 * li))
                  jt = jd + jh
                  ji = 2 * jd
                  jo = jd * ko

                  do j = jo, ju, ji
                        j2 = j + 2
                        iu = j + m

                        select case (k4)
                              case (28)
                                    call case_28(q, p, j2, iu, jt, jh, jd)
                              case (26)
                                    call case_26(q, p, j2, iu, jd)
                              case (24)
                                    call case_24(q, p, j2, iu, jd, jh)
                              case (20)
                                    call case_20(q, p, j2, iu, jd)
                        end select

c Update arrays
                        call update_arrays(g, p, d, temp, c, lo, n, li, m, s, id)

                        do i = j2, iu
                              q(i) = q(i) + p(i - j)
                        end do

                        if (ko == 1) then
                              exit
                        end if
                  end do

                  lo = lo / 2
                  if (lo == 1) then
                        ko = 1
                  end if
            end do

            return
            end subroutine rzpois

            subroutine copy_psi_to_sia(psi, sia, nwb, nhb)
            implicit none
            integer, intent(in) :: nwb, nhb
            double precision, dimension(:), intent(in) :: psi
            double precision, dimension(:), intent(out) :: sia
            integer :: i, j, ii

            do i = 1, nwb
                  ii = (i - 1) * nhb + 1
                  do j = 1, nhb
                        sia(i + j * nwb - nwb) = psi(ii - 1 + j)
                  end do
            end do
            end subroutine copy_psi_to_sia

            subroutine setup_rzpois(sia, ia, ju, nwb, m, shift, dr, s)
            implicit none
            integer, intent(in) :: ia, ju, nwb, m
            double precision, intent(in) :: shift, dr, s
            double precision, dimension(:), intent(inout) :: sia
            integer :: i

            do i = ia, ju, nwb
                  sia(i - m + 1) = sia(i - m + 1) + (.5d0 + .25d0 / (1.d0 + shift / dr)) * sia(i - m) / s
                  sia(i - 1) = sia(i - 1) + (.5d0 - .25d0 / (m - 1 + shift / dr)) * sia(i) / s
            end do
            end subroutine setup_rzpois

            subroutine copy_sia_to_psi(sia, psi, n, nwb, nhb, m)
            implicit none
            integer, intent(in) :: n, nwb, nhb, m
            double precision, dimension(:), intent(in) :: sia
            double precision, dimension(:), intent(out) :: psi
            integer :: i, j, ii

            do i = 2, n
                  ii = (i - 1) * nwb + 1
                  do j = 2, m
                        psi(i + j * nhb - nhb) = sia(ii - 1 + j)
                  end do
            end do
            end subroutine copy_sia_to_psi

            subroutine initialize_arrays(arr1, arr2, size)
            implicit none
            integer, intent(in) :: size
            double precision, dimension(size), intent(out) :: arr1, arr2
            integer :: i

            do i = 1, size
                  arr1(i) = 0.d0
                  arr2(i) = 0.d0
            end do
            end subroutine initialize_arrays

            subroutine compute_temp(temp, m, shftdr)
            implicit none
            integer, intent(in) :: m
            double precision, intent(in) :: shftdr
            double precision, dimension(m), intent(out) :: temp
            integer :: i

            do i = 2, m
                  temp(i) = 1.d0 - .5d0 / (i + shftdr - 1.d0)
            end do
            end subroutine compute_temp

            subroutine compute_c(c, lo, n, s)
            implicit none
            integer, intent(in) :: lo, n
            double precision, intent(in) :: s
            double precision, dimension(n), intent(out) :: c
            integer :: l

            l = lo / 2
            c(l) = sqrt(2.d0 + c(lo))
            lo = l
            c(n - l) = -c(l)
            l = l + 2 * lo
            if ((2 * l / n) * (2 * lo - 3)) then
                  c(l) = (c(l + lo) + c(l - lo)) / c(lo)
            end if
            c(l - 1) = 1.d0 / (2.d0 + s * (2.d0 - c(l - 1)))
            end subroutine compute_c

            subroutine case_28(q, p, j2, iu, jt, jh, jd)
            implicit none
            integer, intent(in) :: j2, iu, jt, jh, jd
            double precision, dimension(:), intent(inout) :: q, p
            double precision :: pi
            integer :: i

            do i = j2, iu
                  pi = q(i) - q(i + jt) - q(i - jt)
                  q(i) = q(i) - q(i + jh) - q(i - jh) + q(i + jd) + q(i - jd)
                  p(i - j2) = pi + q(i)
            end do
            end subroutine case_28

            subroutine case_26(q, p, j2, iu, jd)
            implicit none
            integer, intent(in) :: j2, iu, jd
            double precision, dimension(:), intent(inout) :: q, p
            integer :: i

            do i = j2, iu
                  p(i - j2) = 2.d0 * q(i)
                  q(i) = q(i + jd) + q(i - jd)
            end do
            end subroutine case_26

            subroutine case_24(q, p, j2, iu, jd, jh)
            implicit none
            integer, intent(in) :: j2, iu, jd, jh
            double precision, dimension(:), intent(inout) :: q, p
            integer :: i

            do i = j2, iu
                  p(i - j2) = 2.d0 * q(i) + q(i + jd) + q(i - jd)
                  q(i) = q(i) - q(i + jh) - q(i - jh)
            end do
            end subroutine case_24

            subroutine case_20(q, p, j2, iu, jd)
            implicit none
            integer, intent(in) :: j2, iu, jd
            double precision, dimension(:), intent(inout) :: q, p
            integer :: i

            do i = j2, iu
                  p(i - j2) = 2.d0 * q(i) + q(i + jd) + q(i - jd)
                  q(i) = 0.d0
            end do
            end subroutine case_20

            subroutine update_arrays(g, p, d, temp, c, lo, n, li, m, s, id)
            implicit none
            integer, intent(in) :: lo, n, li, m
            double precision, intent(in) :: s
            double precision, dimension(:), intent(inout) :: g, p, d, temp, c
            integer, intent(inout) :: id
            integer :: l, i, ii, io
            double precision :: a, as

            do l = lo, n, li
                  a = c(l)
                  as = a * s
                  do i = 2, m
                        p(i) = as * p(i)
                        d(i) = a * temp(i)
                        g(i) = 2.d0 * a - d(i)
                  end do
                  g(2) = 0.d0
                  d(m) = 0.d0

                  ii = 2 * id
                  io = ii + 1
                  do i = io, m, ii
                        a = 1.d0 / (1.d0 - d(i) * g(i + id) - g(i) * d(i - id))
                        p(i) = a * (p(i) + d(i) * p(i + id) + g(i) * p(i - id))
                        d(i) = d(i) * d(i + id) * a
                        g(i) = g(i) * g(i - id) * a
                  end do
                  id = ii
                  if (id - m / 2 < 0) then
                        exit
                  end if

                  id = ii / 2
                  io = id + 1
                  do i = io, m, ii
                        p(i) = p(i) + d(i) * p(i + id) + g(i) * p(i - id)
                  end do
                  ii = id
                  if (id > 1) then
                        exit
                  end if
            end do
            end subroutine update_arrays


! PUBFFT is used if nothing else is defined.

!*******************************************************************************
!
!     All of the ewald code and supporting routines were written and 
!     contributed by Tom Darden from the National Institute of 
!     Environmental Health Sciences division of the NIH.  
!     Originally written with a modified version of AMBER 3A, the code 
!     was updated during the summer of 1994 to be compatible with 
!     AMBER 4.1.
!
!     This section contains code necessary to perform 3D FFTs where
!     libraries are not available.  It is based on piecing together a
!     series of 1D FFTs and is probably not super efficient.  The 1D
!     FFT code is a double precision version of fftpack from netlib,
!     written by Paul N. Swartztrauber at NCAR Boulder Colorado.
!
!     The following routines are defined:
!
!     cffti     initialize cfftf and cfftb
!     cffti1
!
!     cfftb     unnormalized inverse of cfftf
!     cfftb1
!     passb4,3,2,5
!     passb
!
!     cfftf     forward transform of a complex periodic sequence
!     cfftf1
!     passf4,3,2,5
!     passf
!
!*******************************************************************************

!*******************************************************************************
! subroutine cffti(n,wsave,ifac)
!
! subroutine cffti initializes the array wsave which is used in
! both cfftf and cfftb. the prime factorization of n together with
! a tabulation of the trigonometric functions are computed and
! stored in wsave.
!
! input parameter
!
! n       the length of the sequence to be transformed
!
! output parameter
!
! wsave   A work array which must be dimensioned at least 4*n.
!         The same work array can be used for both cfftf and cfftb
!         as long as n remains unchanged. different wsave arrays
!         are required for different values of n. the contents of
!         wsave must not be changed between calls of cfftf or cfftb.
! ifac    The saved factors array, dimensioned 30.  This was originally
!         part of wsave, but has been separated to make it possible to
!         put the code into a module.
! 
!*******************************************************************************

subroutine cffti(n, wsave, ifac)
      
  implicit none

! Formal arguments:

  integer               :: n
  double precision      :: wsave(*)
  integer               :: ifac(*)

  if (n .ne. 1) call cffti1(n, wsave(2 * n + 1), ifac)

  return

end subroutine cffti

!*******************************************************************************
!
! Internal Subroutine:  cffti1
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine cffti1(n, wa, ifac)
      
  implicit none

! Formal arguments:

  integer               :: n
  double precision      :: wa(*)
  integer               :: ifac(*)

! Local variables:

  integer               :: ntryh(4), nl, nf, j, ntry, nq, nr, i, ib, l1, &
                           k1, ip, ld, l2, ido, idot, ipm, i1, ii

  double precision      :: tpi, argh, fi, argld, arg

  data ntryh(1), ntryh(2), ntryh(3), ntryh(4) /3, 4, 2, 5/
  data i, j /0, 0/  ! Related to bug #248 in pmemd, fixed by JMS 

  nl = n
  nf = 0
  j = 1
  ntry = ntryh(1)

  do while (nl .ne. 1)

    nq = nl/ntry
    nr = nl - ntry*nq

    if (nr .ne. 0) then

      j = j + 1
      if (j .le. 4) then
        ntry = ntryh(j)
      else
        ntry = ntry + 2
      end if

      cycle

    end if

    nf = nf + 1
    ifac(nf + 2) = ntry
    nl = nq

    if (ntry .eq. 2 .and. nf .ne. 1) then
      do i = 2, nf
        ib = nf - i + 2
        ifac(ib + 2) = ifac(ib + 1)
      end do
      ifac(3) = 2
    end if

  end do

  ifac(1) = n
  ifac(2) = nf
  tpi = 6.28318530717959d0
  argh = tpi/dble(n)
  i = 2
  l1 = 1
  do k1 = 1, nf
    ip = ifac(k1 + 2)
    ld = 0
    l2 = l1*ip
    ido = n/l2
    idot = ido + ido + 2
    ipm = ip - 1
    do j = 1, ipm
      i1 = i
      wa(i - 1) = 1.d0
      wa(i) = 0.d0
      ld = ld+l1
      fi = 0.d0
      argld = dble(ld)*argh
      do ii = 4, idot, 2
        i = i + 2
        fi = fi + 1.d0
        arg = fi*argld
        wa(i - 1) = cos(arg)
        wa(i) = sin(arg)
      end do
      if (ip .le. 5) cycle
      wa(i1 - 1) = wa(i - 1)
      wa(i1) = wa(i)
    end do
    l1 = l2
  end do

  return

end subroutine cffti1

!*******************************************************************************
! subroutine cfftb(n,c,wsave,ifac)
!
! subroutine cfftb computes the backward complex discrete fourier
! transform (the fourier synthesis). equivalently , cfftb computes
! a complex periodic sequence from its fourier coefficients.
! the transform is defined below at output parameter c.
!
! a call of cfftf followed by a call of cfftb will multiply the
! sequence by n.
!
! the array wsave which is used by subroutine cfftb must be
! initialized by calling subroutine cffti(n,wsave,ifac).
!
! input parameters
!
! n      the length of the complex sequence c. the method is
!        more efficient when n is the product of small primes.
!
! c      a complex array of length n which contains the sequence
!
! wsave   a real work array which must be dimensioned at least 4n
!         in the program that calls cfftb. the wsave array must be
!         initialized by calling subroutine cffti(n,wsave) and a
!         different wsave array must be used for each different
!         value of n. this initialization does not have to be
!         repeated so long as n remains unchanged thus subsequent
!         transforms can be obtained faster than the first.
!         the same wsave array can be used by cfftf and cfftb.
!
! ifac    The saved factors array, dimensioned 30.  This was originally
!         part of wsave, but has been separated to make it possible to
!         put the code into a module.  It was created by cffti().
! 
! output parameters
!
! c      for j=1,...,n
!
!            c(j)=the sum from k=1,...,n of
!
!                  c(k)*exp(i*(j-1)*(k-1)*2*pi/n)
!
!                        where i=sqrt(-1)
!
! wsave   contains initialization calculations which must not be
!         destroyed between calls of subroutine cfftf or cfftb
!
! ifac    contains initialization calculations which must not be
!         destroyed between calls of subroutine cfftf or cfftb
!
!*******************************************************************************

subroutine cfftb(n, c, wsave, ifac)
      
  implicit none

! Formal arguments:

  integer               :: n
  double precision      :: c(*)
  double precision      :: wsave(*)
  integer               :: ifac(*)

  if (n .ne. 1) call cfftb1(n, c, wsave, wsave(2 * n + 1), ifac)

  return

end subroutine cfftb

!*******************************************************************************
!
! Subroutine:  cfftb1
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine cfftb1(n, c, ch, wa, ifac)
  
  implicit none

! Formal arguments:

  integer               :: n
  double precision      :: c(*)
  double precision      :: ch(*)
  double precision      :: wa(*)
  integer               :: ifac(*)

! Local variables:

  integer               :: nf, na, l1, iw, k1, ip, l2, ido, idot, idl1, &
                           ix2, ix3, ix4, nac, n2

  nf = ifac(2)
  na = 0
  l1 = 1
  iw = 1

  do k1 = 1, nf

    ip = ifac(k1 + 2)
    l2 = ip*l1
    ido = n/l2
    idot = ido + ido
    idl1 = idot*l1

    if (ip .eq. 4) then

      ix2 = iw + idot
      ix3 = ix2 + idot

      if (na .ne. 0) then
        call passb4(idot, l1, ch, c, wa(iw), wa(ix2), wa(ix3))
      else
        call passb4(idot, l1, c, ch, wa(iw), wa(ix2), wa(ix3))
      end if

      na = 1 - na

    else if (ip .eq. 2) then

      if (na .ne. 0) then
        call passb2(idot, l1, ch, c, wa(iw))
      else
        call passb2(idot, l1, c, ch, wa(iw))
      end if

      na = 1 - na

    else if (ip .eq. 3) then

      ix2 = iw + idot

      if (na .ne. 0) then
        call passb3(idot, l1, ch, c, wa(iw), wa(ix2))
      else
        call passb3(idot, l1, c, ch, wa(iw), wa(ix2))
      end if

      na = 1 - na

    else if (ip .eq. 5) then

      ix2 = iw + idot
      ix3 = ix2 + idot
      ix4 = ix3 + idot

      if (na .ne. 0) then
        call passb5(idot, l1, ch, c, wa(iw), wa(ix2), wa(ix3), wa(ix4))
      else
        call passb5(idot, l1, c, ch, wa(iw), wa(ix2), wa(ix3), wa(ix4))
      end if

      na = 1 - na

    else

      if (na .ne. 0) then
        call passb(nac, idot, ip, l1, idl1, ch, ch, ch, c, c, wa(iw))
      else
        call passb(nac, idot, ip, l1, idl1, c, c, c, ch, ch, wa(iw))
      end if

      if (nac .ne. 0) na = 1 - na

    end if
     
    l1 = l2
    iw = iw + (ip - 1)*idot

  end do

  if (na .eq. 0) return

  n2 = n + n

  c(1:n2) = ch(1:n2)

  return

contains

!*******************************************************************************
!
! Internal Subroutine:  passb4
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine passb4(ido, l1, cc, ch, wa1, wa2, wa3)
      
  implicit none

! Formal arguments:

  integer               :: ido
  integer               :: l1
  double precision      :: cc(ido, 4, l1)           
  double precision      :: ch(ido, l1, 4)           
  double precision      :: wa1(*)
  double precision      :: wa2(*)
  double precision      :: wa3(*)

! Local variables:

  integer               :: k, i
  double precision      :: ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4, &
                           ci2, ci3, ci4, cr2, cr3, cr4

  if (ido .eq. 2) then

    do k = 1, l1
      ti1 = cc(2, 1, k) - cc(2, 3, k)
      ti2 = cc(2, 1, k) + cc(2, 3, k)
      tr4 = cc(2, 4, k) - cc(2, 2, k)
      ti3 = cc(2, 2, k) + cc(2, 4, k)
      tr1 = cc(1, 1, k) - cc(1, 3, k)
      tr2 = cc(1, 1, k) + cc(1, 3, k)
      ti4 = cc(1, 2, k) - cc(1, 4, k)
      tr3 = cc(1, 2, k) + cc(1, 4, k)
      ch(1, k, 1) = tr2 + tr3
      ch(1, k, 3) = tr2 - tr3
      ch(2, k, 1) = ti2 + ti3
      ch(2, k, 3) = ti2 - ti3
      ch(1, k, 2) = tr1 + tr4
      ch(1, k, 4) = tr1 - tr4
      ch(2, k, 2) = ti1 + ti4
      ch(2, k, 4) = ti1 - ti4
    end do

  else
  
    do k = 1, l1
      do i = 2, ido, 2
        ti1 = cc(i, 1, k) - cc(i, 3, k)
        ti2 = cc(i, 1, k) + cc(i, 3, k)
        ti3 = cc(i, 2, k) + cc(i, 4, k)
        tr4 = cc(i, 4, k) - cc(i, 2, k)
        tr1 = cc(i - 1, 1, k) - cc(i - 1, 3, k)
        tr2 = cc(i - 1, 1, k) + cc(i - 1, 3, k)
        ti4 = cc(i - 1, 2, k) - cc(i - 1, 4, k)
        tr3 = cc(i - 1, 2, k) + cc(i - 1, 4, k)
        ch(i - 1, k, 1) = tr2 + tr3
        cr3 = tr2 - tr3
        ch(i, k, 1) = ti2 + ti3
        ci3 = ti2 - ti3
        cr2 = tr1 + tr4
        cr4 = tr1 - tr4
        ci2 = ti1 + ti4
        ci4 = ti1 - ti4
        ch(i - 1, k, 2) = wa1(i - 1)*cr2 - wa1(i)*ci2
        ch(i, k, 2) = wa1(i - 1)*ci2 + wa1(i)*cr2
        ch(i - 1, k, 3) = wa2(i - 1)*cr3 - wa2(i)*ci3
        ch(i, k, 3) = wa2(i - 1)*ci3 + wa2(i)*cr3
        ch(i - 1, k, 4) = wa3(i - 1)*cr4 - wa3(i)*ci4
        ch(i, k, 4) = wa3(i - 1)*ci4 + wa3(i)*cr4
      end do
    end do

  end if

  return

end subroutine passb4

!*******************************************************************************
!
! Internal Subroutine:  passb3
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine passb3(ido, l1, cc, ch, wa1, wa2)
      
  implicit none

! Formal arguments:

  integer               :: ido
  integer               :: l1
  double precision      :: cc(ido, 3, l1)
  double precision      :: ch(ido, l1, 3)
  double precision      :: wa1(*)
  double precision      :: wa2(*)

! Local variables:

  integer               :: k, i

  double precision      :: taur, taui, tr2, cr2, ti2, ci2, cr3, ci3, &
                           dr2, dr3, di2, di3

  data taur, taui / - .5d0, .866025403784439d0/

  if (ido .eq. 2) then

    do k = 1, l1
      tr2 = cc(1, 2, k) + cc(1, 3, k)
      cr2 = cc(1, 1, k) + taur*tr2
      ch(1, k, 1) = cc(1, 1, k) + tr2
      ti2 = cc(2, 2, k) + cc(2, 3, k)
      ci2 = cc(2, 1, k) + taur*ti2
      ch(2, k, 1) = cc(2, 1, k) + ti2
      cr3 = taui*(cc(1, 2, k) - cc(1, 3, k))
      ci3 = taui*(cc(2, 2, k) - cc(2, 3, k))
      ch(1, k, 2) = cr2 - ci3
      ch(1, k, 3) = cr2 + ci3
      ch(2, k, 2) = ci2 + cr3
      ch(2, k, 3) = ci2 - cr3
    end do

  else

    do k = 1, l1
      do i = 2, ido, 2
        tr2 = cc(i - 1, 2, k) + cc(i - 1, 3, k)
        cr2 = cc(i - 1, 1, k) + taur*tr2
        ch(i - 1, k, 1) = cc(i - 1, 1, k) + tr2
        ti2 = cc(i, 2, k) + cc(i, 3, k)
        ci2 = cc(i, 1, k) + taur*ti2
        ch(i, k, 1) = cc(i, 1, k) + ti2
        cr3 = taui*(cc(i - 1, 2, k) - cc(i - 1, 3, k))
        ci3 = taui*(cc(i, 2, k) - cc(i, 3, k))
        dr2 = cr2 - ci3
        dr3 = cr2 + ci3
        di2 = ci2 + cr3
        di3 = ci2 - cr3
        ch(i, k, 2) = wa1(i - 1)*di2 + wa1(i)*dr2
        ch(i - 1, k, 2) = wa1(i - 1)*dr2 - wa1(i)*di2
        ch(i, k, 3) = wa2(i - 1)*di3 + wa2(i)*dr3
        ch(i - 1, k, 3) = wa2(i - 1)*dr3 - wa2(i)*di3
      end do
    end do

  end if

  return

end subroutine passb3

!*******************************************************************************
!
! Internal Subroutine:  passb2
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine passb2(ido, l1, cc, ch, wa1)
      
  implicit none

! Formal arguments:

  integer               :: ido
  integer               :: l1
  double precision      :: cc(ido, 2, l1)
  double precision      :: ch(ido, l1, 2)
  double precision      :: wa1(*)

! Local variables:

  integer               :: k, i
  double precision      :: tr2, ti2

  if (ido .le. 2) then

    do k = 1, l1
      ch(1, k, 1) = cc(1, 1, k) + cc(1, 2, k)
      ch(1, k, 2) = cc(1, 1, k) - cc(1, 2, k)
      ch(2, k, 1) = cc(2, 1, k) + cc(2, 2, k)
      ch(2, k, 2) = cc(2, 1, k) - cc(2, 2, k)
    end do

  else

    do k = 1, l1
      do i = 2, ido, 2
        ch(i - 1, k, 1) = cc(i - 1, 1, k) + cc(i - 1, 2, k)
        tr2 = cc(i - 1, 1, k) - cc(i - 1, 2, k)
        ch(i, k, 1) = cc(i, 1, k) + cc(i, 2, k)
        ti2 = cc(i, 1, k) - cc(i, 2, k)
        ch(i, k, 2) = wa1(i - 1)*ti2 + wa1(i)*tr2
        ch(i - 1, k, 2) = wa1(i - 1)*tr2 - wa1(i)*ti2
      end do
    end do

  end if

  return

end subroutine passb2

!*******************************************************************************
!
! Internal Subroutine:  passb5
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine passb5(ido, l1, cc, ch, wa1, wa2, wa3, wa4)
      
  implicit none

! Formal arguments:

  integer               :: ido
  integer               :: l1
  double precision      :: cc(ido, 5, l1)           
  double precision      :: ch(ido, l1, 5)           
  double precision      :: wa1(*)
  double precision      :: wa2(*)
  double precision      :: wa3(*)
  double precision      :: wa4(*)

! Local variables:

  integer               :: i, k

  double precision      :: ti2, ti3, ti4, ti5, ti11, ti12
  double precision      :: tr2, tr3, tr4, tr5, tr11, tr12
  double precision      :: ci2, ci3, ci4, ci5
  double precision      :: cr2, cr3, cr4, cr5
  double precision      :: di2, di3, di4, di5
  double precision      :: dr2, dr3, dr4, dr5

  data tr11, ti11, tr12, ti12 /.309016994374947d0, &
                               .951056516295154d0, &
                             - .809016994374947d0, &
                               .587785252292473d0/

  if (ido .eq. 2) then

    do k = 1, l1
      ti5 = cc(2, 2, k) - cc(2, 5, k)
      ti2 = cc(2, 2, k) + cc(2, 5, k)
      ti4 = cc(2, 3, k) - cc(2, 4, k)
      ti3 = cc(2, 3, k) + cc(2, 4, k)
      tr5 = cc(1, 2, k) - cc(1, 5, k)
      tr2 = cc(1, 2, k) + cc(1, 5, k)
      tr4 = cc(1, 3, k) - cc(1, 4, k)
      tr3 = cc(1, 3, k) + cc(1, 4, k)
      ch(1, k, 1) = cc(1, 1, k) + tr2 + tr3
      ch(2, k, 1) = cc(2, 1, k) + ti2 + ti3
      cr2 = cc(1, 1, k) + tr11*tr2 + tr12*tr3
      ci2 = cc(2, 1, k) + tr11*ti2 + tr12*ti3
      cr3 = cc(1, 1, k) + tr12*tr2 + tr11*tr3
      ci3 = cc(2, 1, k) + tr12*ti2 + tr11*ti3
      cr5 = ti11*tr5 + ti12*tr4
      ci5 = ti11*ti5 + ti12*ti4
      cr4 = ti12*tr5 - ti11*tr4
      ci4 = ti12*ti5 - ti11*ti4
      ch(1, k, 2) = cr2 - ci5
      ch(1, k, 5) = cr2 + ci5
      ch(2, k, 2) = ci2 + cr5
      ch(2, k, 3) = ci3 + cr4
      ch(1, k, 3) = cr3 - ci4
      ch(1, k, 4) = cr3 + ci4
      ch(2, k, 4) = ci3 - cr4
      ch(2, k, 5) = ci2 - cr5
    end do

  else

    do k = 1, l1
      do i = 2, ido, 2
        ti5 = cc(i, 2, k) - cc(i, 5, k)
        ti2 = cc(i, 2, k) + cc(i, 5, k)
        ti4 = cc(i, 3, k) - cc(i, 4, k)
        ti3 = cc(i, 3, k) + cc(i, 4, k)
        tr5 = cc(i - 1, 2, k) - cc(i - 1, 5, k)
        tr2 = cc(i - 1, 2, k) + cc(i - 1, 5, k)
        tr4 = cc(i - 1, 3, k) - cc(i - 1, 4, k)
        tr3 = cc(i - 1, 3, k) + cc(i - 1, 4, k)
        ch(i - 1, k, 1) = cc(i - 1, 1, k) + tr2 + tr3
        ch(i, k, 1) = cc(i, 1, k) + ti2 + ti3
        cr2 = cc(i - 1, 1, k) + tr11*tr2 + tr12*tr3
        ci2 = cc(i, 1, k) + tr11*ti2 + tr12*ti3
        cr3 = cc(i - 1, 1, k) + tr12*tr2 + tr11*tr3
        ci3 = cc(i, 1, k) + tr12*ti2 + tr11*ti3
        cr5 = ti11*tr5 + ti12*tr4
        ci5 = ti11*ti5 + ti12*ti4
        cr4 = ti12*tr5 - ti11*tr4
        ci4 = ti12*ti5 - ti11*ti4
        dr3 = cr3 - ci4
        dr4 = cr3 + ci4
        di3 = ci3 + cr4
        di4 = ci3 - cr4
        dr5 = cr2 + ci5
        dr2 = cr2 - ci5
        di5 = ci2 - cr5
        di2 = ci2 + cr5
        ch(i - 1, k, 2) = wa1(i - 1)*dr2 - wa1(i)*di2
        ch(i, k, 2) = wa1(i - 1)*di2 + wa1(i)*dr2
        ch(i - 1, k, 3) = wa2(i - 1)*dr3 - wa2(i)*di3
        ch(i, k, 3) = wa2(i - 1)*di3 + wa2(i)*dr3
        ch(i - 1, k, 4) = wa3(i - 1)*dr4 - wa3(i)*di4
        ch(i, k, 4) = wa3(i - 1)*di4 + wa3(i)*dr4
        ch(i - 1, k, 5) = wa4(i - 1)*dr5 - wa4(i)*di5
        ch(i, k, 5) = wa4(i - 1)*di5 + wa4(i)*dr5
      end do
    end do

  end if

  return

end subroutine passb5

!*******************************************************************************
!
! Internal Subroutine:  passb
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine passb(nac, ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa)
      
  implicit none

! Formal arguments:

  integer               :: nac             
  integer               :: ido
  integer               :: ip
  integer               :: l1
  integer               :: idl1
  double precision      :: cc(ido, ip, l1)
  double precision      :: c1(ido, l1, ip)
  double precision      :: c2(idl1, ip)
  double precision      :: ch(ido, l1, ip)
  double precision      :: ch2(idl1, ip)
  double precision      :: wa(*)

! Local variables:

  integer               :: idot, nt, ipp2, ipph, idp, j, jc, k, i, idl, &
                           inc, l, lc, ik, idlj, idij, idj

  double precision      :: war, wai

  idot = ido/2
  nt = ip*idl1
  ipp2 = ip + 2
  ipph = (ip + 1)/2
  idp = ip*ido

  if (ido .ge. l1) then

    do j = 2, ipph
      jc = ipp2 - j
      do k = 1, l1
        do i = 1, ido
          ch(i, k, j) = cc(i, j, k) + cc(i, jc, k)
          ch(i, k, jc) = cc(i, j, k) - cc(i, jc, k)
        end do
      end do
    end do

    do k = 1, l1
      do i = 1, ido
        ch(i, k, 1) = cc(i, 1, k)
      end do
    end do

  else
  
    do j = 2, ipph
      jc = ipp2 - j
      do i = 1, ido
        do k = 1, l1
          ch(i, k, j) = cc(i, j, k) + cc(i, jc, k)
          ch(i, k, jc) = cc(i, j, k) - cc(i, jc, k)
        end do
      end do
    end do

    do i = 1, ido
      do k = 1, l1
        ch(i, k, 1) = cc(i, 1, k)
      end do
    end do

  end if
  
  idl = 2 - ido
  inc = 0

  do l = 2, ipph
    lc = ipp2 - l
    idl = idl + ido
    do ik = 1, idl1
      c2(ik, l) = ch2(ik, 1) + wa(idl - 1)*ch2(ik, 2)
      c2(ik, lc) = wa(idl)*ch2(ik, ip)
    end do
    idlj = idl
    inc = inc + ido
    do j = 3, ipph
      jc = ipp2 - j
      idlj = idlj + inc
      if (idlj .gt. idp) idlj = idlj - idp
      war = wa(idlj - 1)
      wai = wa(idlj)
      do ik = 1, idl1
        c2(ik, l) = c2(ik, l) + war*ch2(ik, j)
        c2(ik, lc) = c2(ik, lc) + wai*ch2(ik, jc)
      end do
    end do
  end do

  do j = 2, ipph
    do ik = 1, idl1
      ch2(ik, 1) = ch2(ik, 1) + ch2(ik, j)
    end do
  end do

  do j = 2, ipph
    jc = ipp2 - j
    do ik = 2, idl1, 2
      ch2(ik - 1, j) = c2(ik - 1, j) - c2(ik, jc)
      ch2(ik - 1, jc) = c2(ik - 1, j) + c2(ik, jc)
      ch2(ik, j) = c2(ik, j) + c2(ik - 1, jc)
      ch2(ik, jc) = c2(ik, j) - c2(ik - 1, jc)
    end do
  end do
  nac = 1

  if (ido .ne. 2) then

    nac = 0
    do ik = 1, idl1
      c2(ik, 1) = ch2(ik, 1)
    end do
    do j = 2, ip
      do k = 1, l1
        c1(1, k, j) = ch(1, k, j)
        c1(2, k, j) = ch(2, k, j)
      end do
    end do
  
    if (idot .le. l1) then
  
      idij = 0
      do j = 2, ip
        idij = idij + 2
        do i = 4, ido, 2
          idij = idij + 2
          do k = 1, l1
            c1(i - 1, k, j) = wa(idij - 1)*ch(i - 1, k, j) - &
                              wa(idij)*ch(i, k, j)
            c1(i, k, j) = wa(idij - 1)*ch(i, k, j) + wa(idij)*ch(i - 1, k, j)
          end do
        end do
      end do
  
    else
    
      idj = 2 - ido
      do j = 2, ip
        idj = idj + ido
        do k = 1, l1
          idij = idj
          do i = 4, ido, 2
            idij = idij + 2
            c1(i - 1, k, j) = wa(idij - 1)*ch(i - 1, k, j) - &
                              wa(idij)*ch(i, k, j)
            c1(i, k, j) = wa(idij - 1)*ch(i, k, j) + wa(idij)*ch(i - 1, k, j)
          end do
        end do
      end do
  
    end if

  end if

  return

end subroutine passb

end subroutine cfftb1

!*******************************************************************************
! subroutine cfftf(n,c,wsave,ifac)
! 
! subroutine cfftf computes the forward complex discrete fourier
! transform (the fourier analysis). equivalently , cfftf computes
! the fourier coefficients of a complex periodic sequence.
! the transform is defined below at output parameter c.
!
! the transform is not normalized. to obtain a normalized transform
! the output must be divided by n. otherwise a call of cfftf
! followed by a call of cfftb will multiply the sequence by n.
!
! the array wsave which is used by subroutine cfftf must be
! initialized by calling subroutine cffti(n,wsave,ifac).
!
! input parameters
!
! n      the length of the complex sequence c. the method is
!       more efficient when n is the product of small primes. n
!
! c      a complex array of length n which contains the sequence
!
! wsave   a real work array which must be dimensioned at least 4n
!         in the program that calls cfftf. the wsave array must be
!         initialized by calling subroutine cffti(n,wsave) and a
!         different wsave array must be used for each different
!         value of n. this initialization does not have to be
!         repeated so long as n remains unchanged thus subsequent
!         transforms can be obtained faster than the first.
!         the same wsave array can be used by cfftf and cfftb.
!
! ifac    The saved factors array, dimensioned 30.  This was originally
!         part of wsave, but has been separated to make it possible to
!         put the code into a module.  It was created by cffti().
! 
! output parameters
!
! c      for j=1,...,n
!
!            c(j)=the sum from k=1,...,n of
!
!                  c(k)*exp(-i*(j-1)*(k-1)*2*pi/n)
!
!                        where i=sqrt(-1)
!
! wsave   contains initialization calculations which must not be
!         destroyed between calls of subroutine cfftf or cfftb
!
! ifac    contains initialization calculations which must not be
!         destroyed between calls of subroutine cfftf or cfftb
!
!*******************************************************************************

subroutine cfftf(n, c, wsave, ifac)
      
  implicit none

! Formal arguments:

  integer           :: n
  double precision  :: c(*)
  double precision  :: wsave(*)
  integer           :: ifac(*)

  if (n .ne. 1) call cfftf1(n, c, wsave, wsave(2 * n + 1), ifac)

  return

end subroutine cfftf

!*******************************************************************************
!
! Subroutine:  cfftf1
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine cfftf1(n, c, ch, wa, ifac)
      
  implicit none

! Formal arguments:

  integer               :: n
  double precision      :: c(*)
  double precision      :: ch(*)
  double precision      :: wa(*)
  integer               :: ifac(*)

! Local variables:

  integer       :: nf, na, l1, iw, k1, ip, l2, ido, idot, idl1, &
                   ix2, ix3, ix4, nac, n2
      
  nf = ifac(2)
  na = 0
  l1 = 1
  iw = 1

  do k1 = 1, nf

    ip = ifac(k1 + 2)
    l2 = ip*l1
    ido = n/l2
    idot = ido + ido
    idl1 = idot*l1
    
    if (ip .eq. 4) then

      ix2 = iw + idot
      ix3 = ix2 + idot

      if (na .ne. 0) then
        call passf4(idot, l1, ch, c, wa(iw), wa(ix2), wa(ix3))
      else
        call passf4(idot, l1, c, ch, wa(iw), wa(ix2), wa(ix3))
      end if

      na = 1 - na

    else if (ip .eq. 2) then
         
      if (na .ne. 0) then
        call passf2(idot, l1, ch, c, wa(iw))
      else
        call passf2(idot, l1, c, ch, wa(iw))
      end if

      na = 1 - na

    else if (ip .eq. 3) then

      ix2 = iw + idot

      if (na .ne. 0) then
        call passf3(idot, l1, ch, c, wa(iw), wa(ix2))
      else
        call passf3(idot, l1, c, ch, wa(iw), wa(ix2))
      end if
  
      na = 1 - na
  
    else if (ip .eq. 5) then

      ix2 = iw + idot
      ix3 = ix2 + idot
      ix4 = ix3 + idot

      if (na .ne. 0) then
        call passf5(idot, l1, ch, c, wa(iw), wa(ix2), wa(ix3), wa(ix4))
      else
        call passf5(idot, l1, c, ch, wa(iw), wa(ix2), wa(ix3), wa(ix4))
      end if

      na = 1 - na

    else

      if (na .ne. 0) then
        call passf(nac, idot, ip, l1, idl1, ch, ch, ch, c, c, wa(iw))
      else
        call passf(nac, idot, ip, l1, idl1, c, c, c, ch, ch, wa(iw))
      end if

      if (nac .ne. 0) na = 1 - na

    end if

    l1 = l2
    iw = iw + (ip - 1)*idot

  end do

  if (na .eq. 0) return

  n2 = n + n

  c(1:n2) = ch(1:n2)

  return

contains

!*******************************************************************************
!
! Internal Subroutine:  passf4
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine passf4(ido, l1, cc, ch, wa1, wa2, wa3)
  
  implicit none

! Formal arguments:

  integer               :: ido
  integer               :: l1
  double precision      :: cc(ido, 4, l1)           
  double precision      :: ch(ido, l1, 4)           
  double precision      :: wa1(*)
  double precision      :: wa2(*)
  double precision      :: wa3(*)

! Local variables:

  integer               :: i, k
  double precision      :: ci2, ci3, ci4
  double precision      :: cr2, cr3, cr4
  double precision      :: ti1, ti2, ti3, ti4
  double precision      :: tr1, tr2, tr3, tr4

  if (ido .eq. 2) then

    do k = 1, l1
      ti1 = cc(2, 1, k) - cc(2, 3, k)
      ti2 = cc(2, 1, k) + cc(2, 3, k)
      tr4 = cc(2, 2, k) - cc(2, 4, k)
      ti3 = cc(2, 2, k) + cc(2, 4, k)
      tr1 = cc(1, 1, k) - cc(1, 3, k)
      tr2 = cc(1, 1, k) + cc(1, 3, k)
      ti4 = cc(1, 4, k) - cc(1, 2, k)
      tr3 = cc(1, 2, k) + cc(1, 4, k)
      ch(1, k, 1) = tr2 + tr3
      ch(1, k, 3) = tr2 - tr3
      ch(2, k, 1) = ti2 + ti3
      ch(2, k, 3) = ti2 - ti3
      ch(1, k, 2) = tr1 + tr4
      ch(1, k, 4) = tr1 - tr4
      ch(2, k, 2) = ti1 + ti4
      ch(2, k, 4) = ti1 - ti4
    end do

  else
  
    do k = 1, l1
      do i = 2, ido, 2
        ti1 = cc(i, 1, k) - cc(i, 3, k)
        ti2 = cc(i, 1, k) + cc(i, 3, k)
        ti3 = cc(i, 2, k) + cc(i, 4, k)
        tr4 = cc(i, 2, k) - cc(i, 4, k)
        tr1 = cc(i - 1, 1, k) - cc(i - 1, 3, k)
        tr2 = cc(i - 1, 1, k) + cc(i - 1, 3, k)
        ti4 = cc(i - 1, 4, k) - cc(i - 1, 2, k)
        tr3 = cc(i - 1, 2, k) + cc(i - 1, 4, k)
        ch(i - 1, k, 1) = tr2 + tr3
        cr3 = tr2 - tr3
        ch(i, k, 1) = ti2 + ti3
        ci3 = ti2 - ti3
        cr2 = tr1 + tr4
        cr4 = tr1 - tr4
        ci2 = ti1 + ti4
        ci4 = ti1 - ti4
        ch(i - 1, k, 2) = wa1(i - 1)*cr2 + wa1(i)*ci2
        ch(i, k, 2) = wa1(i - 1)*ci2 - wa1(i)*cr2
        ch(i - 1, k, 3) = wa2(i - 1)*cr3 + wa2(i)*ci3
        ch(i, k, 3) = wa2(i - 1)*ci3 - wa2(i)*cr3
        ch(i - 1, k, 4) = wa3(i - 1)*cr4 + wa3(i)*ci4
        ch(i, k, 4) = wa3(i - 1)*ci4 - wa3(i)*cr4
      end do
    end do

  end if

  return
  
end subroutine passf4

!*******************************************************************************
!
! Internal Subroutine:  passf3
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine passf3(ido, l1, cc, ch, wa1, wa2)
  
  implicit none

! Formal arguments:

  integer               :: ido
  integer               :: l1
  double precision      :: cc(ido, 3, l1)           
  double precision      :: ch(ido, l1, 3)           
  double precision      :: wa1(*)
  double precision      :: wa2(*)

! Local variables:

  integer               :: i, k

  double precision      :: taur, taui, tr2, cr2, ti2, ci2, cr3, ci3, &
                           dr2, dr3, di2, di3

  data taur, taui / - .5d0, - .866025403784439d0/

  if (ido .eq. 2) then

    do k = 1, l1
      tr2 = cc(1, 2, k) + cc(1, 3, k)
      cr2 = cc(1, 1, k) + taur*tr2
      ch(1, k, 1) = cc(1, 1, k) + tr2
      ti2 = cc(2, 2, k) + cc(2, 3, k)
      ci2 = cc(2, 1, k) + taur*ti2
      ch(2, k, 1) = cc(2, 1, k) + ti2
      cr3 = taui*(cc(1, 2, k) - cc(1, 3, k))
      ci3 = taui*(cc(2, 2, k) - cc(2, 3, k))
      ch(1, k, 2) = cr2 - ci3
      ch(1, k, 3) = cr2 + ci3
      ch(2, k, 2) = ci2 + cr3
      ch(2, k, 3) = ci2 - cr3
    end do

  else
  
    do k = 1, l1
      do i = 2, ido, 2
        tr2 = cc(i - 1, 2, k) + cc(i - 1, 3, k)
        cr2 = cc(i - 1, 1, k) + taur*tr2
        ch(i - 1, k, 1) = cc(i - 1, 1, k) + tr2
        ti2 = cc(i, 2, k) + cc(i, 3, k)
        ci2 = cc(i, 1, k) + taur*ti2
        ch(i, k, 1) = cc(i, 1, k) + ti2
        cr3 = taui*(cc(i - 1, 2, k) - cc(i - 1, 3, k))
        ci3 = taui*(cc(i, 2, k) - cc(i, 3, k))
        dr2 = cr2 - ci3
        dr3 = cr2 + ci3
        di2 = ci2 + cr3
        di3 = ci2 - cr3
        ch(i, k, 2) = wa1(i - 1)*di2 - wa1(i)*dr2
        ch(i - 1, k, 2) = wa1(i - 1)*dr2 + wa1(i)*di2
        ch(i, k, 3) = wa2(i - 1)*di3 - wa2(i)*dr3
        ch(i - 1, k, 3) = wa2(i - 1)*dr3 + wa2(i)*di3
      end do
    end do

  end if

  return

end subroutine passf3

!*******************************************************************************
!
! Internal Subroutine:  passf2
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine passf2(ido, l1, cc, ch, wa1)
  
  implicit none

! Formal arguments:

  integer               :: ido
  integer               :: l1
  double precision      :: cc(ido, 2, l1)           
  double precision      :: ch(ido, l1, 2)           
  double precision      :: wa1(*)

! Local variables:

  integer               :: i, k
  double precision      :: tr2, ti2

  if (ido .le. 2) then

    do k = 1, l1
      ch(1, k, 1) = cc(1, 1, k) + cc(1, 2, k)
      ch(1, k, 2) = cc(1, 1, k) - cc(1, 2, k)
      ch(2, k, 1) = cc(2, 1, k) + cc(2, 2, k)
      ch(2, k, 2) = cc(2, 1, k) - cc(2, 2, k)
    end do
  
  else
  
    do k = 1, l1
      do i = 2, ido, 2
        ch(i - 1, k, 1) = cc(i - 1, 1, k) + cc(i - 1, 2, k)
        tr2 = cc(i - 1, 1, k) - cc(i - 1, 2, k)
        ch(i, k, 1) = cc(i, 1, k) + cc(i, 2, k)
        ti2 = cc(i, 1, k) - cc(i, 2, k)
        ch(i, k, 2) = wa1(i - 1)*ti2 - wa1(i)*tr2
        ch(i - 1, k, 2) = wa1(i - 1)*tr2 + wa1(i)*ti2
      end do
    end do

  end if

  return

end subroutine passf2

!*******************************************************************************
!
! Internal Subroutine:  passf5
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine passf5(ido, l1, cc, ch, wa1, wa2, wa3, wa4)
  
  implicit none

! Formal arguments:

  integer               :: ido
  integer               :: l1
  double precision      :: cc(ido, 5, l1)           
  double precision      :: ch(ido, l1, 5)           
  double precision      :: wa1(*)
  double precision      :: wa2(*)
  double precision      :: wa3(*)
  double precision      :: wa4(*)

! Local variables:

  integer               :: i, k
  double precision      :: ci2, ci3, ci4, ci5
  double precision      :: cr2, cr3, cr4, cr5
  double precision      :: di2, di3, di4, di5
  double precision      :: dr2, dr3, dr4, dr5
  double precision      :: ti2, ti3, ti4, ti5
  double precision      :: tr2, tr3, tr4, tr5
  double precision      :: tr11, ti11, tr12, ti12
  
  data tr11, ti11, tr12, ti12 /.309016994374947d0, &
                             - .951056516295154d0, &
                             - .809016994374947d0, &
                             - .587785252292473d0/

  if (ido .eq. 2) then

    do k = 1, l1
      ti5 = cc(2, 2, k) - cc(2, 5, k)
      ti2 = cc(2, 2, k) + cc(2, 5, k)
      ti4 = cc(2, 3, k) - cc(2, 4, k)
      ti3 = cc(2, 3, k) + cc(2, 4, k)
      tr5 = cc(1, 2, k) - cc(1, 5, k)
      tr2 = cc(1, 2, k) + cc(1, 5, k)
      tr4 = cc(1, 3, k) - cc(1, 4, k)
      tr3 = cc(1, 3, k) + cc(1, 4, k)
      ch(1, k, 1) = cc(1, 1, k) + tr2 + tr3
      ch(2, k, 1) = cc(2, 1, k) + ti2 + ti3
      cr2 = cc(1, 1, k) + tr11*tr2 + tr12*tr3
      ci2 = cc(2, 1, k) + tr11*ti2 + tr12*ti3
      cr3 = cc(1, 1, k) + tr12*tr2 + tr11*tr3
      ci3 = cc(2, 1, k) + tr12*ti2 + tr11*ti3
      cr5 = ti11*tr5 + ti12*tr4
      ci5 = ti11*ti5 + ti12*ti4
      cr4 = ti12*tr5 - ti11*tr4
      ci4 = ti12*ti5 - ti11*ti4
      ch(1, k, 2) = cr2 - ci5
      ch(1, k, 5) = cr2 + ci5
      ch(2, k, 2) = ci2 + cr5
      ch(2, k, 3) = ci3 + cr4
      ch(1, k, 3) = cr3 - ci4
      ch(1, k, 4) = cr3 + ci4
      ch(2, k, 4) = ci3 - cr4
      ch(2, k, 5) = ci2 - cr5
    end do

  else
  
    do k = 1, l1
      do i = 2, ido, 2
        ti5 = cc(i, 2, k) - cc(i, 5, k)
        ti2 = cc(i, 2, k) + cc(i, 5, k)
        ti4 = cc(i, 3, k) - cc(i, 4, k)
        ti3 = cc(i, 3, k) + cc(i, 4, k)
        tr5 = cc(i - 1, 2, k) - cc(i - 1, 5, k)
        tr2 = cc(i - 1, 2, k) + cc(i - 1, 5, k)
        tr4 = cc(i - 1, 3, k) - cc(i - 1, 4, k)
        tr3 = cc(i - 1, 3, k) + cc(i - 1, 4, k)
        ch(i - 1, k, 1) = cc(i - 1, 1, k) + tr2 + tr3
        ch(i, k, 1) = cc(i, 1, k) + ti2 + ti3
        cr2 = cc(i - 1, 1, k) + tr11*tr2 + tr12*tr3
        ci2 = cc(i, 1, k) + tr11*ti2 + tr12*ti3
        cr3 = cc(i - 1, 1, k) + tr12*tr2 + tr11*tr3
        ci3 = cc(i, 1, k) + tr12*ti2 + tr11*ti3
        cr5 = ti11*tr5 + ti12*tr4
        ci5 = ti11*ti5 + ti12*ti4
        cr4 = ti12*tr5 - ti11*tr4
        ci4 = ti12*ti5 - ti11*ti4
        dr3 = cr3 - ci4
        dr4 = cr3 + ci4
        di3 = ci3 + cr4
        di4 = ci3 - cr4
        dr5 = cr2 + ci5
        dr2 = cr2 - ci5
        di5 = ci2 - cr5
        di2 = ci2 + cr5
        ch(i - 1, k, 2) = wa1(i - 1)*dr2 + wa1(i)*di2
        ch(i, k, 2) = wa1(i - 1)*di2 - wa1(i)*dr2
        ch(i - 1, k, 3) = wa2(i - 1)*dr3 + wa2(i)*di3
        ch(i, k, 3) = wa2(i - 1)*di3 - wa2(i)*dr3
        ch(i - 1, k, 4) = wa3(i - 1)*dr4 + wa3(i)*di4
        ch(i, k, 4) = wa3(i - 1)*di4 - wa3(i)*dr4
        ch(i - 1, k, 5) = wa4(i - 1)*dr5 + wa4(i)*di5
        ch(i, k, 5) = wa4(i - 1)*di5 - wa4(i)*dr5
      end do
    end do

  end if

  return

end subroutine passf5

!*******************************************************************************
!
! Internal Subroutine:  passf
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine passf(nac, ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa)
      
  implicit none

! Formal arguments:

  integer               :: nac
  integer               :: ido
  integer               :: ip
  integer               :: l1
  integer               :: idl1
  double precision      :: cc(ido, ip, l1)          
  double precision      :: c1(ido, l1, ip)          
  double precision      :: c2(idl1, ip)
  double precision      :: ch(ido, l1, ip)          
  double precision      :: ch2(idl1, ip)
  double precision      :: wa(*)

! Local variables:

  integer               :: idot, nt, ipp2, ipph, idp, j, jc, k, i, idl, inc, &
                           l, lc, ik, idlj, idij, idj

  double precision      :: war, wai
      
  idot = ido/2
  nt = ip*idl1
  ipp2 = ip + 2
  ipph = (ip + 1)/2
  idp = ip*ido

  if (ido .ge. l1) then

    do j = 2, ipph
      jc = ipp2 - j
      do k = 1, l1
        do i = 1, ido
          ch(i, k, j) = cc(i, j, k) + cc(i, jc, k)
          ch(i, k, jc) = cc(i, j, k) - cc(i, jc, k)
        end do
      end do
    end do
    do k = 1, l1
      do i = 1, ido
        ch(i, k, 1) = cc(i, 1, k)
      end do
    end do

  else
  
    do j = 2, ipph
      jc = ipp2 - j
      do i = 1, ido
        do k = 1, l1
          ch(i, k, j) = cc(i, j, k) + cc(i, jc, k)
          ch(i, k, jc) = cc(i, j, k) - cc(i, jc, k)
        end do
      end do
    end do

    do k = 1, l1
      do i = 1, ido
        ch(i, k, 1) = cc(i, 1, k)
      end do
    end do

  end if
  
  idl = 2 - ido
  inc = 0

  do l = 2, ipph
    lc = ipp2 - l
    idl = idl + ido
    do ik = 1, idl1
      c2(ik, l) = ch2(ik, 1) + wa(idl - 1)*ch2(ik, 2)
      c2(ik, lc) = - wa(idl)*ch2(ik, ip)
    end do
    idlj = idl
    inc = inc + ido
    do j = 3, ipph
      jc = ipp2 - j
      idlj = idlj + inc
      if (idlj .gt. idp) idlj = idlj - idp
      war = wa(idlj - 1)
      wai = wa(idlj)
      do ik = 1, idl1
        c2(ik, l) = c2(ik, l) + war*ch2(ik, j)
        c2(ik, lc) = c2(ik, lc) - wai*ch2(ik, jc)
      end do
    end do
  end do

  do j = 2, ipph
    do ik = 1, idl1
      ch2(ik, 1) = ch2(ik, 1) + ch2(ik, j)
    end do
  end do

  do j = 2, ipph
    jc = ipp2 - j
    do ik = 2, idl1, 2
      ch2(ik - 1, j) = c2(ik - 1, j) - c2(ik, jc)
      ch2(ik - 1, jc) = c2(ik - 1, j) + c2(ik, jc)
      ch2(ik, j) = c2(ik, j) + c2(ik - 1, jc)
      ch2(ik, jc) = c2(ik, j) - c2(ik - 1, jc)
    end do
  end do

  nac = 1

  if (ido .ne. 2) then

    nac = 0

    c2(1:idl1, 1) = ch2(1:idl1, 1)

    do j = 2, ip
      do k = 1, l1
        c1(1, k, j) = ch(1, k, j)
        c1(2, k, j) = ch(2, k, j)
      end do
    end do

    if (idot .le. l1) then

      idij = 0
      do j = 2, ip
        idij = idij + 2
        do i = 4, ido, 2
          idij = idij + 2
          do k = 1, l1
            c1(i - 1, k, j) = wa(idij - 1)*ch(i - 1, k, j) + &
                              wa(idij)*ch(i, k, j)
            c1(i, k, j) = wa(idij - 1)*ch(i, k, j) - wa(idij)*ch(i - 1, k, j)
          end do
        end do
      end do

    else
  
      idj = 2 - ido
      do j = 2, ip
        idj = idj + ido
        do k = 1, l1
          idij = idj
          do i = 4, ido, 2
            idij = idij + 2
            c1(i - 1, k, j) = wa(idij - 1)*ch(i - 1, k, j) + &
                              wa(idij)*ch(i, k, j)
            c1(i, k, j) = wa(idij - 1)*ch(i, k, j) - wa(idij)*ch(i - 1, k, j)
          end do
        end do
      end do

    end if

  end if

  return

end subroutine passf

end subroutine cfftf1

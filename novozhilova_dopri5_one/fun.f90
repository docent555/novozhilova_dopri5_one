module fun
    use, intrinsic :: iso_c_binding
    use ifcore

    integer(c_int) ne, nt, nz, neqf, neqp, lworkf, lworkp, liworkf, liworkp, nrdf, nrdp, iparf, iparp, ioutf, ioutp, ididf, ididp, itolf, itolp
    real(c_double) zex, dz, tend, dtr(2), q(3), i(2), th(2), a(2), dcir(2), r(2), f0(3), dt, &
        pitch, f10, f20, f30, rtolf, rtolp, atolf, atolp, rparf, rparp, ftol, ptol
    complex(c_double_complex) fp(2)

    integer(c_int) breaknum(3)
    real(c_double) phitmp0(3), phitmp1(3)
    complex(c_double_complex) fc, fcomp(3)

    integer(c_int), allocatable, target :: idxre(:, :), idxim(:, :), iworkf(:), iworkp(:)
    complex(c_double_complex), allocatable, target :: mean(:)
    real(c_double), allocatable, target :: tax(:), zax(:), u(:), eta(:, :), etag(:, :), w(:, :), f(:, :), w1(:, :), p(:, :), &
                                           phi(:, :), phios(:, :), wos(:, :), workf(:), workp(:)

    complex(c_double_complex), parameter :: ic = (0.0d0, 1.0d0)
    real(c_double), parameter :: pi = 2.0D0*dacos(0.0D0)

    private freq_out, zex, tend, dtr, q, i, th, a, dcir, r, f0, pitch

contains

    subroutine init()
        implicit none

        integer(c_int) ii

        call read_param()

        nt = tend/dt + 1
        nz = zex/dz + 1

        neqp = 4*ne
        nrdp = 4*ne
        lworkp = 8*neqp + 5*nrdp + 21
        liworkp = nrdp + 21

        !lenwrkp = 32*neqp
        !zstart = 0.0d0
        !errass = .false.
        !hstart = 0.0d0
        !ptol = 1.0d-7
        !method = 2
        !ifail = 0

        neqf = 6
        nrdf = 6
        lworkf = 11*neqf + 8*nrdf + 21
        liworkf = nrdf + 21

        call allocate_arrays()

        !do l = 1, neqp
        !    thres(l) = 1.0d-8
        !end do

        f(1, 1) = f10
        f(3, 1) = f20
        f(5, 1) = f30

        do ii = 1, nt
            tax(ii) = (ii - 1)*dt
        end do

        do ii = 1, nz
            zax(ii) = (ii - 1)*dz
        end do

        call calc_u(u, zex, nz, zax)

        do ii = 1, 2
            idxre(ii, :) = (/2*(ii - 1)*ne + 1:(2*ii - 1)*ne/)
            idxim(ii, :) = (/(2*ii - 1)*ne + 1:2*ii*ne/)
        end do

        !open (1, file='test.dat')
        !do ii = 1, ne
        !    write (1,'(4i)') idxre(1,ii), idxim(1,ii), idxre(2,ii), idxim(2,ii)
        !end do
        !close (1)
        !stop

    end subroutine init

    subroutine allocate_arrays()
        use, intrinsic :: iso_c_binding
        implicit none

        integer(c_int) err_alloc

        allocate (f(6, nt), p(4*ne, nz), u(nz), tax(nt), zax(nz), mean(nz), eta(2, nt), etag(2, nt), w(3, nt - 1), w1(3, nt - 1), &
                  idxre(2, ne), idxim(2, ne), workf(lworkf), iworkf(liworkf), workp(lworkp), iworkp(liworkp), &
                  wos(3, nt - 1), phi(3, nt), phios(3, nt), &
                  !pgot(neqp), ppgot(neqp), pmax(neqp), thres(neqp), workp(lenwrkp), &
                  stat=err_alloc)

        if (err_alloc /= 0) then
            print *, "allocation error"
            pause
            stop
        end if
    end subroutine allocate_arrays

    subroutine deallocate_arrays()
        use, intrinsic :: iso_c_binding
        implicit none

        integer(c_int) err_dealloc

        deallocate (f, p, u, tax, zax, mean, eta, etag, w, w1, stat=err_dealloc)

        if (err_dealloc /= 0) then
            print *, "deallocation error"
            pause
            stop
        end if
    end subroutine deallocate_arrays

    subroutine read_param() bind(c, name='read_param')
        use, intrinsic :: iso_c_binding
        import
        implicit none

        namelist /param/ ne, tend, zex, q1, q2, q3, i1, i2, th1, th2, a1, a2, dtr1, dtr2, &
            dcir1, dcir2, r1, r2, f10, f20, f30, dt, dz, pitch, ftol, ptol

        real(c_double) q1, q2, q3, i1, i2, th1, th2, a1, a2, dtr1, dtr2, dcir1, dcir2, r1, r2, tol

        open (unit=1, file='input_fortran.dat', status='old', err=101)
        read (unit=1, nml=param, err=102)
        close (unit=1)

        q(1) = q1
        q(2) = q2
        q(3) = q3
        i(1) = i1
        i(2) = i2
        th(1) = th1
        th(2) = th2
        a(1) = a1
        a(2) = a2
        dtr(1) = dtr1
        dtr(2) = dtr2
        dcir(1) = dcir1
        dcir(2) = dcir2
        r(1) = r1
        r(2) = r2

        write (*, nml=param)

        return
101     print *, "error of file open"; pause; stop
102     print *, 'error of reading file "input_fortran.dat"'; pause; stop
    end subroutine read_param

    subroutine ode4f()
        import
        implicit none

        integer(c_int) i, j, aiparf(1), aiparp(1), itp, itf
        real(c_double) :: t, z, artolf(1), aatolf(1), arparf(1), artolp(1), aatolp(1), arparp(1), xoutp, xoutf
        real(c_double) yf(6), yp(neqp), pex(neqp)
        logical(4) pressed
        character(1) key
        integer(c_int), parameter :: esc = 27
        common/internp/xoutp, itp
        common/internf/xoutf, itf

        !call d02pvf(neqp, zstart, p(:, 1), zex, ptol, thres, method, 'usual task', errass, hstart, workp, lenwrkp, ifail)

        !solve eq. at t=0
        fp(1) = f(1, 1)*exp(ic*f(2, 1))
        fp(2) = f(3, 1)*exp(ic*f(4, 1))

        rparp = 0.0
        iparp = 0
        itolp = 0
        !rtolp = 1.0d-7
        rtolp = ptol
        atolp = rtolp
        ioutp = neqp
        z = zax(1)
        xoutp = z
        itp = 0
        yp = p(:, 1)
        iworkp(:) = 0
        workp(:) = 0.0d0
        iworkp(5) = neqp

        artolp(1) = rtolp
        aatolp(1) = atolp
        arparp(1) = rparp
        aiparp(1) = iparp

        call dopri5_p(neqp, dpdz, z, yp, zex, artolp, aatolp, itolp, soloutp, ioutp, &
                      workp, lworkp, iworkp, liworkp, arparp, aiparp, ididp)

        p(:, nz) = yp(:)

        !open (1, file='test.dat')
        !do i = 1, nz
        !    write (1,'(5f17.8)') zax(i), p(1,i), p(ne+1,i), p(10,i), p(ne+10,i)
        !end do
        !close (1)
        !stop

        !do i = 1, nz - 1
        !zwant = i*dz
        !call d02pcf(dpdz, zwant, zgot, pgot, ppgot, pmax, workp, ifail)
        !    p(:, i + 1) = p(:, i)
        !end do

        eta(:, 1) = eff(p(:, nz))
        etag(:, 1) = pitch**2/(pitch**2 + 1)*eta(:, 1)

        !open (1, file='test.dat')
        !do j = 1, nz
        !    write (1, '(5f17.8)') (j-1)*dz, p(14, j), p(ne+14, j), p(128, j), p(ne+128, j)
        !end do
        !close (1)
        !stop

        rparf = 0.0
        iparf = 0
        t = tax(1)
        xoutf = t
        itf = 0
        yf = f(:, 1)
        itolf = 0
        !rtolf = 1.0d-7
        rtolf = ftol
        atolf = rtolf
        ioutf = 6

        artolf(1) = rtolf
        aatolf(1) = atolf
        arparf(1) = rparf
        aiparf(1) = iparf

        iworkf(:) = 0
        workf(:) = 0.0d0

        iworkf(5) = 6

        call dopri5_f(6, dfdt, t, yf, tend, artolf, aatolf, itolf, soloutf, ioutf, &
                      workf, lworkf, iworkf, liworkf, arparf, aiparf, ididf)

        do j = 1, neqf
            f(j, nt) = yf(j)
        end do
        call calcpex(f(:, nt), pex)
        eta(:, nt) = eff(pex)
        !eta(:, nt) = eff(p(:, nz))
        etag(:, nt) = pitch**2/(pitch**2 + 1)*eta(:, nt)

    end subroutine ode4f

    function eff(pex) result(eta)
        use, intrinsic :: iso_c_binding, only: c_double, c_int
        import, only:ne, idxre, idxim

        implicit none

        integer(c_int) i
        real(c_double) eta(2)
        real(c_double), intent(in) :: pex(:)

        do i = 1, 2
            eta(i) = 1.0d0 - sum(pex(idxre(i, :))*pex(idxre(i, :)) + pex(idxim(i, :))*pex(idxim(i, :)))/ne
        end do
    end function eff

    subroutine dpdz(neqp, z, p, prhs, rparp, iparp)
        import :: ne, zex, f, ic, dtr
        implicit none

        real(c_double) z, p(*), prhs(*)

        integer(c_int) i, reidx(ne), imidx(ne), neqp, iparp
        real(c_double) u, rparp
        complex(c_double_complex) s(ne), ptmp(ne)

        u = dexp(-3*((z - zex/2)/(zex/2))**2)

        do i = 1, 2
            ptmp = dcmplx(p(idxre(i, :)), p(idxim(i, :)))

            s = ic*(fp(i)*u - (dtr(i) + abs(ptmp)**2 - 1)*ptmp)

            prhs(idxre(i, :)) = dreal(s)
            prhs(idxim(i, :)) = dimag(s)
        end do
    end subroutine dpdz

    complex(c_double_complex) function xi(p, num)
        use, intrinsic :: iso_c_binding, only: c_int, c_double, c_double_complex
        import, only:ne, nz, mean, u, dz, idxre, idxim

        implicit none

        integer(c_int) i, num
        real(c_double) p(:, :)

        do i = 1, nz
            mean(i) = sum(dcmplx(p(idxre(num, :), i), p(idxim(num, :), i)), 1)/ne
        end do

        mean = u*mean
        !mean = dconjg(u)*mean

        xi = (0.5d0*(mean(1) + mean(2)) + sum(mean(2:nz - 1)))*dz
    end function

    subroutine dfdt(neqf, t, f, s, rparf, iparf)
        !use odep, only : dopri5
        implicit none

        integer(c_int) :: ii, jj, neqf, iparf, aiparp(1), itp
        real(c_double) t, z, f(neqf), yp(neqp), s(neqf), xoutp, &
            x1r, x1i, q31, i1, r1, th1, dcir1, cos1, sin1, &
            x2r, x2i, q32, i2, r2, th2, dcir2, cos2, sin2, q3, &
            f1, f2, f3, phi1, phi2, phi3, a1, a2, rparf, artolp(1), aatolp(1), arparp(1)
        complex(c_double_complex) x1, x2
        common/internp/xoutp, itp

        !call d02pvf(neqp, zstart, p(:, 1), zex, ptol, thres, method, 'usual task', errass, hstart, workp, lenwrk, ifail)

        fp(1) = f(1)*exp(ic*f(2))
        fp(2) = f(3)*exp(ic*f(4))

        rparp = 0.0
        iparp = 0
        itolp = 0
        !rtolp = 1.0d-7
        rtolp = ptol
        atolp = rtolp
        ioutp = neqp
        z = zax(1)
        xoutp = z
        itp = 0
        yp = p(:, 1)
        iworkp(:) = 0
        workp(:) = 0.0d0
        iworkp(5) = neqp

        artolp(1) = rtolp
        aatolp(1) = atolp
        arparp(1) = rparp
        aiparp(1) = iparp

        call dopri5_p(neqp, dpdz, z, yp, zex, artolp, aatolp, itolp, soloutp, ioutp, &
                      workp, lworkp, iworkp, liworkp, arparp, aiparp, ididp)

        p(:, nz) = yp(:)

        !open (1, file='test.dat')
        !do ii = 1, nz
        !    write (1,'(5f17.8)') zax(ii), p(1,ii), p(ne+1,ii), p(10,ii), p(ne+10,ii)
        !end do
        !close (1)
        !stop

        x1 = xi(p(1:2*ne, :), 1)
        x2 = xi(p(2*ne + 1:4*ne, :), 1)

        x1r = real(x1)
        x1i = imag(x1)
        x2r = real(x2)
        x2i = imag(x2)

        f1 = f(1)
        phi1 = f(2)
        f2 = f(3)
        phi2 = f(4)
        f3 = f(5)
        phi3 = f(6)

        q31 = q(3)/q(1)
        i1 = i(1)
        r1 = r(1)
        th1 = th(1)
        dcir1 = dcir(1)
        cos1 = cos(f(2))
        sin1 = sin(f(2))

        q32 = q(3)/q(2)
        i2 = i(2)
        r2 = r(2)
        th2 = th(2)
        dcir2 = dcir(2)
        cos2 = cos(f(4))
        sin2 = sin(f(4))

        q3 = q(3)
        a1 = a(1)
        a2 = a(2)

        s(1) = (-f1 + i1*(-x1i*cos1 + x1r*sin1) + 2*r1*f3*cos(phi3 - phi1 - th1))*q31
        s(2) = -2*dcir1*q3 + (i1/f1*(x1r*cos1 + x1i*sin1) + 2*r1*(f3/f1)*sin(phi3 - phi1 - th1))*q31

        s(3) = 0.0d0
        s(4) = 0.0d0

        !s(3) = (-f2 + i2*(-x2i*cos2 + x2r*sin2) + 2*r2*f3*cos(phi3 - phi2 - th2))*q32
        !s(4) = -2*dcir2*q3 + (i2/f2*(x2r*cos2 + x2i*sin2) + 2*r2*(f3/f2)*sin(phi3 - phi2 - th2))*q32

        s(5) = -f3 + a1*f1*cos(phi1 - phi3)
        s(6) = a1*f1/f3*sin(phi1 - phi3)

        !s(5) = -f3 + a1*f1*cos(phi1 - phi3) + a2*f2*cos(phi2 - phi3)
        !s(6) = a1*f1/f3*sin(phi1 - phi3) + a2*f2/f3*sin(phi2 - phi3)

        !if (mod(iter_num, 4) .eq. 0) then
        !    w1(1, time_num) = s(2)
        !    w1(2, time_num) = s(4)
        !    w1(3, time_num) = s(6)
        !    time_num = time_num + 1
        !end if
        !iter_num = iter_num + 1
    end subroutine dfdt

    subroutine calc_u(u, zex, nz, zax)
        import
        implicit none

        integer(c_int), intent(in) :: nz
        real(c_double), intent(in) :: zex, zax(nz)
        real(c_double), intent(out) :: u(:)

        integer(c_int) i

        do i = 1, nz
            u(i) = exp(-3*((zax(i) - zex/2)/(zex/2))**2)
        end do

    end subroutine

    subroutine soloutf(nr, xold, x, y, n, con, icomp, nd, rparf, iparf, irtrn)
        implicit none

        interface
            function contd5_f(ii, x, con, icomp, nd)
                implicit double precision(a - h, o - z)
                dimension con(5*nd), icomp(nd)
            end
        end interface

        integer(c_int) nr, n, nd, icomp(nd), iparf, irtrn, j, itf
        real(c_double) xold, x, con(5*nd), rparf, y(neqf), xoutf, pex(neqp), yy(neqf)
        logical(4) pressed
        character(1) key
        integer(c_int), parameter :: esc = 27
        common/internf/xoutf, itf

        if (nr .eq. 1) then
            itf = 1
            do j = 1, neqf
                f(j, itf) = y(j)
            end do
            call calcpex(y, pex)
            eta(:, itf) = eff(pex)
            !eta(:, itf) = eff(p(:, nz))
            etag(:, itf) = pitch**2/(pitch**2 + 1)*eta(:, itf)
            write (*, '(a,f10.5,a,f10.6,a,f10.6,a,f10.6,a,f10.6,a,f10.6,\,a)') 'Time = ', xoutf, '   |F1| = ', abs(f(1,itf)), '   |F2| = ', abs(f(3,itf)), &
                '   |F3| = ', abs(f(5, itf)), '   Eff1 = ', eta(1, itf), '   Eff2 = ', eta(2, itf), char(13)
            xoutf = x + dt
        else
10          continue
            if (x .ge. xoutf) then
                itf = itf + 1
                do j = 1, neqf
                    f(j, itf) = contd5_f(j, xoutf, con, icomp, nd)
                end do
                call calcpex(f(:, itf), pex)
                eta(:, itf) = eff(pex)
                !eta(:, itf) = eff(p(:, nz))
                etag(:, itf) = pitch**2/(pitch**2 + 1)*eta(:, itf)
                write (*, '(a,f10.5,a,f10.6,a,f10.6,a,f10.6,a,f10.6,a,f10.6,\,a)') 'Time = ', xoutf, '   |F1| = ', abs(f(1,itf)), '   |F2| = ', abs(f(3,itf)), &
                    '   |F3| = ', abs(f(5, itf)), '   Eff1 = ', eta(1, itf), '   Eff2 = ', eta(2, itf), char(13)
                xoutf = xoutf + dt
                goto 10
            end if
        end if

        pressed = peekcharqq()
        if (pressed) then
            key = getcharqq()
            if (ichar(key) .eq. esc) then
                write (*, '(/,a)') 'Quit?'
                key = getcharqq()
                if (ichar(key) .eq. 121 .or. ichar(key) .eq. 89) then
                    nt = itf
                    irtrn = -1
                    !return
                end if
            end if
        end if
        return
    end subroutine soloutf

    subroutine soloutp(nr, xold, x, y, n, con, icomp, nd, rparp, iparp, irtrn)
        implicit none

        interface
            function contd5_p(ii, x, con, icomp, nd)
                implicit double precision(a - h, o - z)
                dimension con(5*nd), icomp(nd)
            end
        end interface

        integer(c_int) nr, n, nd, icomp(nd), iparp, irtrn, j, itp
        real(c_double) xold, x, con(5*nd), rparp, y(neqp), xoutp
        logical(4) pressed
        character(1) key
        integer(c_int), parameter :: esc = 27
        common/internp/xoutp, itp

        if (nr .eq. 1) then
            itp = 1
            do j = 1, neqp
                p(j, itp) = y(j)
            end do
            xoutp = x + dz
        else
10          continue
            if (x .ge. xoutp) then
                itp = itp + 1
                do j = 1, neqp
                    p(j, itp) = contd5_p(j, xoutp, con, icomp, nd)
                end do
                xoutp = xoutp + dz
                goto 10
            end if
        end if
        return
    end subroutine soloutp

    subroutine calcpex(f, yp)
        implicit none

        real(c_double), intent(in) :: f(neqf)
        real(c_double) z, yp(neqp), artolp(1), aatolp(1), arparp(1), xoutp
        integer(c_int) i, aiparp(1), itp
        common/internp/xoutp, itp

        fp(1) = f(1)*exp(ic*f(2))
        fp(2) = f(3)*exp(ic*f(4))

        rparp = 0.0
        iparp = 0
        itolp = 0
        !rtolp = 1.0d-7
        rtolp = ptol
        atolp = rtolp
        ioutp = 0
        z = zax(1)
        xoutp = z
        itp = 0
        yp = p(:, 1)
        iworkp(:) = 0
        workp(:) = 0.0d0
        !iworkp(5) = neqp

        artolp(1) = rtolp
        aatolp(1) = atolp
        arparp(1) = rparp
        aiparp(1) = iparp

        call dopri5_p(neqp, dpdz, z, yp, zex, artolp, aatolp, itolp, solout_fiction, 0, &
                      workp, lworkp, iworkp, liworkp, arparp, aiparp, ididp)
    end subroutine calcpex

    subroutine solout_fiction
    end subroutine solout_fiction

end module fun

module log1p_custom
    use iso_c_binding
    implicit none

    integer, parameter :: dp = selected_real_kind(15, 307) ! double (64 bit)
    integer, parameter :: lp = selected_real_kind(18, 4931) ! long double (80 bit)
    integer, parameter :: qp = selected_real_kind(33, 4931) ! quadruple (128 bit)

    interface
        real(c_double) function log1pd_libmath(x) bind(c, name="log1p")
            import c_double
            real(c_double), value :: x
        end function 

        real(c_long_double) function log1pl_libmath(x) bind(c, name="log1pl")
            import c_long_double
            real(c_long_double), value :: x
        end function 

        real(c_float128) function log1pq_libquadmath(x) bind(c, name="log1pq")
            import c_float128
            real(c_float128), value :: x
        end function 
    end interface

    contains
        
        function log1p(x) result(y)
            real(dp), intent(in) :: x
            real(dp) :: y
            y = real(log1pd_libmath(real(x, c_double)), dp)
        end function

        function log1pl(x) result(y)
            real(lp), intent(in) :: x
            real(lp) :: y
            y = real(log1pl_libmath(real(x, c_long_double)), lp)
        end function

        function log1pq(x) result(y)
            real(qp), intent(in) :: x
            real(qp) :: y
            y = real(log1pq_libquadmath(real(x, c_float128)), qp)
        end function

        real(dp) function log1p_alt(x)
            real(dp), intent(in) :: x
            real(dp) :: y, z

            y = 1._dp + x
            z = y - 1._dp

            if (z == 0._dp) then
                log1p_alt = x
            else
                log1p_alt = x*log(y)/z
            end if
        end function

        real(qp) function log1pq_alt(x)
            real(qp), intent(in) :: x
            real(qp) :: y, z

            y = 1._qp + x
            z = y - 1._qp

            if (z == 0._qp) then
                log1pq_alt = x
            else
                log1pq_alt = x*log(y)/z
            end if
        end function
end module

program test_log1p
    use log1p_custom
    implicit none

    real(dp) :: x
    real(lp) :: xl
    real(qp) :: xq

    x = 1e-7_dp
    xl = 1e-7_lp
    xq = 1e-7_qp

    write(6, *) "Testing log(1+x)"
    write(6, *)

    write(6, *) "log(1+x) = x - 1/2*x**2 + 1/3*x**3 - 1/4*x**4 + ..."
    write(6, *)

    write(6, *) "Evaluate following quantities"
    write(6, *)

    write(6, *) "c1 = log(1+x)/x"
    write(6, *) "c2 = (log(1+x)-x)/x^2"
    write(6, *) "c3 = (log(1+x) - x + 1/2 x^2)/x^3"
    write(6, *) "c4 = (log(1+x) - x + 1/2 x^2 - 1/3 x^3)/x^4"
    write(6, *)

    write(6, *) "for x = ", x
    write(6, *)

    write(6, *) "note that log1p(x) = log(1+x) is defined in the standard c libmath."
    write(6, *)

    write(6, *) "1. Using fortran intrinsic log(1+x) vs. taylor expansion"
    write(6, *) "c1", log(1._dp + x )/x, 1._dp - 0.5_dp*x + 1._dp/3._dp*x**2
    write(6, *) "c2", (log(1._dp + x) - x)/x**2, -0.5_dp + 1._dp/3._dp*x - 0.25_dp*x**2
    write(6, *) "c3", (log(1._dp + x) - x + 0.5_dp*x**2)/x**3, 1._dp/3._dp - 0.25_dp*x + 0.2_dp*x**2
    write(6, *) "c4", (log(1._dp + x) - x + 0.5_dp*x**2 - 1._dp/3._dp*x**3)/x**4, -0.25_dp + 0.2_dp*x + 1._dp/6._dp*x**2
    write(6, *)

    write(6, *) "2. Using custom log(1+x) (double, 64 bit) vs. taylor expansion"
    write(6, *) "c1", log1p_alt(x)/x, 1._dp - 0.5_dp*x + 1._dp/3._dp*x**2
    write(6, *) "c2", (log1p_alt(x) - x)/x**2, -0.5_dp + 1._dp/3._dp*x - 0.25_dp*x**2
    write(6, *) "c3", (log1p_alt(x) - x + 0.5_dp*x**2)/x**3, 1._dp/3._dp - 0.25_dp*x + 0.2_dp*x**2
    write(6, *) "c4", (log1p_alt(x) - x + 0.5_dp*x**2 - 1._dp/3._dp*x**3)/x**4, -0.25_dp + 0.2_dp*x + 1._dp/6._dp*x**2
    write(6, *)

    write(6, *) "3. Using log1p (double, 64 bit) in libmath vs. taylor expansion"
    write(6, *) "c1", log1p(x)/x, 1._dp - 0.5_dp*x + 1._dp/3._dp*x**2
    write(6, *) "c2", (log1p(x) - x)/x**2, -0.5_dp + 1._dp/3._dp*x - 0.25_dp*x**2
    write(6, *) "c3", (log1p(x) - x + 0.5_dp*x**2)/x**3, 1._dp/3._dp - 0.25_dp*x + 0.2_dp*x**2
    write(6, *) "c4", (log1p(x) - x + 0.5_dp*x**2 - 1._dp/3._dp*x**3)/x**4, -0.25_dp + 0.2_dp*x + 1._dp/6._dp*x**2
    write(6, *)

    write(6, *) "4. Using log1pl (long double, 80 bit) in libmath vs. taylor expansion"
    write(6, *) "c1", log1pl(xl)/xl, 1._lp - 0.5_lp*xl + 1._lp/3._lp*xl**2
    write(6, *) "c2", (log1pl(xl) - xl)/xl**2, -0.5_lp + 1._lp/3._lp*xl - 0.25_lp*xl**2
    write(6, *) "c3", (log1pl(xl) - xl + 0.5_lp*xl**2)/xl**3, 1._lp/3._lp - 0.25_lp*xl + 0.2_lp*xl**2
    write(6, *) "c4", (log1pl(xl) - xl + 0.5_lp*xl**2 - 1._lp/3._lp*xl**3)/xl**4, -0.25_lp + 0.2_lp*xl + 1._lp/6._lp*xl**2
    write(6, *)

    write(6, *) "5. Using log1pq (quadruple, 128 bit) in libquadmath vs. taylor expansion"
    write(6, *) "c1", log1pq(xq)/xq, 1._qp - 0.5_qp*xq + 1._qp/3._qp*xq**2
    write(6, *) "c2", (log1pq(xq) - xq)/xq**2, -0.5_qp + 1._qp/3._qp*xq - 0.25_qp*xq**2
    write(6, *) "c3", (log1pq(xq) - xq + 0.5_qp*xq**2)/xq**3, 1._qp/3._qp - 0.25_qp*xq + 0.2_qp*xq**2
    write(6, *) "c4", (log1pq(xq) - xq + 0.5_qp*xq**2 - 1._qp/3._qp*xq**3)/xq**4, -0.25_qp + 0.2_qp*xq + 1._dp/6._qp*xq**2
    write(6, *)

    write(6, *) "6. Using custom log(1+x) (quadruple, 128 bit) vs. taylor expansion"
    write(6, *) "c1", log1pq_alt(xq)/xq, 1._qp - 0.5_qp*xq + 1._qp/3._qp*xq**2
    write(6, *) "c2", (log1pq_alt(xq) - xq)/xq**2, -0.5_qp + 1._qp/3._qp*xq - 0.25_qp*xq**2
    write(6, *) "c3", (log1pq_alt(xq) - xq + 0.5_qp*xq**2)/xq**3, 1._qp/3._qp - 0.25_qp*xq + 0.2_qp*xq**2
    write(6, *) "c4", (log1pq_alt(xq) - xq + 0.5_qp*xq**2 - 1._qp/3._qp*xq**3)/xq**4, -0.25_qp + 0.2_qp*xq + 1._dp/6._qp*xq**2
    write(6, *)
    
end program test_log1p

module log1p_custom
    use iso_c_binding
    implicit none

    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, parameter :: qp = selected_real_kind(33, 4931)

    interface
        real(c_float128) function log1pq_libquadmath(x) bind(c, name="log1pq")
            import c_float128
            real(c_float128) :: x
        end function 

        real(c_double) function log1p_libmath(x) bind(c, name="log1p")
            import c_double
            real(c_double) :: x
        end function 
    end interface

    contains
        
        function log1pq(x) result(y)
            real(qp), intent(in) :: x
            real(qp) :: y
            y = real(log1pq_libquadmath(real(x, c_float128)), qp)
        end function

        function log1p(x) result(y)
            real(dp), intent(in) :: x
            real(dp) :: y
            y = real(log1p_libmath(real(x, c_double)), dp)
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
end module

program test_log1p
    use log1p_custom
    implicit none

    real(dp) :: x
    real(qp) :: y

    x = 1e-7_dp
    y = 1e-7_qp

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

    write(6, *) "1. fortran intrinsic log vs. taylor expansion"
    write(6, *) "c1", log(1._dp + x )/x, 1._dp - 0.5_dp*x + 1._dp/3._dp*x**2
    write(6, *) "c2", (log(1._dp + x) - x)/x**2, -0.5_dp + 1._dp/3._dp*x - 0.25_dp*x**2
    write(6, *) "c3", (log(1._dp + x) - x + 0.5_dp*x**2)/x**3, 1._dp/3._dp - 0.25_dp*x + 0.2_dp*x**2
    write(6, *) "c4", (log(1._dp + x) - x + 0.5_dp*x**2 - 1._dp/3._dp*x**3)/x**4, -0.25_dp + 0.2_dp*x + 1._dp/6._dp*x**2
    write(6, *)

    write(6, *) "2. custom log(1+x) vs. taylor expansion"
    write(6, *) "c1", log1p_alt(x)/x, 1._dp - 0.5_dp*x + 1._dp/3._dp*x**2
    write(6, *) "c2", (log1p_alt(x) - x)/x**2, -0.5_dp + 1._dp/3._dp*x - 0.25_dp*x**2
    write(6, *) "c3", (log1p_alt(x) - x + 0.5_dp*x**2)/x**3, 1._dp/3._dp - 0.25_dp*x + 0.2_dp*x**2
    write(6, *) "c4", (log1p_alt(x) - x + 0.5_dp*x**2 - 1._dp/3._dp*x**3)/x**4, -0.25_dp + 0.2_dp*x + 1._dp/6._dp*x**2
    write(6, *)

    write(6, *) "3. doubel precision log(1+x) in libmath vs. taylor expansion"
    write(6, *) "c1", log1p(x)/x, 1._dp - 0.5_dp*x + 1._dp/3._dp*x**2
    write(6, *) "c2", (log1p(x) - x)/x**2, -0.5_dp + 1._dp/3._dp*x - 0.25_dp*x**2
    write(6, *) "c3", (log1p(x) - x + 0.5_dp*x**2)/x**3, 1._dp/3._dp - 0.25_dp*x + 0.2_dp*x**2
    write(6, *) "c4", (log1p(x) - x + 0.5_dp*x**2 - 1._dp/3._dp*x**3)/x**4, -0.25_dp + 0.2_dp*x + 1._dp/6._dp*x**2
    write(6, *)

    write(6, *) "4. quadruple precision log(1+x) in libquadmath vs. taylor expansion"
    write(6, *) "c1", log1pq(y)/y, 1._qp - 0.5_qp*y + 1._qp/3._qp*y**2
    write(6, *) "c2", (log1pq(y) - y)/y**2, -0.5_qp + 1._qp/3._qp*y - 0.25_qp*y**2
    write(6, *) "c3", (log1pq(y) - y + 0.5_qp*y**2)/y**3, 1._qp/3._qp - 0.25_qp*y + 0.2_qp*y**2
    write(6, *) "c4", (log1pq(y) - y + 0.5_qp*y**2 - 1._qp/3._qp*y**3)/y**4, -0.25_qp + 0.2_qp*y + 1._dp/6._qp*y**2
    write(6, *)

    
end program test_log1p
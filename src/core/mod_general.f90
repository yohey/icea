module mod_general
  implicit none

contains

  subroutine gauss_elimination(G, X, Msing)
    real(8), intent(inout):: G(:, :)
    real(8), intent(out):: X(:)
    integer, intent(out), optional:: Msing

    integer:: i, j, k, n, itmp(1)
    real(8), allocatable:: tmp(:)
    real(8), parameter:: bigNo = 1d25

    n = size(X)

    if (any(shape(G) /= [n, n+1])) then
       write(0, *) '[ERROR] Matrix ranks mismatch.'
       error stop
    end if

    allocate(tmp(n+1))

    X(:) = 0
    if (present(Msing)) Msing = 0

    do k = 1, n
       if (k /= n) then
          ! SEARCH FOR MAXIMUM COEFFICIENT IN EACH ROW
          do i = k, n
             tmp(i) = bigNo

             if (G(i, k) /= 0) then
                tmp(i) = maxval(abs(G(i, k+1:n+1)))

                if (tmp(i) < bigNo * abs(G(i, k))) then
                   tmp(i) = tmp(i) / abs(G(i, k))
                else
                   tmp(i) = bigNo
                end if
             end if
          end do

          if (all(tmp(k:n) == bigNo)) then
             if (present(Msing)) Msing = k
             return
          end if

          ! LOCATE ROW WITH SMALLEST MAXIMUM COEFFICIENT
          itmp(1:1) = minloc(tmp(k:n))
          i = itmp(1) + k - 1

          ! INDEX I LOCATES EQUATION TO BE USED FOR ELIMINATING THE NTH
          ! VARIABLE FROM THE REMAINING EQUATIONS
          ! INTERCHANGE EQUATIONS I AND K
          if (i /= k) then
             tmp(k:n+1) = G(i, k:n+1)
             G(i, k:n+1) = G(k, k:n+1)
             G(k, k:n+1) = tmp(k:n+1)
          end if
       end if

       if (G(k, k) == 0) then
          if (present(Msing)) Msing = k
          return
       end if

       ! DIVIDE NTH ROW BY NTH DIAGONAL ELEMENT AND ELIMINATE THE NTH
       ! VARIABLE FROM THE REMAINING EQUATIONS
       do concurrent (j = k+1:n+1)
          G(k, j) = G(k, j) / G(k, k)
       end do

       if (k /= n) then
          !DIR$ IVDEP
          do concurrent (i = k+1:n, j = k+1:n+1)
             G(i, j) = G(i, j) - G(i, k) * G(k, j)
          end do
       end if
    end do

    deallocate(tmp)

    ! BACKSOLVE FOR THE VARIABLES
    do k = n, 1, -1
       X(k) = G(k, n+1)
       if (n >= k+1) then
          X(k) = X(k) - sum(G(k, k+1:n) * X(k+1:n))
       end if
    end do

    return
  end subroutine gauss_elimination

end module mod_general

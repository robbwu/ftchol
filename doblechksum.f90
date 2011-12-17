program doublechksum
implicit none

double precision :: A(4,4) = &
   reshape( (/ 2, 1, 3, 4, 1, 2, 3, 5, 3, 3, 6, 0, 4, 5, 0, 14 /), (/ 4, 4 /) )
integer info

call MYDPOTF2('l', 4, A, 4, info)

write (*, *) 'info=', info
call pmat(A, 4, 4, 4)

contains
subroutine pmat(A, LDA, M, N)
implicit none
integer                 :: LDA, M, N
double precision        :: A(LDA, *)

integer                 :: i,j


do i=1,M
    do j = 1,N
        write (*,10, advance="no") A(i, j)
        !print *, A(i,j)
    end do
    write (*,*) " "
end do

10 format(F6.4, 2x)
end subroutine pmat

end program doublechksum

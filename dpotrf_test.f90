program  dpotrf_test
integer                 :: NB, N, i, j, k, INFO
parameter               ( NB = 3, N = 12)

double precision        :: A(N, N), B(N,N), CHK(NB-1)


! get block size for dpotrf
!NB = ILAENV( 1, 'DPOTRF', 'L', 1000, -1, -1, -1 )
call random_number(B)
do i=1,N
    do j=1,N
        if (i.lt.j) then
            B(i,j) = 0.d0
        end if
    end do
end do


print *, 'B:'
call pmat(B, N, N, N)
A = matmul(B, transpose(B))
!do i=1,N,NB
    !do j=1,N,NB
        !A(i+NB-1,j:j+NB-2) = sum(A(i:i+NB-2, j:j+NB-2), 1)
        !A(i:i+NB-2, j+NB-1) = sum( A( i:i+NB-2, j:j+NB-2), 2)
        !A(i+NB-1, j+NB-1) = sum(A(i+NB-1, j:j+NB-2))
    !end do
!end do 
!A(N-NB+1:N, 1:N) = ZERO
!A(1:N, N-NB+1:N) = ZERO
!do i = 1, N-NB, NB
    !A(N-NB+1:N, 1:N) = A(N-NB+1:N, 1:N) + A(i:i+NB-1, 1:N) 
!end do
!do j = 1, N-NB, NB
    !A(1:N, N-NB+1:N) = A(1:N, N-NB+1:N) + A(1:N, j:j+NB-1)
!end do
call RANDOM_NUMBER(CHK)
call BLDCHK(3, A, N )

print *, 'before DPOTRF'
call pmat(A, N, N, N)

call ftdpotrf('l', N, A, N, INFO)
print *, 'after DPOTRF'
print *, 'INFO', INFO
call pmat(A, N, N, N)


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
end program dpotrf_test
    


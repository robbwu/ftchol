program  dpotrf_test
implicit none
!include 'f90papi.h'
integer                 :: NB, N, i, j, k, INFO, N1, II, JJ, NB1
parameter               ( NB = 112, N = 10000/NB*NB)

double precision        :: A(N, N), B(N,N), CHK(NB-2), A1(N/NB*(NB-2), N/NB*(NB-2)), ONE,ZERO
parameter               ( ONE = 1.0D+0, ZERO = 0.0D+0 )
real                    :: mfl,T1, T2, rt, pt
integer                 ::  chkflg, ncnt
integer(kind=8)         :: fl
external ILAENV


!print *, 'ncnt ', ncnt
!NB1 = ILAENV( 1, 'DPOTRF', 'L', 1000, -1, -1, -1 )
!print *, 'NB1', NB1
N1 = N/NB*(NB-2)
! get block size for dpotrf
call random_number(B)

!print *, 'B:'
!call pmat(B, N, N, N)
!A1 = matmul(B(1:N1, 1:N1), transpose(B(1:N1,1:N1)))
call dgemm('n','t', N1, N1, N1, ONE, B, N, B, N, ZERO, A1, N1)
!call RANDOM_NUMBER(CHK)
!print *, 'chksum', CHK

DO I = 1, N, NB
   DO J = 1, N, NB
      II = (I-1)/NB*(NB-2)+1
      JJ = (J-1)/NB*(NB-2)+1
      A(I:I+NB-3, J:J+NB-3) = A1(II:II+NB-3, JJ:JJ+NB-3)
   END DO
END DO
!print *, 'A1'
!call pmat(A1, N1, N1, N1)
call CPU_TIME(T1)
call DPOTRF('l', N1, A1, N1, INFO)
call CPU_TIME(T2)
print *, 'DPOTRF time', T2-T1, 'seconds', 'with flops:', (N1/1.0D+3)**3/3/(T2-T1)
if (INFO.NE.0) print *, 'mydpotf3 info', INFO
!print *, 'factorized A1'
!call pmat(A1, N1, N1, N1)
      



!print *, 'before DPOTRF'
!call pmat(A, N, N, N)

call CPU_TIME(T1)
!call PAPIF_flops(rt, pt,fl, mfl, chkflg)
call ftdpotrf('l', N, A, N, INFO, NB)
!call PAPIF_flops(rt, pt,fl, mfl, chkflg)
call CPU_TIME(T2)
print *, 'FTDPOTRF time', T2-T1, 'seconds', 'with flops:', (N/1.0D+3)**3/3/(T2-T1), 'effective flops', (N1/1.0D+3)**3/(T2-T1)/3
!print *, 'PAPI mflops ', mfl
if (INFO.NE.0) print *, 'ftdpotrf info', INFO
!print *, 'after DPOTRF'
!print *, 'INFO', INFO
!call pmat(A, N, N, N)



DO I = 1, N, NB
   DO J = 1, N, NB
      II = (I-1)/NB*(NB-2)+1
      JJ = (J-1)/NB*(NB-2)+1
      A1(II:II+NB-3, JJ:JJ+NB-3) = A1(II:II+NB-3, JJ:JJ+NB-3) - A(I:I+NB-3, J:J+NB-3)
   END DO
END DO
print *, 'infnorm of |A-A1|= ', infnorm(A1, N1, N1)


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

10 format(F8.3, 2x)
end subroutine pmat
double precision function infnorm(A, M, N)
implicit none
integer :: M, N
double precision :: A(M, N)

infnorm = MAXVAL( ABS(A) )

end function infnorm


end program dpotrf_test
    


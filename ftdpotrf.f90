!  Authors:                                                             
!  ========                                                             
!                                                                       
!> \author Univ. of Tennessee                                           
!> \author Univ. of California Berkeley                                 
!> \author Univ. of Colorado Denver                                     
!> \author NAG Ltd.                                                     

!                                                                       
!> \date November 2011                                                  
!                                                                       
!> \ingroup doublePOcomputational                                       
!                                                                       
!
!   Add Fault tolerant functionality by
!   Panruo Wu(pwu@mines.edu)
!   December 2011
!  =====================================================================
      SUBROUTINE FTDPOTRF( UPLO, N, A, LDA, INFO , NB) 
      IMPLICIT NONE 
!                                                                       
!  -- LAPACK computational routine (version 3.4.0) --                   
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
!     November 2011                                                     
!                                                                       
!     .. Scalar Arguments ..                                            
      CHARACTER          UPLO 
      INTEGER            INFO, LDA, N 
!     ..                                                                
!     .. Array Arguments ..                                             
      DOUBLE PRECISION   A( LDA, * ) 
!     ..                                                                
!                                                                       
!  =====================================================================
!                                                                       
!     .. Parameters ..                                                  
      DOUBLE PRECISION   ONE, ZERO, EPS
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0, EPS = 1.0D-7 ) 
!     ..                                                                
!     .. Local Scalars ..                                               
      LOGICAL            UPPER 
      INTEGER            J, JB, NB, II, JJ , ERR
!     ..                                                                
!     .. External Functions ..                                          
      LOGICAL            LSAME 
      INTEGER            ILAENV 
      EXTERNAL           LSAME, ILAENV 
!     ..                                                                
!     .. External Subroutines ..                                        
      EXTERNAL           DGEMM, MYDPOTF2, DSYRK, DTRSM, XERBLA 
!     ..                                                                
!     .. Intrinsic Functions ..                                         
      INTRINSIC          MAX, MIN 

      DOUBLE PRECISION, ALLOCATABLE  ::  CHKSUM(:)
      REAL                          :: T1, T2
!     ..                                                                
!     .. Executable Statements ..                                       
!                                                                       
!     Test the input parameters.                                        
!                                                                       
      INFO = 0 
      UPPER = LSAME( UPLO, 'U' ) 
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN 
         INFO = -1 
      ELSE IF( N.LT.0 ) THEN 
         INFO = -2 
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN 
         INFO = -4 
      END IF 
      IF( INFO.NE.0 ) THEN 
         CALL XERBLA( 'DPOTRF', -INFO ) 
         RETURN 
      END IF 
!                                                                       
!     Quick return if possible                                          
!                                                                       
      IF( N.EQ.0 )                                                      &
     &   RETURN                                                         
!                                                                       
!     Determine the block size for this environment.                    
!                                                                       
      !NB = ILAENV( 1, 'DPOTRF', UPLO, N, -1, -1, -1 )                  
      !write (*,*) "bs",ILAENV( 1, 'DPOTRF', UPLO, N, -1, -1, -1 ) 
                                                                        
      !NB = 96
      CALL CPU_TIME(T1)
      ALLOCATE(CHKSUM(NB-2))
      CALL RANDOM_NUMBER(CHKSUM)
      CALL BLDCHK3(NB, A, N, CHKSUM )
      CALL CPU_TIME(T2)
      PRINT *, 'BLDCHK3 takes', T2-T1

      IF( NB.LE.1 .OR. NB.GE.N ) THEN 
!                                                                       
!        Use unblocked code.                                            
!                                                                       
         CALL MYDPOTF3( UPLO, N, A, LDA, INFO ) 
      ELSE 
!                                                                       
!        Use blocked code.                                              
!                                                                       
         IF( UPPER ) THEN 
!                                                                       
!           Compute the Cholesky factorization A = U**T*U.              
!                                                                       
            DO 10 J = 1, N, NB 
!                                                                       
!              Update and factorize the current diagonal block and test 
!              for non-positive-definiteness.                           
!                                                                       
               JB = MIN( NB, N-J+1 ) 
               CALL DSYRK( 'Upper', 'Transpose', JB, J-1, -ONE,         &
     &                     A( 1, J ), LDA, ONE, A( J, J ), LDA )        
               CALL MYDPOTF2( 'Upper', JB, A( J, J ), LDA, INFO ) 
               IF( INFO.NE.0 )                                          &
     &            GO TO 30                                              
               IF( J+JB.LE.N ) THEN 
!                                                                       
!                 Compute the current block row.                        
!                                                                       
                  CALL DGEMM( 'Transpose', 'No transpose', JB, N-J-JB+1,&
     &                        J-1, -ONE, A( 1, J ), LDA, A( 1, J+JB ),  &
     &                        LDA, ONE, A( J, J+JB ), LDA )             
                  CALL DTRSM( 'Left', 'Upper', 'Transpose', 'Non-unit', &
     &                        JB, N-J-JB+1, ONE, A( J, J ), LDA,        &
     &                        A( J, J+JB ), LDA )                       
               END IF 
   10       CONTINUE 
!                                                                       
         ELSE 
!                                                                       
!           Compute the Cholesky factorization A = L*L**T.              
!                                                                       
            DO 20 J = 1, N, NB 
!                                                                       
!              Update and factorize the current diagonal block and test 
!              for non-positive-definiteness.                           
!                                                                       
       200     JB = MIN( NB, N-J+1 ) 
               CALL DSYRK( 'Lower', 'No transpose', JB, J-1, -ONE,      &
     &                     A( J, 1 ), LDA, ONE, A( J, J ), LDA )        
               CALL MYDPOTF3( 'Lower', JB, A( J, J ), LDA, INFO ) 
                                                                        
               !CALL CHK1(NB, A, N, J, ERR)
               !IF ( ERR.EQ.-1 ) THEN
                  !GOTO 200
               !END IF
!                                                                       
!               if the right-bottom corner is (almost) zero,            
!               set it to 1. Threshold needs reviewing.                 
!               make sure the head block factorizes correctly           
!                                                                       
               IF ( ABS( A( J+JB-1, J+JB-1) ) .LT. 1.0D-7 ) THEN 
                   A( J+JB-1, J+JB-1) = ONE 
               END IF 
               IF ( ABS( A( J+JB-2, J+JB-2) ) .LT. 1.0D-7 ) THEN 
                   A( J+JB-2, J+JB-2) = ONE 
               END IF 
               IF( INFO.NE.0 )                                          &
     &            GO TO 30                                              
               IF( J+JB.LE.N ) THEN 
!                                                                       
!                 Compute the current block column.                     
!                                                                       
                  ! Here we should use ftdgemm instead..
                  CALL DGEMM( 'No transpose', 'Transpose', N-J-JB+1, JB,&
     &                        J-1, -ONE, A( J+JB, 1 ), LDA, A( J, 1 ),  &
     &                        LDA, ONE, A( J+JB, J ), LDA )             
                  CALL DTRSM( 'Right', 'Lower', 'Transpose', 'Non-unit',&
     &                        N-J-JB+1, JB, ONE, A( J, J ), LDA,        &
     &                        A( J+JB, J ), LDA )                       
                  !CALL CHK2(NB, A, N, J, ERR)
               END IF 
               IF ( A(J+JB-2, J+JB-2).EQ.ONE .AND. A(J+JB-1, J+JB-1).EQ.ONE) THEN
                  A(J+JB-2, J+JB-2) = ZERO
                  A(J+JB-1, J+JB-1) = ZERO
               END IF
               !WRITE (*,*) "ITERATION J=" , J
               !!WRITE (*,*) J 
                !DO II=1,N 
                    !DO JJ = 1,N 
                        !WRITE (*,100, ADVANCE="NO") A(II, JJ) 
                    !END DO 
                    !WRITE (*,*) " " 
                !END DO 
  !100           FORMAT(F6.4, 2X) 
   20       CONTINUE 
         END IF 
      END IF 
      GO TO 40 
!                                                                       
   30 CONTINUE 
      INFO = INFO + J - 1 
!                                                                       
   40 CONTINUE 
      DEALLOCATE(CHKSUM)
      RETURN 
!                                                                       
!     End of DPOTRF                                                     
!                                                                       
      END                                           
                                                                        
                                                                        
                                                                        
!                                                                       
!   Build local and global unit checksum matrix of A                         
!   assuming that additional memory space are allocated                 
!                                                                       
      SUBROUTINE BLDCHK(NB, A, N) 
      IMPLICIT NONE 
                                                                        
!       N must be divisible by block size NB                            
      INTEGER           NB, N 
      DOUBLE PRECISION  A(N, N) 
      !=====================                                            
                                                                        
      DOUBLE PRECISION  S, ZERO, ONE 
      INTEGER           I, J, K, L, II, JJ 
                                                                        
      PARAMETER         ( ZERO = 0.0D+0, ONE = 0.0D+0 ) 
                                                                        
      DO J = 1, N-NB, NB 
        DO I = J, N-NB, NB 
!           make sure the blocks on diagonal is symmetric               
            IF ( I.EQ.J ) THEN 
                DO II = I, I+NB-1 
                    DO JJ = I+1, J+NB-1 
                        A(II, JJ) = A(JJ, II) 
                    END DO 
                END DO 
                DO II = I, I+NB-1 
                    A(II, J+NB-1) = A(J+NB-1, II) 
                END DO 
            END IF 
            DO K = J, J+NB-2 
                S = ZERO 
                DO L = I, I+NB-2 
                    S = S + A(L, K) 
                END DO 
                A(I+NB-1, K) = S 
            END DO 
            DO K = I, I+NB-1 
                S = ZERO 
                DO L = J, J+NB-2 
                    S = S + A(K, L) 
                END DO 
                A(K, J+NB-1) = S 
            END DO 
        END DO 
      END DO 
                                                                        
!       Now build the global checksum matrix                            
      A(N-NB+1:N, 1: N-NB) = ZERO
                                                                        
      DO J = 1, N-NB, NB 
        DO I = 1, N-NB, NB 
          DO JJ = 0, NB-1 
            DO II = 0, NB-1 
              IF (I.GE.J) THEN 
                A(N-NB+1+II, J+JJ) = A(N-NB+1+II, J+JJ) + A(I+II,J+JJ) 
              ELSE 
                A(N-NB+1+II, J+JJ) = A(N-NB+1+II, J+JJ) + A(J+JJ,I+II) 
              END IF 
            END DO 
          END DO 
        END DO 
      END DO 
                                                                        
      END                                           
                                                                        
!                                                                       
!   Build local and global random checksum matrix of A                         
!   assuming that additional memory space are allocated                 
!                                                                       
      SUBROUTINE BLDCHK2(NB, A, N, CHKSUM) 
      IMPLICIT NONE 
                                                                        
!       N must be divisible by block size NB                            
      INTEGER           NB, N 
      DOUBLE PRECISION  A(N, N), CHKSUM(NB-1)
      !=====================                                            
                                                                        
      DOUBLE PRECISION  S, ZERO, ONE 
      INTEGER           I, J, K, L, II, JJ 
                                                                        
      PARAMETER         ( ZERO = 0.0D+0, ONE = 0.0D+0 ) 
                                                                        
      DO J = 1, N-NB, NB 
        DO I = J, N-NB, NB 
!           make sure the blocks on diagonal is symmetric               
            IF ( I.EQ.J ) THEN 
                DO II = I, I+NB-1 
                    DO JJ = I+1, J+NB-1 
                        A(II, JJ) = A(JJ, II) 
                    END DO 
                END DO 
                DO II = I, I+NB-1 
                    A(II, J+NB-1) = A(J+NB-1, II) 
                END DO 
            END IF 
            DO JJ = J, J+NB-2
               A(I+NB-1, JJ) = SUM( A(I:I+NB-2, JJ) * CHKSUM )
            END DO
            DO II = I, I+NB-1
               A(II, J+NB-1) = SUM( A(II, J:J+NB-2) * CHKSUM )
            END DO
            !DO K = J, J+NB-2 
                !S = ZERO 
                !DO L = I, I+NB-2 
                    !S = S + A(L, K) 
                !END DO 
                !A(I+NB-1, K) = S 
            !END DO 
            !DO K = I, I+NB-1 
                !S = ZERO 
                !DO L = J, J+NB-2 
                    !S = S + A(K, L) 
                !END DO 
                !A(K, J+NB-1) = S 
            !END DO 
             
        END DO 
      END DO 
                                                                        
!       Now build the global checksum matrix                            
      A(N-NB+1:N, 1: N-NB) = ZERO
                                                                        
      DO J = 1, N-NB, NB 
        DO I = 1, N-NB, NB 
          DO JJ = 0, NB-1 
            DO II = 0, NB-1 
              IF (I.GE.J) THEN 
                A(N-NB+1+II, J+JJ) = A(N-NB+1+II, J+JJ) + A(I+II,J+JJ) 
              ELSE 
                A(N-NB+1+II, J+JJ) = A(N-NB+1+II, J+JJ) + A(J+JJ,I+II) 
              END IF 
            END DO 
          END DO 
        END DO 
      END DO 
                                                                        
      END                                           
      ! 2 local checksums no global
      SUBROUTINE BLDCHK3(NB, A, N, CHKSUM) 
      IMPLICIT NONE
      
      DOUBLE PRECISION  A(N,N), CHKSUM(NB-2)
      INTEGER           NB, N

      INTEGER           I, J, II, JJ
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D+0, ONE = 1.0D+0 )

      DO J = 1, N, NB
         DO  I = J, N, NB
            IF (I.EQ.J) THEN
               DO JJ = J, J+NB-3
                  DO II = 1, JJ-1
                     A(II, JJ) = A(JJ, II)
                  END DO
               END DO
            END IF
            A(I+NB-2:I+NB-1, J:J+NB-1) = ZERO
            A(I:I+NB-1, J+NB-2:J+NB-1) = ZERO
            DO II = I, I+NB-3
               A(I+NB-2, J:J+NB-3) = A(I+NB-2, J:J+NB-3) + A(II, J:J+NB-3)
               A(I+NB-1, J:J+NB-3) = A(I+NB-1, J:J+NB-3) + A(II, J:J+NB-3) * CHKSUM(II-I+1)
            END DO
            DO JJ = J, J+NB-3
               A(I:I+NB-3, J+NB-2) = A(I:I+NB-3, J+NB-2) + A(I:I+NB-3, JJ)
               A(I:I+NB-3, J+NB-1) = A(I:I+NB-3, J+NB-1) + A(I:I+NB-3, JJ) * CHKSUM(JJ-J+1)
            END DO
            A(I+NB-2, J+NB-2) = SUM( A(I+NB-2, J:J+NB-3) )
            A(I+NB-1, J+NB-1) = SUM( A(I+NB-1, J:J+NB-3) * CHKSUM )
            A(I+NB-1, J+NB-2) = SUM( A(I:I+NB-3, j+NB-2) * CHKSUM )
            A(I+NB-2, J+NB-1) = SUM( A(I+NB-2, J:J+NB-3) * CHKSUM )
         END DO
      END DO



      END SUBROUTINE BLDCHK3
                                                                        
      SUBROUTINE CHK1(NB, A, N, I, INFO) 
      IMPLICIT NONE 
                                                                        
      DOUBLE PRECISION    A(N, N) 
      INTEGER             N, I, INFO, NB
                                                                        
      !=================                                                
                                                                        
      INTEGER             II, JJ 
      DOUBLE PRECISION     ZERO, EPS 
      PARAMETER            (ZERO = 0.0D+0, EPS = 1.0D-7) 
                                                                        
!     ABS(I+NB-1, I+NB-1) should be ZERO; if not then there's           
!     something wrong. Use global chkmat to recover it and signal       
!     a recomputation                                                   
      INFO = 0 

      IF ( ABS( A(I+NB-1, I+NB-1) ) > EPS ) THEN 
         INFO = -1 
         A(I:I+NB-1, I:I+NB-1) = A(N-NB+1:N, I:I+NB-1)
         DO II = I+NB,N-NB,NB
            A(I:I+NB-1, I:I+NB-1) = A(I:I+NB-1, I:I+NB-1) - A(II:II+NB-1, I:I+NB-1);
         END DO
                 
!        pass the check.                                                
      END IF 
      
                                                                        
      END SUBROUTINE CHK1                                          

      SUBROUTINE CHK2(NB, A, N, I, INFO)
      IMPLICIT NONE 

      DOUBLE PRECISION    A(N, N) 
      INTEGER             N, I, INFO, NB
                                                                        
      !=================                                                
                                                                        
      INTEGER             II, JJ 
      DOUBLE PRECISION     ZERO, EPS 
      PARAMETER            (ZERO = 0.0D+0, EPS = 1.0D-7) 

      DOUBLE PRECISION     WORK(N)

      DO II = I+NB, N
         IF ( ABS( A(II, I+NB-1) ) > EPS ) THEN
            INFO = II
         END IF
      END DO

      !DO II = I+NB, N, NB
         


      END SUBROUTINE CHK2


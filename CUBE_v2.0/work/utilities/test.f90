implicit none

integer i,j,k,kk,l,ll,nrot
real d(3)
real,dimension(3,3) :: dd,dd1,dd2,dd3,v,e,m

e(1,1)=2.345; e(1,2)=2.0; e(1,3)=0.7;
e(2,1)=2.0; e(2,2)=3.456; e(2,3)=3.0;
e(3,1)=0.7; e(3,2)=3.0; e(3,3)=4.567;

call dsyevj3(e,d,v)
call eigsrt(d,v)
print*,'d'
print*,d
dd=0; dd1=0; dd2=0; dd3=0;
dd1(1,1)=d(1); dd2(2,2)=d(2); dd3(3,3)=d(3)
print*,'v'
print*,v(:,1)
print*,v(:,2)
print*,v(:,3)

!v=matinv3(v)
!print*,'inv(v)'
!print*,v(:,1)
!print*,v(:,2)
!print*,v(:,3)
!stop

m=matmul(matmul(v,dd1),transpose(v))&
 +matmul(matmul(v,dd2),transpose(v))&
 +matmul(matmul(v,dd3),transpose(v));
print*,'m'
print*,m(:,1)
print*,m(:,2)
print*,m(:,3)
stop


call eigsrt(d,v)
print*,'d'
print*,d
print*,'v'
print*,v(:,1)
print*,v(:,2)
print*,v(:,3)



contains

      SUBROUTINE dsyevj3(A,W,Q)
      REAL A(3,3)
      REAL Q(3,3)
      REAL W(3)

      INTEGER          N
      PARAMETER        ( N = 3 )

      REAL SD, SO
      REAL S, C, T
      REAL G, H, Z, THETA
      REAL THRESH
      INTEGER          I, X, Y, R

      DO 10 X = 1, N
        Q(X,X) = 1.0D0
        DO 11, Y = 1, X-1
          Q(X, Y) = 0.0D0
          Q(Y, X) = 0.0D0
   11   CONTINUE
   10 CONTINUE

      DO 20 X = 1, N
        W(X) = A(X, X)
   20 CONTINUE

      SD = 0.0D0
      DO 30 X = 1, N
        SD = SD + ABS(W(X))
   30 CONTINUE
      SD = SD**2

      DO 40 I = 1, 50
        SO = 0.0D0
        DO 50 X = 1, N
          DO 51 Y = X+1, N
            SO = SO + ABS(A(X, Y))
   51     CONTINUE
   50   CONTINUE
        IF (SO .EQ. 0.0D0) THEN
          RETURN
        END IF

        IF (I .LT. 4) THEN
          THRESH = 0.2D0 * SO / N**2
        ELSE
          THRESH = 0.0D0
        END IF

        DO 60 X = 1, N
          DO 61 Y = X+1, N
            G = 100.0D0 * ( ABS(A(X, Y)) )
            IF ( I .GT. 4 .AND. ABS(W(X)) + G .EQ. ABS(W(X)) &
                         .AND. ABS(W(Y)) + G .EQ. ABS(W(Y)) ) THEN
              A(X, Y) = 0.0D0
            ELSE IF (ABS(A(X, Y)) .GT. THRESH) THEN
              H = W(Y) - W(X)
              IF ( ABS(H) + G .EQ. ABS(H) ) THEN
                T = A(X, Y) / H
              ELSE
                THETA = 0.5D0 * H / A(X, Y)
                IF (THETA .LT. 0.0D0) THEN
                  T = -1.0D0 / (SQRT(1.0D0 + THETA**2) - THETA)
                ELSE
                  T = 1.0D0 / (SQRT(1.0D0 + THETA**2) + THETA)
                END IF
              END IF

              C = 1.0D0 / SQRT( 1.0D0 + T**2 )
              S = T * C
              Z = T * A(X, Y)

              A(X, Y) = 0.0D0
              W(X)    = W(X) - Z
              W(Y)    = W(Y) + Z
              DO 70 R = 1, X-1
                T       = A(R, X)
                A(R, X) = C * T - S * A(R, Y)
                A(R, Y) = S * T + C * A(R, Y)
   70         CONTINUE
              DO 80, R = X+1, Y-1
                T       = A(X, R)
                A(X, R) = C * T - S * A(R, Y)
                A(R, Y) = S * T + C * A(R, Y)
   80         CONTINUE
              DO 90, R = Y+1, N
                T       = A(X, R)
                A(X, R) = C * T - S * A(Y, R)
                A(Y, R) = S * T + C * A(Y, R)
   90         CONTINUE

              DO 100, R = 1, N
                T       = Q(R, X)
                Q(R, X) = C * T - S * Q(R, Y)
                Q(R, Y) = S * T + C * Q(R, Y)
  100         CONTINUE
            END IF
   61     CONTINUE
   60   CONTINUE
   40 CONTINUE

      PRINT *, "DSYEVJ3: No convergence."

      END SUBROUTINE

      SUBROUTINE eigsrt(d,v)
      INTEGER,PARAMETER:: n=3
      REAL d(n),v(n,n)
      INTEGER i,j,k
      REAL p
      do 13 i=1,n-1
        k=i
        p=d(i)
        do 11 j=i+1,n
          if(d(j).ge.p)then
            k=j
            p=d(j)
          endif
11      continue
        if(k.ne.i)then
          d(k)=d(i)
          d(i)=p
          do 12 j=1,n
            p=v(j,i)
            v(j,i)=v(j,k)
            v(j,k)=p
12        continue
        endif
13    continue
      return
      END SUBROUTINE

			pure function matinv3(A) result(B)
		    !! Performs a direct calculation of the inverse of a 3×3 matrix.
		    real(4), intent(in) :: A(3,3)   !! Matrix
		    real(4)             :: B(3,3)   !! Inverse matrix
		    real(4)             :: detinv

		    ! Calculate the inverse determinant of the matrix
		    detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
		              - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
		              + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

		    ! Calculate the inverse of the matrix
		    B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
		    B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
		    B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
		    B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
		    B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
		    B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
		    B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
		    B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
		    B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
		  end function
end

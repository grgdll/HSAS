      INCLUDE "foq.h"
      INTEGER I,J,K,L,M
      REAL FOQ(N_W,N_S,N_C,N_N,N_A)
      REAL F(N_W,N_S,N_C)
      REAL DUM1,DUM2,DUM3

      OPEN(10,FILE="foq.dat",STATUS="OLD")
      do 100 I = 1,N_W
      do 100 J = 1,N_S
      do 100 K = 1,N_C
      do 100 L = 1,N_N
 100  read(10,*) (FOQ(I,J,K,L,M),M=1,N_A)
      CLOSE(10)
      OPEN(10,FILE="f.dat",STATUS="OLD")
      do 150 I = 1,N_W
      do 150 J = 1,N_S
      do 150 K = 1,N_C
 150  read(10,*) DUM1,DUM2,DUM3,F(I,J,K)
      CLOSE(10)

      OPEN(10,FILE="values.dat")
      do 200 I = 1,N_W
      do 200 J = 1,N_S
      do 200 K = 1,N_C
      WRITE(10,'(f9.2,f7.1,f7.2,f10.6)')WAVE(I),SUN(J),CHL(K),F(I,J,K)
      do 200 L = 1,N_N
      do 210 M = 1,N_A
 210  WRITE(10,'(f8.5,$)') FOQ(I,J,K,L,M)
 200  WRITE(10,'("")')
      CLOSE(10)

      END

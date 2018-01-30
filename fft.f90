Module FFT_Mod
  implicit none
  Integer , parameter :: DP = Selected_Real_Kind( 15 )
  Integer , parameter :: FFT_Forward = 1
  Integer , parameter :: FFT_Inverse = -1
contains

  Subroutine fcFFT( x , forback )
    !//Subroutine FFT , Cooley-Tukey , radix-2
    !// www.fcode.cn
    Real(Kind=DP) , parameter :: PI = 3.141592654_DP
    Complex(Kind=DP) , Intent(INOUT) :: x(:)
    Integer , Intent(IN) :: forback
    Integer :: n
    integer :: i , j , k , ncur , ntmp , itmp
    real(Kind=DP) :: e
    complex(Kind=DP) :: ctmp
    n = size(x)
    ncur = n
    Do
      ntmp = ncur
      e = 2.0_DP * PI / ncur
      ncur = ncur / 2
      if ( ncur < 1 ) exit
      Do j = 1 , ncur
        Do i = j , n , ntmp
          itmp = i + ncur
          ctmp = x(i) - x(itmp)
          x(i) = x(i) + x(itmp)
          x(itmp) = ctmp * exp( forback * cmplx( 0.0_DP , e*(j-1) ) )
        End Do   
      End Do
    End Do
    j = 1
    Do i = 1, n - 1
      If ( i < j ) then
        ctmp = x(j)
        x(j) = x(i)
        x(i) = ctmp
      End If
      k = n/2
      Do while( k < j )
        j = j - k
        k = k / 2
      End Do
      j = j + k
    End Do
    Return
  End Subroutine fcFFT

End Module FFT_Mod

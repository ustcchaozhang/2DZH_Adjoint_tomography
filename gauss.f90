Subroutine gauss( x , NSTEP , dist , deltat , T )

    use FFT_Mod
    use constants
    Implicit None
    complex(Kind=DP) , Intent(INOUT) :: x(NSTEP)
    integer , Intent(IN) :: NSTEP
    real(Kind=DP) , Intent(IN) :: dist
    real(Kind=DP) , Intent(IN) :: deltat
    real(Kind=DP) , Intent(IN) :: T

    integer :: i
    real log_result , alpha , cenfre
    integer p_num
    complex(Kind=DP) , allocatable :: x_pro(:),x_pro1(:)
    complex(Kind=DP) , allocatable :: gauss_win(:)

    log_result = LOG(REAL(NSTEP))/LOG(2.0)
    p_num = 2**(CEILING(log_result))
    allocate(x_pro(p_num))
    allocate(x_pro1(p_num))
    x_pro(:) = 0.0
    do i = 1,NSTEP
        x_pro(i) = x(i)
    enddo

    alpha = (85.89*dist + 18420)/(dist + 3261)
    cenfre=deltat*p_num/T+1
    allocate(gauss_win(p_num))
    do i =1,p_num/2
        gauss_win(i) = EXP(-alpha*(i-cenfre)**2/(cenfre)**2)
    enddo
		
    call fcFFT( x_pro , FFT_Forward )

    do i =1,p_num/2+1
        x_pro1(i) = x_pro(i) * gauss_win(i)
    enddo
		
    do i =2,p_num/2
        x_pro1(p_num-i+2)=conjg(x_pro1(i))
    enddo
		
    call fcFFT( x_pro1 , FFT_Inverse )

    do i = 1,NSTEP
        x(i) = x_pro1(i)
    enddo

    Return
End Subroutine gauss

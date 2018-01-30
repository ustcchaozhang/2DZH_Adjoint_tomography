program adj_seismogram_AH
use FFT_Mod
use constants
use m_hilbert_transform
! This program cuts a certain portion of the seismograms and convert it
! into the adjoint source for generating banana-dougnut kernels

implicit none

integer, parameter :: NSTEP = 25000
integer, parameter :: nrec = 1
integer, parameter :: NT = 10000    ! labeled as point number ,T*5/0.02 
double precision, parameter :: T = 40.0    ! labeled as zhouqi
double precision, parameter :: t0 = 0.0    ! labeled as 'time delay'
double precision, parameter :: dist = 700     ! labeled as distance
double precision, parameter :: deltat = 0.02
double precision, parameter :: EPS = 1.d-40

integer :: itime,icomp,istart,iend,nlen,irec,NDIM,NDIMr,adj_comp
double precision :: time,tstart(nrec),tend(nrec)
character(len=150), dimension(nrec) :: station_name
double precision, dimension(NSTEP) :: time_window,time_window_cut,time_window_cut2
double precision :: seism(NSTEP,3),Nnorm,seism_win(NSTEP),gauss_win(NSTEP),Nnorm_x(nrec),Nnorm_z(nrec)
double precision :: seism_veloc(NSTEP),seism_accel(NSTEP),ft_bar(NSTEP)
character(len=3) :: compr(2),comp(3)
character(len=150) :: filename
real seism_win_temp(NSTEP),seism_win_h(NSTEP),seism_win_f(NSTEP),seism_win_fg(NSTEP)
complex(Kind=DP) :: seism_win_cut(NSTEP)
real,allocatable :: taper(:)
real :: taper2(NT)
integer max_center    ! the amplitude max point
real :: max_value,cenfre     ! the max point value

NDIM=3
comp = (/"BXX","BXY","BXZ"/)

! number of components
NDIMr=2  ! P-SV
! NDIMr=1  ! SH (membrane)
compr = (/"BXX","BXZ"/)    ! P-SV
! compr = (/"BXY","tmp"/)  ! SH (membrane)

! list of stations
station_name =(/"S0001"/)

! chose the component for the adjoint source (adj_comp = 1:X, 2:Y, 3:Z)
adj_comp = 1

  ! KEY: 'absolute' time interval for window
  tstart =(/700/)/4.0
  tend = (/700/)/2.0

do irec =1,nrec
    do icomp = 1, NDIMr
        filename = 'OUTPUT_FILES/'//'AA.'//trim(station_name(irec))//'.'// compr(icomp) // '.semd'
        open(unit = 10, file = trim(filename))
        do itime = 1,NSTEP
            read(10,*) time , seism(itime,icomp)
        enddo
    enddo

    if(NDIMr==2)then
        seism(:,3) = seism(:,2)
        seism(:,2) = 0.d0
    else
        seism(:,2) = seism(:,1)
        seism(:,1) = 0.d0
        seism(:,3) = 0.d0
    endif

    close(10)

    istart = max(floor(tstart(irec)/deltat),1)
    iend = min(floor(tend(irec)/deltat),NSTEP)
    print *,'istart =',istart, 'iend =', iend
    print *,'tstart =',istart*deltat, 'tend =', iend*deltat
    if(istart >= iend) stop 'check istart,iend'
    nlen = iend - istart +1

    allocate(taper(nlen))
    call window_taper(nlen,0.3,'hann',taper(:))

    do icomp = 1, NDIM
        print *,comp(icomp)

        filename = 'SEM/'//'AA.'//trim(station_name(irec))//'.'// comp(icomp) // '.adj'
        open(unit = 11, file = trim(filename))      ! the file contain adjoint source

        filename = 'SEM/'//'AA.'//trim(station_name(irec))//'.'// comp(icomp) // '.wae'
        open(unit = 12, file = trim(filename))      ! the original waveform

        filename = 'SEM/'//'AA.'//trim(station_name(irec))//'.'// comp(icomp) // '.tap'
        open(unit = 13, file = trim(filename))      ! the taper used in the waveform

        filename = 'SEM/'//'AA.'//trim(station_name(irec))//'.'// comp(icomp) // '.tpd'
        open(unit = 14, file = trim(filename))      ! the tapered waveform

        filename = 'SEM/'//'AA.'//trim(station_name(irec))//'.'// comp(icomp) // '.ifr'
        open(unit = 16, file = trim(filename))      ! the inverse Fourier result

        filename = 'SEM/'//'AA.'//trim(station_name(irec))//'.'// comp(icomp) // '.env'
        open(unit = 17, file = trim(filename))      ! the gaussed waveform envelope

        filename = 'SEM/'//'AA.'//trim(station_name(irec))//'.'// comp(icomp) // '.tp2'
        open(unit = 18, file = trim(filename))      ! the taper used in the gaussed waveform

        filename = 'SEM/'//'AA.'//trim(station_name(irec))//'.'// comp(icomp) // '.tp2d'
        open(unit = 19, file = trim(filename))      ! the gaussed and tapered waveform

        filename = 'ZHratio.txt'
        open(unit = 20, file = trim(filename))      ! the zh ratio output

        time_window(:) = 0.d0
        time_window_cut(:) = 0.d0
        time_window_cut2(:) = 0.d0
        seism_win(:) = seism(:,icomp)
        seism_veloc(:) = 0.d0
        seism_accel(:) = 0.d0

        do itime =1,NSTEP      ! output the waveform
            write(12,*) (itime-1)*deltat - t0, seism_win(itime)
        enddo

        do itime =1,nlen      ! creat a hamming time window used in the waveform
            time_window_cut(itime + istart) = taper(itime)
        enddo

        do itime =1,NSTEP      ! output the hamming time window
            write(13,*) (itime-1)*deltat - t0, time_window_cut(itime)
        enddo

        seism_win_cut(:) = seism_win(:) * time_window_cut(:)   ! cut the waveform use the hamming window

        do itime =1,NSTEP      ! output the taped waveform
            write(14,*) (itime-1)*deltat - t0, real(seism_win_cut(itime))
        enddo

        call gauss(seism_win_cut,NSTEP,dist,deltat,T)

        seism_win_f(:) = real(seism_win_cut(:))
        seism_win_temp(:) = seism_win_f(:)

        do itime =1,NSTEP      ! print the gaussed waveform
            write(16,*) (itime-1)*deltat - t0, seism_win_f(itime)
        enddo

        call hilbert(seism_win_temp,NSTEP)    ! hilbert transform
        seism_win_h(:)=sqrt(seism_win_f(:)**2+seism_win_temp(:)**2)

        do itime =1,NSTEP      ! print the envelope of gaussed waveform
            write(17,*) (itime-1)*deltat - t0, seism_win_h(itime)
        enddo

        max_value = 0.0
        do itime =1,NSTEP      ! find the biggest value
            if (seism_win_h(itime) .GT. max_value) then
                max_value = seism_win_h(itime)
                max_center = itime
            endif
        enddo

        print *,'max_center =', max_center

        call window_taper(NT,0.3,'cose',taper2)
        do itime =1,NT      ! add a hamming time window for the second time
            time_window_cut2(itime + max_center - NT/2) = taper2(itime)
        enddo

        do itime =1,NSTEP      ! print the hanning time window for the second time
            write(18,*) (itime-1)*deltat - t0, time_window_cut2(itime)
        enddo

        do itime =1,NSTEP      ! apply hamming to the gaussed waveform
        !    seism_win_fg(itime) = seism_win_f(itime) * time_window_cut2(itime)
            seism_win_fg(itime) = seism_win_f(itime) * 1.0
        enddo

        do itime =1,NSTEP      ! print the tapered gaussed waveform
            write(19,*) (itime-1)*deltat - t0, seism_win_fg(itime)
        enddo

        !do itime =istart,iend  	   ! creat welch/cosine time window
            !time_window(itime) = 1.d0 - cos(pi*(itime-1)/NSTEP+1)**10   ! cosine window
            !time_window(itime) = 1.d0 - (2* (dble(itime) - istart)/(iend-istart) -1.d0)**2  ! Welch window
        !enddo

        !do itime = 2,NSTEP-1      !calculate derivation , velocity
            !seism_veloc(itime) = (seism_win_fg(itime+1) - seism_win_fg(itime-1))/(2*deltat)
        !enddo
        !seism_veloc(1) = (seism_win_fg(2) - seism_win_fg(1))/deltat
        !seism_veloc(NSTEP) = (seism_win_fg(NSTEP) - seism_win_fg(NSTEP-1))/deltat

        !do itime = 2,NSTEP-1
            !seism_accel(itime) = (seism_veloc(itime+1) - seism_veloc(itime-1))/(2*deltat)
        !enddo
        !seism_accel(1) = (seism_veloc(2) - seism_veloc(1))/deltat
        !seism_accel(NSTEP) = (seism_veloc(NSTEP) - seism_veloc(NSTEP-1))/deltat

        !do itime =1,NSTEP      ! print the accel series
            !write(20,*) (itime-1)*deltat - t0, seism_accel(itime)
        !enddo

        ! cross-correlation traveltime adjoint source
        Nnorm = deltat * sum(seism_win_fg(:) * seism_win_fg(:))
        !Nnorm = deltat * sum(time_window(:) * seism_win_h(:) * seism_win_h(:))
        !Nnorm = deltat * sum(time_window(:) * seism_veloc(:) * seism_veloc(:))
        if(abs(Nnorm) > EPS) then
            !ft_bar(:) = -seism_veloc(:) * time_window(:) / Nnorm
            !ft_bar(:) = seism_win_h(:) * time_window(:) / Nnorm
            ft_bar(:) = seism_win_fg(:) / Nnorm
            print *,'Norm =', Nnorm
            if(icomp==1) then
                Nnorm_x(irec)=Nnorm
            endif
            if(icomp==3) then
                Nnorm_z(irec)=Nnorm
            endif
        else
            print *, 'norm < EPS for file '
            print *,'Norm =', Nnorm
            ft_bar(:) = 0.d0
        endif

        do itime =1,NSTEP
           if(icomp == adj_comp) then
              write(11,*) (itime-1)*deltat - t0, ft_bar(itime)
           else
              write(11,*) (itime-1)*deltat - t0, 0.d0
           endif
        enddo
    enddo      ! end the loop of comp
    close(11)
    close(12)
    close(13)
    close(14)
    close(16)
    close(17)
    close(18)
    close(19)
    deallocate(taper)
    write(20,*) irec,SQRT(Nnorm_z(irec)/Nnorm_x(irec))
enddo      ! end the loop of rec
print *,'*************************'
print *,'The input files (S****.AA.BXX/BXY/BXZ.adj) needed to run the adjoint simulation are in SEM'
print *,'*************************'
end program adj_seismogram_AH

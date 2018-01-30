subroutine window_taper(npts,taper_percentage,taper_type,tas)
    use constants
    implicit none

    ! input parameters
    integer, intent(in) :: npts
    real, intent(in) :: taper_percentage
    character(len=4) :: taper_type
    real, dimension(*), intent(out) :: tas

    integer :: taper_len
    integer :: i

    ! initialization
    tas(1:npts)=1.0

    if (taper_percentage <= 0.0 .or. taper_percentage >= 1.0) then
        taper_len = int(npts*taper_percentage / 2.0)
    else
        taper_len = int(npts*taper_percentage / 2.0 + 0.5)
    endif

    print *,'taper_len=',taper_len    
    do i=1, taper_len
    if (trim(taper_type) == 'boxc') then
        tas(i)=1.0
    elseif (trim(taper_type) == 'hann') then
        tas(i)=0.5 - 0.5 * cos(2.0 * PI * (i-1) / (2 * taper_len - 1))
    elseif (trim(taper_type) == 'hamm') then
        tas(i)=0.54 - 0.46 * cos(2.0 * PI * (i-1) / (2 * taper_len - 1))
    elseif (trim(taper_type) == 'cose') then
        tas(i)=cos(PI * (i-1) / (2 * taper_len - 1) - PI / 2.0) ** ipwr_t
    elseif (trim(taper_type) == 'cosp') then
        tas(i)=1.0 - cos(PI * (i-1) / (2 * taper_len - 1)) ** ipwr_t
    else
        print*,'taper_type must be among "boxc"/"hann"/"hamm"/"cose"/"cosp"!'
    endif
    tas(npts-i+1)=tas(i)
    enddo

end subroutine window_taper

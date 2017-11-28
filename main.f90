
        program main
        use var !mod file containing all variables
        implicit none
            real (kind=8) :: tms
            integer :: i, j
            
            !declare in file
            in_file = "nmr.in"
            
            !read in nmr file ***ALL VARIABLES READ IN CORRECTLY
            call nmr()

            !read in data file containg x and y points to be analyzed
            !***ALL DATA READ IN AND STORED PROPERLY***
            num_pts = 0
            call read_datfile()

            !find TMS peak, return shift variable
            !assumes data is read in in descending xpt values
            ! ***SHIFT LOCATION FOUND PROPERLY***
            shift = tms()

            !shifts x points over so tms is at 0 on x-axis
            !shifts y points over so baseline is new 0
            do i=1, num_pts
                !xpt(i) = xpt(i) + shift
                !ypt(i) = ypt(i) - inf%base
            end do

            !filters data with user designated filter
            !***ALL FILTERS AND INDEXES WORK PROPERLY***
            select case(inf%filter)
            case (0)
                continue
            case (1)
                call boxcar()
            case (2)
                call SGFilter()
            end select

            !creates coefficents for a spline function between every
            !data point
            !***SPLINE FUNCTION WORKS***!
            call cubic_spline()

        end program

        subroutine find_root()
        use var
        real(kind=8), dimension(100000) :: spline_pts
        real (kind=8) :: next_x, spline_x, y_val, x_val, s, f
        integer :: genpts

            genpts = 80000
            i=1
            j=1
            s = (xpt(1) - xpt(num_pts))/genpts
            next_x = xpt(j+1)
            spline_x = xpt(j)

            do while(i<genpts+1)
                ! if step size has put x at another point, change the
                ! spline function being called
                if(spline_x .lt. next_x) then
                    j = j + 1
                    next_x = xpt(j+1)
                end if
                x_val = spline_x - xpt(j)
                y_val=f(A(j), B(j), C(j), D(j), x_val)
                i = i + 1
                spline_x = spline_x - s
            end do

        end subroutine 

        !subroutine to read in nmr data and place in fortran struct
        subroutine nmr()
        use var

            open (unit=1, file= in_file, status = "old")
            read(1, *) inf%dat_file
            read(1, *) inf%base
            read(1,*) inf%tol
            read(1,*) inf%filter
            read(1,*) inf%size
            read(1,*) inf%passes
            read(1,*) inf%integrate
            read(1,*) inf%out_file

        end subroutine nmr

        !subroutine to read in data file containing x and y points
        subroutine read_datfile()
        use var
            integer :: errval
            real(kind=8) :: x, y

            open(unit=15, file=trim(inf%dat_file), status="old")
            errval=0
            do while(errval /= -1)
                read(15,*,iostat=errval) x, y
                if(errval == 0) then
                    num_pts = num_pts + 1
                    xpt(num_pts) = x
                    ypt(num_pts) = y
                else
                    close(15)
                end if
            end do

        end subroutine

        !subroutine finds the most positive peak with the largest
        ! x-value
        !assumes data is read in in descending xpt values
        real *8 function tms()
        use var
        real(kind=8) :: temp_x, temp_y
        integer :: i, x

            !loops until it finds max peak, assuming x points go from
            !largest to smallest in array
            do i=2, num_pts
                temp_x = xpt(i)
                temp_y = ypt(i)

                if (temp_y .gt. inf%base) then
                    !start searching for peak
                    x = i
                    do while(ypt(x) .lt. ypt(x+1))
                        x = x + 1
                    end do
                    exit
                end if
            end do
            
            tms = 0.0 - xpt(x)
    
        return
        end function

        !returns correct index to run filter across
        subroutine filt_index(curr_pt, index_arr)
        use var
            integer, dimension(*) :: index_arr
            integer :: curr_pt, i, start

            !gets starting index for summation
            start = curr_pt - ((inf%size)-1)/2
    
            i = 1
            do while(i .le. inf%size)
                if(start <= 0) then
                    index_arr(i) = num_pts + start
                    start = start + 1
                else
                    index_arr(i) = start
                    start = start + 1
                end if
                i = i + 1
            end do

        end subroutine

        subroutine boxcar()
        use var
            real (kind=8) :: sum, fsize
            integer, allocatable :: index_arr(:)
            integer :: test_pts, i, j, k

            allocate(index_arr(1:inf%size))
            fsize = inf%size

            do i=1, inf%passes
                do j=1, num_pts
                    call filt_index(j, index_arr)
                    sum = 0
                    do k=1, inf%size
                        sum = sum + ypt(index_arr(k))
                    end do
                    ypt(j) = sum/fsize
                end do
            end do

        end subroutine

        subroutine SGFilter()
        use var
            real(kind=8):: norm_fact, sum
            integer, allocatable :: W(:), index_arr(:)
            integer:: i, j, k
    
            allocate( W(1:inf%size))
            allocate(index_arr(1:inf%size))       
 
            call SGinf(norm_fact, W)

            do i=1, inf%passes
                do j=1, num_pts
                   call filt_index(j, index_arr)
                   sum = 0
                    do k=1, inf%size
                        sum = sum + W(k) * ypt(index_arr(k))
                    end do
                    ypt(j) = sum/norm_fact
                end do
            end do 

        end subroutine

        subroutine cubic_spline()
        use var
        real(kind=8), dimension(10000) :: alpha, iota, mu, zeta, h
        integer :: i, j

            do i=1, num_pts
                A(i) = ypt(i)
            end do

            do i=1, num_pts-1
                h(i) = xpt(i+1) - xpt(i)
            end do

            do i=2, num_pts-1
                alpha(i)=(3.0*(A(i+1)-A(i)))/h(i)-(3.0*(A(i)-A(i-1)))/h(i-1)
            end do

            iota(1) = 1.0
            mu(1) = 0.0
            zeta(1) = 0.0

            do i=2, num_pts-1
                iota(i)=2.0*(xpt(i+1)-xpt(i-1))-h(i-1)*mu(i-1)
                mu(i) = h(i)/iota(i)
                zeta(i) = (alpha(i)-h(i-1)*zeta(i-1))/iota(i)
            end do

            iota(num_pts) = 1.0
            zeta(num_pts) = 0.0
            C(num_pts) = 0.0

            do j=num_pts-1, 1, -1
                C(j) = zeta(j)-mu(j)*C(j+1)
                B(j)=(A(j+1)-A(j))/h(j) - h(j)*(C(j+1)+2.0*C(j))/3.0
                D(j) = (C(j+1) - C(j))/(3*h(j))
            end do

        end subroutine

        !assignes correct normalization factor and  populates W array
        subroutine SGinf(norm_fact, W)
        use var
            real(kind = 8) :: norm_fact
            integer, dimension(*) :: W
           
            if(inf%size == 5) then
                norm_fact = 35.000
                W(1:inf%size) = (/-3, 12, 17, 12, -3 /)
            else if(inf%size == 11) then
                norm_fact = 429.000
                W (1:inf%size)= (/-36, 9, 44, 69, 84, 89, 84, 69, 44, 9,-36/)
            else if(inf%size == 17) then
                norm_fact = 323.000
                W(1:inf%size)=(/-21,-6,7,18,27,34,39,42,43,42,39,34,27,18,7,-6,-21/)
            end if

        end subroutine SGinf

        !returns y-value from spline function
        real *8 function f(Ai, Bi, Ci, Di, x_val)
        real (kind=8) :: Ai, Bi, Ci, Di, x_val

        f= Ai + Bi*(x_val) + Ci*(x_val**2) + Di*(x_val**3)

        return
        end function

        


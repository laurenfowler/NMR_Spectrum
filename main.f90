
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
            do i=1, num_pts
                xpt(i) = xpt(i) + shift
                ypt(i) = ypt(i) - inf%base
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

            call find_root()

            !call Fourier_Transform()

            !Trapezoid method works
            call trapezoid()

            !Crazy numbers..
            !call romberg()

        end program

        !bisection algorithm
        real *8 function bisection(aa,bb,s_loc)
        use var
        real(kind=8) :: aa, bb, FA, FP, a_end, b_end, f, p
        integer :: i, x, s_loc
  
        x = s_loc
        a_end = aa
        b_end = bb

        i=1
        FA = f(A(x),B(x),C(x),D(x),a_end-xpt(x))        
        do while(i<5000)
            p = a_end + (b_end-a_end)/2.0
            FP=f(A(x),B(x),C(x),D(x),p-xpt(x))
            if((FP==0.0) .or. ((b_end-a_end)/2.0 .le. inf%tol)) then
                bisection = p
                exit
            else
                i=i+1
                if((FA*FP) .gt. 0) then
                    a_end = p
                    FA = FP
                else
                    b_end = p
                end if
            end if
        end do
        return 
        end function

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

        subroutine Fourier_Transform()
        use var
        complex *16, dimension(10000) :: y
        complex *16, dimension(10000,10000) :: cc, z, g
        complex *16 :: cnum, cmpx
        integer, dimension(10000) :: IPIV
        integer :: j, k, LDA, INFO, N
        real(kind=8) :: tmp


            N = num_pts

            LDB = num_pts
            cnum = cmplx(cos((2.0*pi)/num_pts), -sin((2.0*pi)/num_pts))
            y = ypt
            cmpx = cnum/num_pts

            !solve for z matrix
            do j=1, num_pts
                do k=1, num_pts
                    tmp = j*k
                    z(j,k) = ((exp(cmpx))**tmp)/(num_pts**(1.0/2.0))
                end do
            end do

            !assign z to cc 
            !cc = z

            !use linear solver to calcuate c matrix
            !call ZGESV(N, 1, cc, 1, IPIV, y, LDB, INFO)
            !print *, INFO  

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

        subroutine find_root()
        use var
        real (kind=8) :: next_x, spline_x, y_val, x_val, s, gp, f
        real (kind=8) :: bisection
        integer :: genpts, x, s_loc, i, j, next_loc

            !generate points to begin searching for roots
            genpts = 80000
            gp = genpts
            i=1
            j=1
            s = (xpt(1) - xpt(num_pts))/gp
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
                
                !assign interpolated values to array
                s_xpts(i) = spline_x
                s_ypts(i) = y_val
                i = i + 1
                spline_x = spline_x - s
            end do

            s_loc = 1
            next_loc = 2
            x = 1

            !begin searching for roots
            do i=1, genpts-1
                if(s_xpts(i) .lt. xpt(next_loc)) then
                    s_loc = s_loc + 1
                    next_loc = next_loc + 1
                end if

                if((s_ypts(i).gt.0).and.(s_ypts(i+1).lt.0))then
                    roots(x)=bisection(s_xpts(i),s_xpts(i+1),s_loc)
                    s_func(x) = s_loc
                    x = x + 1
                else if((s_ypts(i).lt.0).and.(s_ypts(i+1).gt.0))then
                    roots(x)=bisection(s_xpts(i),s_xpts(i+1),s_loc)
                    s_func(x) = s_loc
                    x = x + 1
                else
                    continue
                end if
            end do
    
            num_rts = x-1

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

        subroutine trapezoid()
        use var
        real(kind=8) :: h, next_pt, from_pt, curr_area, fa, fb, f
        real(kind=8) :: spline_x
        integer :: c_spl, indx

        !index for area matrix
        indx = 1
        
        do i=num_rts, 4, -2
            !stop condition for integration
            stop_x = roots(i-1)

            !inital do loop conditions
            !starts at first root location
            from_pt = roots(i)
            
            !first spline function
            c_spl = s_func(i)             

            !ends at next xpt from root location
            to_pt = xpt(c_spl-1)

            !step size
            h = abs(from_pt-to_pt)

            !spline x
            spline_x = xpt(c_spl)

            curr_area = 0
            
            do while(to_pt .lt. stop_x)
                fa = f(A(c_spl),B(c_spl),C(c_spl),D(c_spl), spline_x-from_pt) 
                fb = f(A(c_spl),B(c_spl),C(c_spl),D(c_spl), spline_x-to_pt)

                curr_area = curr_area + ((fa+fb) * (h/2.0))

                c_spl = c_spl - 1
                from_pt = to_pt
                to_pt = xpt(c_spl-1)
                h = abs(from_pt-to_pt)
                spline_x = xpt(c_spl)
            end do
            area(indx) = curr_area
            indx = indx + 1
        end do

        do i=1, indx-1
            print *, area(i)
        end do

        end subroutine

        subroutine romberg()
        use var
        real (kind=8) :: f, from_pt, to_pt, fa, fb
        real(kind=8) :: spline_x, summat, sum_val
        real(kind=8), allocatable :: R(:,:)
        integer :: c_spl, indx, n, e
    
        !index for area matrix
        indx = 1
        n=5
        allocate(R(1:5,1:5))

        do i=num_rts, 21, -2
            !stop condition for integration
            stop_x = roots(i-1)

            !inital do loop conditions
            !starts at first root location
            from_pt = roots(i)
            
            !first spline function
            c_spl = s_func(i)             

            !ends at next xpt from root location
            to_pt = xpt(c_spl-1)

            !spline x
            spline_x = xpt(c_spl)
            
            !step size
            h = abs(from_pt-to_pt)

            do while(to_pt .lt. stop_x)
                print *, "points", from_pt, to_pt
                print *, " "

                fa = f(A(c_spl), B(c_spl), C(c_spl),D(c_spl),spline_x-from_pt)
                fb = f(A(c_spl), B(c_spl), C(c_spl),D(c_spl),spline_x-to_pt)

                R(1,1) = (h/2.0) * (fa + fb)
                curr_area = R(1,1)
                print *, R(1,1)

                do j=2, n
                    sum_val = summat(j,h,from_pt, c_spl)
                    R(j,1)=(1.0/2.0)*(R(j-1,1)+(h*sum_val))
                    do e=2, j
                        R(j,e)=R(j,e-1)+((R(j,e-1)-R(j-1,e-1))/((4**(e-1))-1))
                    end do
                print*, R(j, 1:j)
                h = h/2.0
                curr_area = R(j,e)
                

                !print *, "R(j-1,j-1)", R(j-1,j-1)
                !print *, "R(j,j)", R(j,j)
                !print *, abs(R(j-1,j-1)-R(j,j))

                if(abs(R(j-1,j-1)-R(j,j)) < inf%tol)then
                    print *, "exit condition"
                    exit
                end if

                end do

                c_spl = c_spl - 1
                from_pt = to_pt
                to_pt = xpt(c_spl-1)
                h = abs(from_pt-to_pt)
                spline_x = xpt(c_spl)

                print *, " "
            end do
            area(indx) = curr_area
            indx = indx + 1
        end do        

        end subroutine

        !returns y-value from spline function
        real *8 function f(Ai, Bi, Ci, Di, x_val)
        real (kind=8) :: Ai, Bi, Ci, Di, x_val

        f= Ai + Bi*(x_val) + Ci*(x_val**2) + Di*(x_val**3)

        return
        end function

        !for rhomberg
        real *8 function summat(j, h, aa, c_spl)
        use var
        real(kind=8) :: sum, h, aa, fa
        integer :: upper_index, c_spl, j

        upper_index = 2 * (j-2)
        sum = 0

        do k=1, upper_index
            ak = aa + ((k - 0.5)*h)
            fa = f(A(c_spl),B(c_spl),C(c_spl),D(c_spl),xpt(c_spl)-ak)
            sum = sum + fa
        end do

        summat = sum

        return
        end function

        



        module var
        implicit none

            !used for reading nmr data
            Type nmr_info 
                character (len=20) :: dat_file, out_file
                double precision :: base
                real (kind=8) :: tol
                integer :: filter, size, passes, integrate
            end type

            !x and y point arrays
            real (kind=8), dimension(10000) :: xpt, ypt
            
            !declare type 
            Type (nmr_info) :: inf

            !in_file variable
            character(len=6) :: in_file
            
            !num points variable
            integer :: num_pts 

            !variable used to shift all x-values so tms peak is on 0
            real(kind=8) :: shift  

            !variables that hold spline coefficents
            real(kind=8), dimension(10000) :: A, B, C, D

            !variable that holds all root locations
            real(kind=8), dimension(200) :: rt_loc

        end module

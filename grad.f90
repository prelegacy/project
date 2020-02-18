module grad
    implicit none
    
contains

    subroutine grad_a(rlength,tlength)
        integer, intent(in) :: rlength, tlength
        integer :: N, J

        !Import the lengths of vector r and t to make sure that they correspond to the accretion and is assigned to N and J
        ! print*,' rlength is ', rlength
        ! print*, ' tlength is', tlength
        N = rlength
        J = tlength


    end subroutine grad_a
    
end module grad
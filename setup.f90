module setup 
    implicit none
    contains
    subroutine setupinitial(k, rho, c, fal, ffe, al_ab, fe_ab, tau_al, tau_fe, E_al, E_fe,r,t,dt,dr,al,fe,P,init,bdry,Hin,M)
        real, allocatable, dimension(:), intent(out) :: k, rho, fal, ffe,c ,r,t,dt,dr,P
        real,allocatable,dimension(:,:),intent(out)::Hin, M
        ! integer, intent(out) :: n, ntotal
        real, intent(out) :: al_ab, fe_ab,tau_al, tau_fe, E_al, E_fe,al,fe,init,bdry
        integer :: i,n, endtime, starttime
        !real :: xmin, xmax, space, gamma
        real, parameter::  pi=4*atan(1.), stab = 0.01

        !allocate arrays for our values
       ! print*, 'allocating arrays of size 3, (1) normal condrite, (2) regolith, (3) post melt condrite'
        allocate(k(3), rho(3), fal(3), ffe(3),c(4),P(4),Hin(4,5),M(4,5))
        
            !print*,'allocating thermal conductivity, k' 
            k = (/3.2e7,3.2e6,12.6e7/)

           ! print*,'allocating bulk density, rho'
            rho = (/3440,3400,3290/)

            !print*,'allocating fal'
            fal = (/2.53e23,2.66e23,2.65e23/)

            !print*,'allocating ffe'
            ffe = (/2.96e24,2.36e24,2.12e24/)

            !print*,'allocating specific heat capacities, c'
            !note specific heat capacities array is set up as (1) silicates, (2) metals, (3) sulfide, (4) conjoined grains
            c = (/892,598,699,616/)

            !print*,'allocating weight fraction of chondrite phases, P'
            !silicates, metal, sulfide, conjoined grains
            P =(/0.76,0.05,0.03,0.16/)


            n = (100e3)/500
            allocate(r(n),dr(n))
                r= (/(i,i=0,INT(100e3),500)/)
            do i = 1, n
                if ( i < n ) then
                    dr(i) = r(i+1) -r(i)
                else if (i == n) then
                    dr(i)= r(i) - r(i-1)
                end if                
                
            end do


            endtime = 60e6!500e6
            starttime = 2.85e6
            !original values are endtime = 240e6, 12000
            n = INT((endtime-starttime)/((dr(1)**2)*stab))
            !n = INT((endtime-2.85e6)/(12000))
            print*,'n =',n
            ! n = INT((60e6-2.98e6)/12000)
            ! print*, 'n =',n
            allocate(t(n),dt(n))
            print*,Size(t),Size(dt)

            ! t = (/(i,i=INT(2.98e6),INT(60e6),n)/) ! code version
            ! paper version

            ! t = (/(i,i=INT(2.85e6),INT(endtime),12000)/)
            t = (/(i,i=INT(starttime),INT(endtime),(INT(endtime)-INT(starttime))/(n))/)
            print*,Size(t),Size(dt)
            do i = 1, SIZE(dt)
                if ( i < SIZE(dt) ) then
                    dt(i) = t(i+1) -t(i)
                    
                else if (i == SIZE(dt)) then
                    dt(i)= t(i) - t(i-1)
                end if  
                
            end do

            !Latent heat fusion for phases (Jkg^-1)
            !The order is Silicates, metal, sulfide, conjoined grains
            Hin= reshape((/20000,252000,358000,3797,20000,0,0,4250,40000,0,0,4647,120000,0,0,9308,200000,0,0,214817/),shape(Hin)) !(/20000,20000,40000,120000,200000,252000,0,0,0,0,358000,0,0,0,0,3797,4250,4647,9308,214817/)    
            
            !Temperature at each melting step, in K
            !The order is Silicates, metal, sulfide, conjoined grains
            M = reshape((/1353,1809,1463,1236,1393,0,0,1400,1483,0,0,1500,1753,0,0,1600,1913,0,0,1702/),shape(M))

        !abundance of Al in kg
        al_ab = 2.53e23
        !abundance of Fe in kg
        fe_ab = 2.96e24
        !Mean life of Al in yr
        tau_al = 1.07e6
        !Mean life of Fe in yr
        tau_fe = 3.49e6
        !Decay energy per atom in J
        E_al = 6.4154e-13
        !Decay energy per atom in J
        E_fe = 4.87e-13
        !Al ratio
        al = 5e-5
        !Fe ratio
        fe = 6e-7
        !initial temp of asteroid K
        init = 180
        !External temp of env K
        bdry = 180


    end subroutine setupinitial
        
        
    end module setup
    
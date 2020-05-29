module setup 
    use functions
    implicit none

    contains
    subroutine setupinitial(k, rho, c, fal, ffe, al_ab, fe_ab, tau_al, tau_fe, E_al, E_fe,r,t,dt,dr,al,fe,P,init,bdry,Hin,M,melting)
        real, allocatable, dimension(:), intent(out) :: k, rho, fal, ffe,c ,r,t,dt,dr,P
        real,allocatable,dimension(:,:),intent(out)::Hin, M
        ! integer, intent(out) :: n, ntotal
        real, intent(out) :: al_ab, fe_ab,tau_al, tau_fe, E_al, E_fe,al,fe,init,bdry
        integer :: i,n, endtime, starttime
        integer,intent(out):: melting
        real, parameter::  pi=4*atan(1.), stab = 0.01

        !allocate arrays for our values
       ! print*, 'allocating arrays of size 3, (1) normal condrite, (2) regolith, (3) post melt condrite'
        allocate(k(3), rho(3), fal(3), ffe(3),c(4),P(4),Hin(4,5),M(4,5))
        
            
        print*,'What type of melting?'
        print*,'(1) Regular'
        print*,'(2) Exclue metal sulfide melting'
        print*,'(3) Exclude Conjoined melting'
        print*,'(4) Use alternate C'
        read*, melting


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
        

        endtime = 120e6!500e6
        select case(melting)
        case(1) !regular case
            starttime = 2.98e6
        case(2) !tacc for no metal sulfide 
            starttime = 2.353e6
        case(3) !No conjoined melting
            M = reshape((/1353,1809,1463,2000,1393,0,0,2000,1483,0,0,2001,1753,0,0,2002,1913,0,0,2003/),shape(M))
        case(4) !Alt C values
            c = (/892.93,598.21,698.89,790.64/)
        end select
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
        if (melting == 2) then
            !do nothing
        else
            starttime = 2.98e6
        endif
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
        Hin= reshape((/20000,252000,358000,37978,20000,0,0,4250,40000,0,0,4647,120000,0,0,9308,200000,0,0,214817/),shape(Hin)) !(/20000,20000,40000,120000,200000,252000,0,0,0,0,358000,0,0,0,0,3797,4250,4647,9308,214817/)    
        
        !Temperature at each melting step, in K
        !The order is Silicates, metal, sulfide, conjoined grains
        if (melting ==3) then
            !do nothing
        else 
            M = reshape((/1353,1809,1463,1236,1393,0,0,1400,1483,0,0,1500,1753,0,0,1600,1913,0,0,1702/),shape(M))
        endif
        ! M = reshape((/1353,1809,1463,2000,1393,0,0,2000,1483,0,0,2001,1753,0,0,2002,1913,0,0,2003/),shape(M))
        
        if (melting ==4) then
            !do nothing
        else
            c = (/892,598,699,616/)
        endif
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
    
    subroutine gradinitial(k,reg,rho,c,P,init,bdry,Hstart,Hstart_imp,M,Z,final_rad,rvals, rstep_tot,t_acc,t_dur,tfin,tvals &
        ,tstep_dur,tstep_fin,N,J,tstep_tot,temps_time,rad,tac,delt,delx,delxx,tcounter,rcounter,deltt,tac_final,melting)
        real,allocatable,dimension(:),intent(out):: k,rho, c,P,Hstart_imp,rvals,tvals,tcounter, rcounter,delxx,deltt
        ! real(kind=8), allocatable,dimension(:),intent(out) ::
        real,allocatable,dimension(:,:),intent(out):: Hstart,M,N,J,temps_time,rad,tac,delt,delx,tac_final
        integer, intent(out):: reg,Z,rstep_tot,tstep_dur,tstep_fin,tstep_tot,melting
        real, intent(out):: init, bdry,final_rad,t_acc,t_dur,tfin
        integer ::nz, zval,i,nN,iu,nJ,ncounter
        character(len=25) :: filename
        real :: number, stab1,stab2
        ! real,intent(out):: rho

        print*,'gradinitial online'

        print*,'What type of melting?'
        print*,''
        print*,'(1) Regular'
        print*,''
        print*,'(2) Exclue metal sulfide melting'
        print*,''
        print*,'(3) Exclude Conjoined melting'
        print*,''
        print*,'(4) Use alternate C'
        print*,''
        read*, melting

       
        ! case(3) !No conjoined melting
        !     M = reshape((/1353,1809,1463,2000,1393,0,0,2000,1483,0,0,2001,1753,0,0,2002,1913,0,0,2003/),shape(M))
        ! case(4) !Alt C values
        !     c = (/892.93,598.21,698.89,790.64/)
        ! end select

        !Thermal conductivity in J/(yr*m*K)
        allocate(k(3))
        k = (/3.2e7,3.2e6,12.6e7/)
        
        !Regolith thickness (in number of radial space steps)
        reg = 2

        !Average bulk density of chondrite, previous setup allocated 3 values, but only one was used anyway
        allocate(rho(3))
        !Values of density for H, L, and LL bodies
        rho = (/3440,3400,3290/)

        !Specific heat capacities [silicates, metal, sulfide, conjoined grains]
        if (melting == 4) then
            !Alt C values
            c = (/892.93,598.21,698.89,790.64/)
        else 
            !Regular case
            c = (/892,598,699,616/)
        endif

        !Weight fraction of phases in chondrite [silicates, metal, sulfide, conjoined grains]
        P = (/0.76,0.05,0.03,0.16/)

        !Initial temperature throughout asteroid and initial temperature of newly accreting material
        init = 180

        !External temperature throughout time
        bdry = 180

        !Latent heat of fusion for phases in J kg^-1, in steps
        !Each step is Htot value multiplied by the percentage of the phase that melts
        allocate(Hstart(4,5))
        Hstart=reshape((/20000,252000,358000,136282,20000,0,0,23313,40000,0,0,27815,120000,0,0,65331,300000,0,0,52260/),&
        shape(Hstart))

        !Latent heat of fusion, differnt format of Hstart to setup Hin values
        allocate(Hstart_imp(20))
        Hstart_imp=(/20000,20000,40000,120000,200000,252000,358000,136282,23313,27815,65331,52260/)

        !Temperature (K) at each melting step
        allocate(M(4,5))
        if (melting == 3) then
            !No conjoined melting case
            M = reshape((/1353,1809,1463,2000,1393,0,0,2000,1483,0,0,2001,1753,0,0,2002,1913,0,0,2003/),shape(M))
        else 
            !Regular case
            M = reshape((/1353,1809,1463,1236,1393,0,0,1400,1483,0,0,1500,1753,0,0,1600,1913,0,0,1702/),shape(M))
        endif   
        
        !Set up accretion conditions
        !----------------------------------
        !# of accretion steps
        Z = 50
        
        !FInal radius
        final_rad = 100e3
        !----------------------------------

        !Values of radius (m) at each accretion step
        allocate(rvals(INT(z)))
        rvals = (/(i,i=INT(final_rad/z),INT(final_rad),INT(final_rad/z))/)
       
        !Final number of space steps once radius is at full size
        rstep_tot = INT(final_rad/(1e3))+1
       
        !Inputs
        !Accretion time (Myr after CAIs)
        if(melting == 2) then
            !No metal sulfide
            t_acc = 2.353e6
        else
            !Regular case
            t_acc = 2.55e6 !non-accretion code used t_acc value 2.98e6
        endif

        !Accretion duration (Myr)
        t_dur = 0.35e6

        !Final time (MYR after cais) to compute temp out to
        tfin = 500e6!40e6

        !Values of time (yr) at each accretion step
        allocate(tvals(INT(z)))
        tvals = (/(i,i=INT(t_acc),INT(t_acc+t_dur),INT((t_dur)/Z))/)
        
        !Number of timesteps used between accretion steps
        tstep_dur = 201

        !# of timesteps used after accretion finishes
        tstep_fin = INT(tfin/1e4+1)
        
        !Space steps at each accretion step
        allocate(N(1,Z))
     
        !Time steps at each accretion step
        allocate(J(1,Z))

        !Total number of time steps used
        tstep_tot = tstep_dur*z +tstep_fin
        
        !Creates matrix that tracks timestep in column 1 and temperatures for remaining radial positions
        allocate(temps_time(tstep_tot,rstep_tot+1))        

        !Create a matrix to hold the values of rad and tac
        allocate(rad(Z,INT(rvals(Z)/1000+1)))
        ! allocate(tac(Z,tstep_fin)) !Old values
        allocate(tac(Z,tstep_dur))
        allocate(tac_final(Z,tstep_fin))
        
        ! print*,'shape of temps_time',shape(temps_time)
        ! print*,'shape of rad', shape(rad)
        ! print*,'shape of tac', shape(tac)
        ! print*,'shape of tac_final', shape(tac_final)

        !For all accretion steps
        do nz =1,Z
            zval = nz
           
            !Space steps for this accretion step - set up to check for stability
            do nn = 1,INT(rvals(nz)/1000)+1
                rad(nz,nn)= nn*1000 -1000
            enddo
            !if accretion is finished, time steps go out to tfin (Myr)
                do nn = 1, tstep_dur
                    tac(nz,nn) = INT(tvals(nz)+((tvals(nz+1)-tvals(nz))/tstep_dur)*nn)
                enddo
        enddo

        !Computing remaining time steps after accretion has ended
            do nn = 1,tstep_fin
                tac_final(Z,nn) = INT(tvals(50))+INT(((tfin-tvals(50))/tstep_fin)*nn)
            enddo
        !If accretion is continuing, time steps go between accretion step and next

        !Set length of accretion steps in time and space
        nN = 50
        nJ = 50 

        !Find dt and dx values 
        allocate(delx(Z,INT(rvals(Z)/1000+1)))
        allocate(delxx(Z),deltt(Z))
        allocate(delt(Z,tstep_fin))
        
        do i = 1, SIZE(rad(:,1))
            do nj = 1,SIZE(rad(1,:))-1
                !Sets the term to 0 instead of giving a negative difference in the term
                if (rad(i,nj) ==0) then
                    delx(i,nj) = 0 
                else if (rad(i,nj+1) ==0 .and. rad(i,nj) /= 0) then
                    delx(i,nj) = rad(i,nj) - rad(i,nj-1)
                else
                    !Find the difference in between each step
                    delx(i,nj) = rad(i,nj+1) - rad(i,nj)
                endif
            enddo
        enddo   
        
        do i = 1, SIZE(tac(:,1))
            do nj = 1,SIZE(tac(1,:))-1
                if (tac(i,nj+1) ==0) then
                    delt(i,nj) = 0 
                else
                    delt(i,nj) = tac(i,nj+1) - tac(i,nj)
                endif
            enddo
        enddo

        allocate(rcounter(SIZE(rad(:,1))),tcounter(SIZE(tac(:,1))))

        do nz =1,Z
            Ncounter = 0

            !Note that the shape of the accretion time steps is uniform until accretion has ended, steps 1-49 is 202 steps, and step 50 is 4001 steps
            do nj = 1, SIZE(tac(1,:))
                Ncounter = Ncounter + 1
                if (tac(nz,nj) == 0 .and. nj >1) then
                    tcounter(nz)=Ncounter
                    exit
                else if (nz == Z .and. nj == SIZE(tac(1,:))) then
                    tcounter(nz) = SIZE(tac(1,:))
                endif 
            enddo 
     
    
            !Because the total length is set to max length of final accretion step, we have to find what value the current accretion step reaches, by finding where the radius value becomes 0 after r=0
            Ncounter = 0
            do nj = 1, SIZE(rad(1,:))
                Ncounter = Ncounter + 1
                if (rad(nz,nj) == 0 .and. nj >1) then
                    rcounter(nz)=Ncounter
                    exit
                else if (nz == Z .and. nj == SIZE(rad(1,:))) then
                    rcounter(nz) = SIZE(rad(1,:))
                    !     !Keeps counting
                endif  
                    
            enddo
        enddo

        do i = 1, SIZE(delx(:,1))
            do nj = 1, size(delx(1,:))
                if(nj == rcounter(i)) then
                 delxx(i) = (rcounter(i)*delx(i,2))/rcounter(i)
                endif
            enddo
        enddo

        print*,'size of tac_final is', SIZE(tac_final)

    end subroutine gradinitial
        
    end module setup
    
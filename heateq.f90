module heateq
    use functions
    use grad
    implicit none
    
contains

    subroutine heateqn(temp,r,t,dr,dt,init,bdry,Hin,Hsil,Hmet,Hsulf,Hconj,P,c,k,bulkk,q,fal, al, E_al, tau_al, ffe,fe,E_fe,tau_fe&
        , rho, M)
        real, dimension(:), intent(inout) ::r,t ,dr,dt,P,c,k,fal,ffe,rho
        real, allocatable, dimension(:,:),intent(out) :: temp,Hsil,Hmet,Hsulf,Hconj,bulkk
        real, intent(inout):: al, E_al, tau_al, fe,E_fe,tau_fe
        real,allocatable,dimension(:,:),intent(in)::Hin,M
        real,intent(in)::init,bdry
        real :: bulkc, firstterm,secondterm,thirdterm,stab,factor
        real,dimension(2)::th
        real, intent(out)::q
        integer, parameter:: Reg = 2
        integer:: N, J,I,nN,nJ,L,melting
        
      
        

        !The lengths we will be using to step through each radial distance at each timestep
        N = SIZE(r)
        J = SIZE(t)

        allocate(temp(J,N))
        
        !inputing initial conditions 
        !First row of matrix T
        temp(1,:)=init 
        !last column of matrix T
        temp(:,N)=bdry

        !Allocating residual heat of fusion at each space step (N) - assumes initial chondritic material is solid

        !Silicates
        allocate(Hsil(5,N))
        do I = 1,5
            Hsil(i,:) = Hin(1,i)
        enddo

        !Metals
        allocate(Hmet(5,N))
        do I = 1,5
            Hmet(i,:) = Hin(2,i)
        enddo
        !Sulfides
        allocate(Hsulf(5,N))
        do I = 1,5
            Hsulf(i,:) = Hin(3,i)
        enddo
        
        !Conjoined grains
        allocate(Hconj(5,N))
        do I = 1,5
            Hconj(i,:) = Hin(4,i)
            
        enddo

        !Compute bulk C (specific heat capactiy) as weighted average of specific capacities of each phase
        bulkc = dot_product(p,c)

        !Compute matrix of k (thermal conductivity) at each space step (N)
        allocate(bulkk(1,N))
        bulkk(1,N-Reg:N)=k(1)
        bulkk(1,1:N-Reg)=k(2)
        

        !Computing T at each time step
        !Computing heat term with H abundance of Al and Fe
        do nJ = 1,J-1
            ! print*,'timestep ',nJ
            do nN = 2,N-1

                q = heat(fal(1), al, E_al, tau_al, ffe(1),fe,E_fe,tau_fe, t(nJ+1))
                
                factor = dt(1)/(rho(1)*bulkc*dr(1)**2)

                firstterm = (2*(bulkk(1,nN)*factor)/nN)*(temp(nJ,nN+1)-temp(nJ,nN-1))

                secondterm =  ((bulkk(1,nN)*factor)/nN)*(temp(nJ,nN+1)-(2*temp(nJ,nN))+temp(nJ,nN-1))

                thirdterm = temp(nJ,nN)+((dt(1)/bulkc)*q)

                temp(nJ+1,nN) = firstterm+secondterm+thirdterm

                temp(nJ+1,1)=temp(nJ+1,2)
               
            enddo     
                
            do i = 1,N
                do L = 1,5
                    Hsil(L,i)=Hsil(L,i)-bulkc*(temp(nJ+1,i)-M(1,L))
                    Hconj(L,i)=Hconj(L,i)-bulkc*(temp(nJ+1,i)-M(4,L))
                enddo
                Hmet(1,i)=Hmet(1,i)-bulkc*(temp(nJ+1,i)-M(2,1))
                Hsulf(1,i)=Hsulf(1,i)-bulkc*(temp(nJ+1,i)-M(3,1))
            ! Silicates        
                !Silicates reynolds algorithm to incorporate melting
                do L =1,5
                    th = reynolds(temp(nJ+1,i),M(1,L),Hsil(L,i),Hin(1,L),c(1),P(1))
                    temp(nJ+1,i) = th(1)
                    Hsil(L,i)=th(2)
                    
                enddo

            !Metals     
                th = reynolds(temp(nJ+1,i),M(2,1),Hmet(1,i),Hin(2,1),c(2),P(2))
                temp(nJ+1,i) = th(1)
                Hmet(1,i) = th(2) 

            !Sulfides

                th = reynolds(temp(nJ+1,i),M(3,1),Hsulf(1,i),Hin(3,1),c(3),P(3))
                temp(nJ+1,i) = th(1)
                Hsulf(1,i) = th(2) 
            
            !Conjoined grains
            
                do L =1,5
                    th = reynolds(temp(nJ+1,i),M(4,L),Hconj(L,i),Hin(4,L),c(4),P(4))
                    temp(nJ+1,i) = th(1)
                    Hconj(L,i)=th(2)
                enddo  
                
                !Adjusts values of thermal conductivity for decreasing pore space after partial silicate melting
                if (temp(nJ+1,i) > temp(nJ,i)) then
                    if (i > N - reg) then
                    !does nothing
                    else   
                        if (temp(nJ+1,i) > M(1,5)) then
                            bulkk(1,i) = k(3)
                        else !if T is below silicate solidus
                        !does nothing
                        endif
                    endif
                endif   
            enddo    
        
        enddo
    end subroutine heateqn
    
    subroutine heateqn_a(k,Z,rad,reg,tac,deltt,delxx,temp,init,bdry,c,p,Hin,Hstart_imp,init_array,acc_con,rho,tT,temps_time,&
        bulkk)
        real,allocatable,dimension(:,:):: bulkk
        real, intent(inout)::init,bdry
        real,dimension(:,:),intent(inout):: rad,tac
        real,allocatable,dimension(:,:),intent(inout)::temp,tT,temps_time,Hin
        real,dimension(:),intent(inout):: deltt,delxx,c,p,rho
        real,dimension(:),intent(in) :: Hstart_imp
        real, allocatable,dimension(:),intent(inout) :: init_array
        real, dimension(:), intent(in):: k 
        integer, intent(in):: Z,reg
        integer, intent(out) :: acc_con
        integer,allocatable, dimension(:):: rcounter, tcounter
        integer:: nz,nj,Ncounter, i, iu
        character(len=25) :: filename
        
        !Set up k-values to be used in the program, make same length as radius, so the corresponding K-value can be used
        allocate(bulkk(SIZE(rad(:,1)),SIZE(rad(1,:))))
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
            !If this is the first accretion step, set the initial thermal conductivity
            if (nz == 1) then
                bulkk(nz,1:rcounter(nz)) = k(2)
            
            !If this is the last accretion step
            elseif ( nz == Z ) then

                !Set K computed for existing material
                bulkk(nz,1:rcounter(nz-1)) = 5!thk(14,:)
                ! call grad_a(nz,rcounter(nz),tcounter(nz),delxx,deltt,temp,init,bdry)
                
                !Set K for newly accreted material
                bulkk(nz,rcounter(nz-1):rcounter(nz))=k(2)
                !Set K for regolith
                bulkk(nz, rcounter(nz)- Reg: rcounter(nz)) = k(1)

            !if accretion i scontinuing  
            else 
                ! call grad_a(nz,rcounter(nz),tcounter(nz),delxx,deltt,temp,init,bdry)
                bulkk(nz,1:rcounter(nz-1)) = THK(14,:)!Need to find out what thk is , could be K computed for existing material
                !Set normal K for newly accreted material
                bulkk(nz,rcounter(nz-1):rcounter(nz)) = k(2)
            end if   

            !Set up Hin, the full array is a 12x50 matrix but will be reallocated with each step
            if (z == 1) then
                !Allocate Hin for the length of the current accretion step
                allocate(Hin(12,nz))
                !Set initial Hin for all material at 1st accretion step
                do i = 1,12
                    Hin(i,:) = Hstart_imp(i)
                enddo
            !If this is not the first accretion step
            else
                !Reallocate the length of Hin
                deallocate(Hin)
                allocate(Hin(12,nz))
                do i = 1,12
                    !Sets Hin as the last computed accretion step for existing material
                    Hin(i,1:rcounter(nz-1)) = THK(i+1,:)
                    !Sets initial Hin for newly accreted material
                    Hin(i,rcounter(nz-1):rcounter(nz)) = Hstart_imp(i)
                enddo
            endif
            
            !Setup initial T condition values
            if (nz ==1) then 

                allocate(init_array(rcounter(nz)))

                !Set initial T array to init tmep
                init_array(:) = init

            else 

                !Reallocate array
                deallocate(init_array)
                allocate(init_array(rcounter(nz)))

                !Set initial T array to computed T for existing material
                init_array(1:rcounter(nz-1)) = THK(:)

                !Set iniitail T for newly accreted material
                init_array(rcounter(nz-1):rcounter(nz)) = init

            endif

            !Assign a value to accretion condition (acc_con) - 0 if in progressed, 1 if finished

            if (nz == Z) then

                !If accreiton has finished
                acc_con = 1

            else 

                !If accretion is still occuring 
                acc_con = 0

            endif 

            !Runs the model to find the temperature at each space step through time

            !If this is the first accretion step
            if (nz == 1) then

                !Allocate the length of the tT array
                allocate(tT(tcounter(nz),rcounter(nz)))
                
            !If this is not the first accretion step
            else 

                !Reallocate tT
                deallocate(tT)
                allocate(tT(tcounter(nz),rcounter(nz)))
            endif 
            
            !Might need to fix up
            !tT = heateqn_grad_a(count,rlength,tlength, delxx,deltt,temp,init,bdry,Hin,c,p,tac,rho,bulkk,M,Hstart,acc_con,reg,k)
            call grad_a(count,rlength,tlength, delxx,deltt,temp,init,bdry,Hin,c,p,tac,rho,bulkk,M,Hstart,acc_con,reg,k,tT)
    

            !Fils temps_time matrix with the times and temperatures
            if(nz ==1 ) then

                temps_time(1:tcounter(nz),1:(rcounter(nz)+1)) = tT
            endif



        enddo

        

        ! write(filename,"(a)")'koutput.dat'
        ! print "(a)",' writing to '//trim(filename)
        ! open(newunit=iu,file=filename,status='replace',&
        ! action='write')
        ! write(iu,"(a)") '#  k'
        ! do i=1,SIZE(bulkk(1,:))
        !         write(iu,fmt='(50F15.2)') bulkk(:,i) 
        ! enddo
        ! close(iu)

        !THK = call heat eqn grad B
    end subroutine heateqn_a



    ! heateqn_a(Z,rad,rvals,tac,tvals,tfin,tstep_fin,tstep_dur,delt,delx)
    !     integer,intent(in):: Z,tstep_fin,tstep_dur
    !     real, allocatable,dimension(:,:), intent(out):: rad,tac,delt,delx
    !     real,dimension(:),intent(in)::rvals,tvals
    !     real,intent(in):: tfin
    !     real :: stab1, stab2
    !     integer ::nz, zval,i,nN,iu,nJ,j
    !     character(len=25) :: filename
    !     allocate(rad(Z,INT(rvals(Z)/1000+1)))
    !     allocate(tac(Z,tstep_fin))
        
    !     !For all accretion steps
    !     do nz =1,Z
    !         zval = nz
           
    !         !Space steps for this accretion step - set up to check for stability
    !         do nn = 1,INT(rvals(nz)/1000)+1
    !             rad(nz,nn)= nn*1000 -1000
    !         enddo
    !         !if accreiton is finished, time steps go out to tfin (Myr)
    !         if (nz == Z) then
    !             do nn = 1,tstep_fin
    !                 tac(Z,nn) = INT(tvals(z))+INT(((tfin-tvals(z))/tstep_fin)*nn)
    !             enddo
    !         !If accretion is continuing, time steps go between accretion step and next
    !         else 
    !             do nn = 1, tstep_dur
    !                 tac(nz,nn) = tvals(nz)+INT((tvals(nz+1)-tvals(nz))*nn)
    !             enddo
    !         endif

    !     enddo
        
    !     !Set length of accretion steps in time and space
    !     nN = 50
    !     nJ = 50 

    !     !FInd dt and dx values 
    !     allocate(delx(Z,INT(rvals(Z)/1000+1)))
    !     allocate(delt(Z,tstep_fin))
        
    !     do i = 1, SIZE(rad(:,1))
    !         do j = 1,SIZE(rad(1,:))-1
    !             !Sets the term to 0 instead of giving a negative difference in the term
    !             if (rad(i,j+1) ==0) then
    !                 delx(i,j) = 0 
    !             else
    !             delx(i,j) = rad(i,j+1) - rad(i,j)
    !             endif
    !         enddo
    !     enddo
        
    !     do i = 1, SIZE(tac(:,1))
    !         do j = 1,SIZE(tac(1,:))-1
    !             if (tac(i,j+1) ==0) then
    !                 delt(i,j) = 0 
    !             else
    !                 delt(i,j) = tac(i,j+1) - tac(i,j)
    !             endif
    !         enddo
    !     enddo


     

    !     write(filename,"(a)")'toutput.dat'
    !     print "(a)",' writing to '//trim(filename)
    !     open(newunit=iu,file=filename,status='replace',&
    !     action='write')
    !     write(iu,"(a)") '#  t'
    !     do i=1,SIZE(tac(:,1))
    !             write(iu,fmt='(4001F15.2)') tac(i,:) 
    !     enddo
    !     close(iu)

    !     print*,(SUM(delx)/SIZE(delx))*2
    !     print*,(SUM(delt)/SIZE(delt))*2
    !     print*, MAXVAL(delt), MINVAL(delt)
    !     print*,SHAPE(delt)

    !     stab1 = stability(MAXVAL(delt),MAXVAL(delx)) 
    !     print*, stab1
    !     stab2 = stability(MAXVAL(delt),(SUM(delx)/SIZE(delx))*2) 
    !     print*, stab2

    !     !Only continue on if stab1 and stab2 are less than 0.01
    !     if (stab1 < 0.01 .and. stab2 < 0.01) then
            
    !     endif

        ! !check to see if the output is correct
        ! ! write(filename,"(a)")'routput.dat'
        ! print "(a)",' writing to '//trim(filename)
        ! open(newunit=iu,file=filename,status='replace',&
        ! action='write')
        ! write(iu,"(a)") '#  r'
        ! do i=1,SIZE(rad(1,:))
        !         write(iu,fmt='(50F9.2)') rad(:,i) 
        ! enddo
        ! close(iu)
       
        ! print*,'size =', SIZE(tac(:,1))
        ! write(filename,"(a)")'toutput.dat'
        ! print "(a)",' writing to '//trim(filename)
        ! open(newunit=iu,file=filename,status='replace',&
        ! action='write')
        ! write(iu,"(a)") '#  t'
        ! do i=1,SIZE(tac(:,1))
        !         write(iu,fmt='(50F10.0)') tac(:,i) 
        ! enddo
        ! close(iu)




        ! do nz =1,Z
        !     zval = nz

        !     !Space steps for this accretion step - set up to check for stability
        !     allocate(r(INT(rvals(nz)/1000)+1))
        !     r = (/(i,i=0,INT(rvals(nz)), 1000)/)

        !     !if accreiton is finished, time steps go out to tfin (Myr)
        !     if (nz == Z) then
        !         t = 
        !     endif
        !     deallocate(r)

        ! enddo
end module heateq
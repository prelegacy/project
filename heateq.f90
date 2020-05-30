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
        real :: bulkc, firstterm,secondterm,thirdterm,factor
        real,dimension(2)::th
        real, intent(out)::q
        integer, parameter:: Reg = 2
        integer:: N, J,I,nN,nJ,L
        
      
        

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

                secondterm =  ((bulkk(1,nN)*factor))*(temp(nJ,nN+1)-(2*temp(nJ,nN))+temp(nJ,nN-1))

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
        bulkk,THK,m,Hstart,tstep_dur,tac_final,tstep_fin,tstep_tot,fALs,fFes,Altratio,Feratio)
        real,allocatable,dimension(:,:):: bulkk
        real, intent(inout)::init,bdry
        real,dimension(:,:),intent(inout):: rad,tac, m, Hstart,temps_time,tac_final
        real,allocatable,dimension(:,:),intent(inout)::temp,tT,Hin,thk
        real,dimension(:),intent(inout):: deltt,delxx,c,p,rho
        real,dimension(:),intent(in) :: Hstart_imp
        real, allocatable,dimension(:),intent(inout) :: init_array
        real, dimension(:), intent(in):: k 
        real, intent(inout) :: fALs,fFes,Altratio,Feratio
        integer, intent(in):: Z,reg,tstep_tot
        integer, intent(out) :: acc_con
        integer,intent(inout) :: tstep_dur,tstep_fin
        integer,allocatable, dimension(:):: rcounter!, tcounter
        integer:: nz,nj,Ncounter, i, iu,j
        character(len=25) :: filename
        
        !Set up k-values to be used in the program, make same length as radius, so the corresponding K-value can be used
        allocate(bulkk(SIZE(rad(:,1)),SIZE(rad(1,:))))
        allocate(rcounter(SIZE(rad(:,1))))!,tcounter(SIZE(tac(:,1))))

        do nz = 1, Z
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
        do nz =1,Z+1

            if (nz < z+1) then
                Ncounter = 0

                !Note that the shape of the accretion time steps is uniform until accretion has ended, steps 1-49 is 202 steps, and step 50 is 4001 steps
                ! do nj = 1, SIZE(tac(1,:))
                !     Ncounter = Ncounter + 1
                !     if (tac(nz,nj) == 0 .and. nj >1) then
                !         tcounter(nz)=Ncounter
                !         exit
                !     else if (nz == Z .and. nj == SIZE(tac(1,:))) then
                !         tcounter(nz) = SIZE(tac(1,:))
                !     endif 
                ! enddo 
        

                !Because the total length is set to max length of final accretion step, we have to find what value the current accretion step reaches, by finding where the radius value becomes 0 after r=0
                
                !If this is the first accretion step, set the initial thermal conductivity
                if (nz == 1) then
                    bulkk(nz,1:rcounter(nz)) = k(2)
                
                !If this is the last accretion step
                elseif ( nz == Z ) then

                    !Set K computed for existing material
                    bulkk(nz,1:rcounter(nz-1)) = THK(14,:)
            
                    
                    !Set K for newly accreted material
                    bulkk(nz,rcounter(nz-1):rcounter(nz))=k(2)
                    !Set K for regolith
                    bulkk(nz, rcounter(nz)- Reg: rcounter(nz)) = k(1)

                !if accretion i scontinuing  
                else 
                    
                    bulkk(nz,1:rcounter(nz-1)) = THK(14,:)!Need to find out what thk is , could be K computed for existing material
                    !Set normal K for newly accreted material
                    bulkk(nz,rcounter(nz-1):rcounter(nz)) = k(2)
                end if   

                !Set up Hin, the full array is a 12x50 matrix but will be reallocated with each step
                if (nZ == 1) then
                    !Allocate Hin for the length of the current accretion step
                    allocate(Hin(12,rcounter(nz)))
                    !Set initial Hin for all material at 1st accretion step
                    do i = 1,12
                        Hin(i,:) = Hstart_imp(i)
                    enddo
                !If this is not the first accretion step
                else
                    !Reallocate the length of Hin
                    deallocate(Hin)
                    allocate(Hin(12,rcounter(nz)))
                    do i = 1,12
                        !Sets Hin as the last computed accretion step for existing material
                        Hin(i,1:rcounter(nz-1)) = THK(i+1,:)

                        !Sets initial Hin for newly accreted material
                        !print*,'looking at i', i,'looking at rounter',rcounter(nz)
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
                    init_array(1:rcounter(nz-1)) = THK(1,:)

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
                    allocate(tT(SIZE(tac(nZ,:)),rcounter(nz)+1))                
                    
                !If this is not the first accretion step
                else 

                    !Reallocate tT
                    deallocate(tT)
                    allocate(tT(SIZE(tac(nZ,:)),rcounter(nz)+1))

                endif 
                
                deltt = (tac(nz,2)-tac(nz,1))
                !Might need to fix up   
                call grad_a(nZ,rcounter(nz),tstep_dur, delxx,deltt,temp,init,bdry,Hin,c,p,tac,rho,&
                bulkk,M,Hstart,acc_con,reg,k,tT,thk,init_array,fALs,fFes,Altratio,Feratio)
            
                ! print*,'tcounter(nz) is', tcounter(nz)

                print*,'----------------------'
                !Fils temps_time matrix with the times and temperatures
        
                if(nZ == 1) then
                    temps_time(1:tstep_dur*nZ,1:(rcounter(nz)+1))=tT(:,:)
                else 
                    temps_time(tstep_dur*(nZ-1):(tstep_dur*nZ)-1,1:(rcounter(nz)+1))=tT(:,:)
                endif      

            !     write(filename,"(a,i5.5,a)") 'tt_',nz,'.txt'
            ! print "(a)",' writing to '//trim(filename)
            ! open(newunit=iu,file=filename,status='replace',&
            ! action='write')
            ! write(iu,"(a)") '#  t, r'
            !     do i=1,SIZE(tT(:,1))
            !         write(iu,fmt='(10F15.2)') tT(i,:)
            !     enddo
            ! close(iu)
            else 

                !Call output to do just accretion

                
                print*,'Accretion completed---------------------------'
                !Completing the last time step using tac_final rather than tac

                !If this is the last accretion step
                !Set K computed for existing material
                
                bulkk(nz-1,1:rcounter(nz-1)) = THK(14,:)
                
                !Set K for newly accreted material
                bulkk(nz-1,rcounter(nz-2):rcounter(nz-1))=k(2)
                !Set K for regolith
                bulkk(nz-1, rcounter(nz-1)- Reg: rcounter(nz-1)) = k(1)

                !Set up Hin, the full array is a 12x50 matrix but will be reallocated with each step
                !Reallocate the length of Hin
                deallocate(Hin)
                allocate(Hin(12,rcounter(Z)))
                do i = 1,12
                    !Sets Hin as the last computed accretion step for existing material
                    Hin(i,1:rcounter(Z)) = THK(i+1,:)
                enddo
                    
                deallocate(init_array)
                allocate(init_array(rcounter(Z)))
                !Set initial T array to computed T for existing material
                init_array(1:rcounter(Z)) = THK(1,:)

                acc_con = 1

                !Runs the model to find the temperature at each space step through time
                !Reallocate tT
                deallocate(tT)
                allocate(tT(SIZE(tac_final(Z,:)),rcounter(Z)+1))
                
                deltt = (tac_final(Z,2)-tac_final(Z,1))
              
               
                !Might need to fix up   
                call grad_a(nZ-1,rcounter(nz-1),tstep_fin, delxx,deltt,temp,init,bdry,Hin,c,p,tac_final,rho,&
                bulkk,M,Hstart,acc_con,reg,k,tT,thk,init_array,fALs,fFes,Altratio,Feratio)
            
                
                temps_time(tstep_dur*(Z):tstep_tot-1,1:(rcounter(Z)+1))=tT(:,:)
            
            endif
        enddo
        
      


        !Create several outputs to work with


    end subroutine heateqn_a
end module heateq
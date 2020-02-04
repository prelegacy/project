module heateq
    use functions
    implicit none
    
contains
    subroutine heateqn(temp,r,t,dr,dt,init,bdry,Hin,Hsil,Hmet,Hsulf,Hconj,P,c,k,bulkk,q,fal, al, E_al, tau_al, ffe,fe,E_fe,tau_fe&
        , rho, M)
        real, dimension(:), intent(inout) ::r,t ,dr,dt,P,c,k,fal,ffe,rho
        real, allocatable, dimension(:,:),intent(out) :: temp,Hsil,Hmet,Hsulf,Hconj,bulkk
        real, intent(inout):: al, E_al, tau_al, fe,E_fe,tau_fe
        real,allocatable,dimension(:,:),intent(in)::Hin,M
        real,intent(in)::init,bdry
        real :: bulkc, firstterm,secondterm,thirdterm,fourthterm
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
        do nJ = 1,5!J-1
            ! print*,'timestep ',nJ
            do nN = 2,N-1
                ! print*,'Spacestep ',nN
                q = heat(fal(1), al, E_al, tau_al, ffe(1),fe,E_fe,tau_fe, t(nJ+1))
                firstterm = (2*bulkk(1,nN)*dt(nJ)/(rho(1)*bulkc*(dr(nN)**2)))*(temp(nJ,nN+1)-temp(nJ,nN-1))
                secondterm =  (bulkk(1,nN)*dt(nJ)/(rho(1)*bulkc*(dr(nN)**2)))*(temp(nJ,nN+1)-2*temp(nJ,nN)+temp(nJ,nN))
                thirdterm = temp(nJ,nN)+((dt(nJ)/bulkc)*q)
                temp(nJ+1,nN) = firstterm+secondterm+thirdterm
                temp(nJ+1,1)=temp(nJ+1,2)

               
            enddo
            do i = 1,5
                do L = 1,5
                     Hsil(L,i)=Hsil(L,i)-bulkc*(temp(nJ+1,nN)-M(1,L))
                     Hconj(L,i)=Hconj(L,i)-bulkc*(temp(nJ+1,nN)-M(4,L))
                enddo
                Hmet(1,i)=Hmet(1,i)-bulkc*(temp(nJ+1,nN)-M(2,1))
                Hsulf(1,i)=Hsulf(1,i)-bulkc*(temp(nJ+1,nN)-M(3,1))
                
            enddo
            
        enddo

        print*,temp(5,:)
        
       
        !print*,t

        ! do nJ = 1,1!J-1
        !     ! print*,'timestep ',nJ
        !     do nN = 2,N-1
        !         ! print*,'Spacestep ',nN
        !         print*,-temp(nJ,nN-1)
        !         secondterm =  (bulkk(1,nN)*dt(nJ)/(rho(1)*bulkc*(dr(nN)**2)))*(temp(nJ,nN+1)-2*temp(nJ,nN)+temp(nJ,nN))
        !     enddo
        ! enddo
    end subroutine heateqn
end module heateq
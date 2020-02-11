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
end module heateq
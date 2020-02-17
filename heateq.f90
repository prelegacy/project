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
    
    subroutine heateqn_a(Z,rad,rvals,tac,tvals,tfin,tstep_fin,tstep_dur,delt,delx)
        integer,intent(in):: Z,tstep_fin,tstep_dur
        real, allocatable,dimension(:,:), intent(out):: rad,tac,delt,delx
        real,dimension(:),intent(in)::rvals,tvals
        real,intent(in):: tfin
        integer ::nz, zval,i,nN,iu,nJ,j
        character(len=25) :: filename
        allocate(rad(Z,INT(rvals(Z)/1000+1)))
        allocate(tac(Z,tstep_fin))
        
        !For all accretion steps
        do nz =1,Z
            zval = nz
           
            !Space steps for this accretion step - set up to check for stability
            do nn = 1,INT(rvals(nz)/1000)+1
                rad(nz,nn)= nn*1000 -1000
            enddo
            !if accreiton is finished, time steps go out to tfin (Myr)
            if (nz == Z) then
                do nn = 1,tstep_fin
                    tac(Z,nn) = INT(tvals(z))+INT(((tfin-tvals(z))/tstep_fin)*nn)
                enddo
            !If accretion is continuing, time steps go between accretion step and next
            else 
                do nn = 1, tstep_dur
                    tac(nz,nn) = tvals(nz)+INT((tvals(nz+1)-tvals(nz))*nn)
                enddo
            endif

        enddo
        
        !Set length of accretion steps in time and space
        nN = 50
        nJ = 50 

        !FInd dt and dx values 
        allocate(delx(Z,INT(rvals(Z)/1000+1)))
        allocate(delt(Z,tstep_fin))
        
        do i = 1, SIZE(rad(:,1))
            do j = 1,SIZE(rad(1,:))-1
                !Sets the term to 0 instead of giving a negative difference in the term
                if (rad(i,j+1) ==0) then
                    delx(i,j) = 0 
                else
                delx(i,j) = rad(i,j+1) - rad(i,j)
                endif
            enddo
        enddo
        
        do i = 1, SIZE(tac(:,1))
            do j = 1,SIZE(tac(1,:))-1
                print*,'i =',i,'j =', j
                if (tac(i,j+1) ==0) then
                    delt(i,j) = 0 
                else
                    delt(i,j) = tac(i,j+1) - tac(i,j)
                endif
            enddo
        enddo


        print*,(SUM(delx)/SIZE(delx))*2
        print*,(SUM(delt)/SIZE(delt))*2
        print*, MAXVAL(delt), MINVAL(delt)

        !   write(filename,"(a)")'deltoutput.dat'
        ! print "(a)",' writing to '//trim(filename)
        ! open(newunit=iu,file=filename,status='replace',&
        ! action='write')
        ! write(iu,"(a)") '#  t'
        ! do i=1,SIZE(delx(1,:))
        !         write(iu,fmt='(50F9.2)') delt(:,i) 
        ! enddo
        ! close(iu)





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

    end subroutine heateqn_a
end module heateq
module grad
    use functions
    implicit none
    
contains

    subroutine grad_a(count,rlength,tlength, delxx,deltt,temp,init,bdry,Hin,c,p,t)
        integer, intent(in) :: rlength, tlength,count
        real,intent(in)::init,bdry,t
        real,dimension(:),intent(in):: delxx,deltt,c,p
        real,dimension(:,:),intent(in)::Hin
        real, allocatable,dimension(:,:), intent(out):: temp
        real,allocatable,dimension(:,:)::Hsil,Hmet,Hsulf,Hconj
        real :: bulkC, fAL,fFe,Altratio,Feratio,EAl,EFe,LifeAl,LifeFe,q
        integer :: N, J, dr, dt, i, iu,ni,nj
        character(len=25) :: filename



        !Might move heat source setup later to setup.f90
        
        !Abundance of Al (kg^-1)
        fAL = 2.53e23
        !Abundance of Fe (kg^-1)
        fFe = 2.96e24
        !26Al/27Al initial ratio
        Altratio = 5e-5
        !60Fe/56Fe initial ratio
        Feratio = 6e-7
        !26Al decay energy per atom (J)
        EAl = 6.4154e-13
        !60Fe decay energy per atom (J)
        EFe = 4.87e-13
        !26Al mean life (years)
        LifeAl = 1.07e6
        !60Fe mean life (years)
        LifeFe = 3.49e6

        !Import the lengths of vector r and t to make sure that they correspond to the accretion and is assigned to N and J
        ! print*,' rlength is ', rlength
        ! print*, ' tlength is', tlength
        N = rlength
        J = tlength

        !Size of space step and time step
        dr = INT(delxx(count))
        dt = INT(deltt(count))
        
        !Set up matrix of J-N (time-space) 
        !Filled with temperature values for each space step and time s
        allocate(temp(J,N))

        !Fills first row of matrix T with initial condition
        temp(1,:) = init
        
        !Fill last column of matrix Temp with bundary temperature (Dirichlet boundary condition)
        temp(:,N) = bdry !input into sub
        
        !Set up arrays for residual heat of fusion at each space step (N) - assumes initial chondritic material is solid

        allocate(Hsil(5,N))
        do i = 1,5
        Hsil(i,:) = Hin(i,:)
        enddo

        allocate(Hmet(1,N))
        Hmet(1,:) = Hin(6,:)

        allocate(Hsulf(1,N))
        Hsulf(1,:)=Hin(7,:)

        allocate(Hconj(5,N))
        do i = 1,5
            Hconj(i,:)=Hin(i+7,:)
        enddo
    
        !Compute BUlk C (specific heat capacity) as a weighted average of the specific heat capacities of each phase
        bulkC = dot_product(p,c)
        
        !import t
        !For each element of q (i.e each time step), calculate the value of the heat source term
        q = heat(fAL, Altratio, EAl, LifeAl, fFe,Feratio,EFe,LifeFe, t)

        !for each time step
        do nJ = 1,J-1
            print*,'timestep =',nJ
            !Compute T at next time step
            do nN = 2,N-1
                temp(nJ+1,nN) = !insert term

                !Neumann boundary conditions
                temp(nJ+1,1) = temp(nJ+1,2)
            enddo
        enddo

         
       

    end subroutine grad_a
    
end module grad
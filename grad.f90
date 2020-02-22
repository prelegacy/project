module grad
    use functions
    implicit none
    
contains

    subroutine grad_a(count,rlength,tlength, delxx,deltt,temp,init,bdry,Hin,c,p,t)
        integer, intent(in) :: rlength, tlength,count
        real,intent(in)::init,bdry,t
        real,dimension(:),intent(in):: delxx,deltt,c,p
        real,dimension(:,:),intent(in)::Hin
        real,allocatable,dimension(:)::th
        real, allocatable,dimension(:,:), intent(out):: temp
        real,allocatable,dimension(:,:)::Hsil,Hmet,Hsulf,Hconj
        real :: bulkC, fAL,fFe,Altratio,Feratio,EAl,EFe,LifeAl,LifeFe,q
        integer :: N, J, dr, dt, i, iu,ni,nj,nn,nw
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
        

        !for each time step
        do nJ = 1,J-1
            print*,'timestep =',nJ
            !Compute T at next time step
            do nN = 2,N-1
                !Note that our bulkk k-value is a 2D array which is updated at each accretion step
                q = heat(fAL, Altratio, EAl, LifeAl, fFe,Feratio,EFe,LifeFe, t(nJ+1))
                temp(nJ+1,nN) = ((2*bulkk(nz,nN)*dt)/(rho(1)*bulkC*nN*(dr**2)))*(temp(nJ,nN+1)-temp(nJ,nN-1))+((bulkk(nz,nN)*dt)/ &
                (rho(1)*bulkC*(dr**2)))*(temp(nJ,nN+1)-2*temp(nJ,nN)+temp(nJ,nN-1))+temp(nJ,nN)+(dt/bulkC)*q
                !Neumann boundary conditions
                temp(nJ+1,1) = temp(nJ+1,2)
            enddo

            !Computes residual heat of fusion values 
            do nN = 1,N 
                do nW = 1,5
                    Hsil(nW,nN)=Hsil(nW,nN)-bulkC*(temp(nJ+1,nN)-M(1,nW))
                enddo
                Hmet(1,nN)=Hmet(1,nN)-bulkC*(temp(nJ+1,nN)-M(2,1))
                Hsulf(1,nN)=Hsulf(1,nN)-bulkC*(temp(nJ+1,nN)-M(3,1))
                do nW= 1,5
                    Hconj(nW,nN) = Hconj(nW,nN)-C*(temp(nJ+1,nN)-M(4,nW))
                enddo

            

            !Following seciton uses an algorithm modified from Reynolds et al (1966)
            !Incorporates melting, algorithm is a subfunction called Renolds which is presented below
            !The output is a 1x2 matrix TH = [new temp, new residual heat of fusion]
            !Values from TH are then assigned to temp(nj+1,nN) and Hphase(w,n)

            !Silicates
            allocate(th(2))

            do nW = 1,5
                th = reynolds(temp(nJ+1,nN),M(1,nw),Hsil(nW,nN),Hstart(1,nW),c(1),P(1))
                temp(nJ+1,nN) = th(1)
                Hmet(1,nN) = th(2)
            enddo

            !metals-only
            th = reynolds(temp(nJ+1,nN),M(2,1),Hmet(1,nN),Hstart(2,1),c(2),P(2))
            temp(nj+1,nN) = th(1)
            Hmet(1,nN) = th(2)

            enddo   
        enddo

         
       

    end subroutine grad_a
    
end module grad
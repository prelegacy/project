module grad
    use functions
    implicit none
    
contains

    subroutine grad_a(count,rlength,tlength, delxx,deltt,temp,init,bdry,Hin,c,p,tac,rho,bulkk,M,Hstart,acc_con,reg,k,tT,thk,&
        init_array,fALs,fFes,Altratio,Feratio)
        real,allocatable,dimension(:,:),intent(inout) :: tT
        integer, intent(in) :: rlength, tlength,count,reg
        real,intent(in)::init,bdry
        real,dimension(:),intent(in):: delxx,deltt,c,p,rho,k,init_array
        real,dimension(:,:),intent(in)::Hin,tac,M,Hstart
        real,dimension(:,:),intent(inout) :: bulkk
        real,allocatable,dimension(:)::th
        real, allocatable,dimension(:,:), intent(out):: temp,thk
        real,allocatable,dimension(:,:)::Hsil,Hmet,Hsulf,Hconj
        real :: bulkC,EAl,EFe,LifeAl,LifeFe,q,mass
        real,parameter:: core_density=7870,mantle_density=3976
        real,intent(inout)::fALs,fFes,Altratio,Feratio
        integer :: N, J, dr, dt, i, iu,ni,nj,nn,nw,differentiation
        
        integer, intent(in) :: acc_con
        character(len=25) :: filename

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
        temp(1,:) = init_array(:)
        
        !Fill last column of matrix Temp with bundary temperature (Dirichlet boundary condition)
        temp(:,N) = bdry !input into sub
        
        !Set up arrays for residual heat of fusion at each space step (N) - assumes initial chondritic material is solid

        !Might move heat source setup later to setup.f90
        differentiation = 0
        if (N > 200) then
            if (MAXVAL(temp(:,200))>1200)  then
                differentiation=1
            endif
        endif
        
        print*, 'at accretion step ', count,differentiation
        !Abundance of Al (kg^-1) assuming 1.13%
        !Over estimated the abundance?
        mass = abs(1.333*3.14*(abs(rlength*1e3))**3)
        !if (rlength<=100) then

        fALs = rho(1)* mass *8500
        !fALs = 2.623e23

        !endif
        !if ((rlength>100)) then
        !fALs = 2.62e23+(rho(1)*  abs(1.333*3.14*(abs(rlength-100)*1e5)**3) *(0.0169))
           !print*,'second'
        !endif
        !fALs = 2.53e23
        !Abundance of Fe (kg^-1)
        !if (rlength<=100) then

        fFes = rho(1)* mass *182800
        !fFes = 2.41e24

        !endif
        !if ((rlength>100)) then
    
        !    fFes =2.41e24+( rho(1)*  abs(1.333*3.14*(abs(rlength-100)*1e5)**3) *(0.24))
        !endif
        !
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
        
        !Creates a 1D array to store the temperature and H-values further down
        allocate(th(2))

        !Compute BUlk C (specific heat capacity) as a weighted average of the specific heat capacities of each phase
        bulkC = dot_product(p,c)
        
        !import t
        !For each element of q (i.e each time step), calculate the value of the heat source term
        

        !for each time step
       
        do nJ = 1,J-1

            !print*,'timestep =',nJ
            !Compute T at next time step
            do nN = 2,N-1
                !Note that our bulkk k-value is a 2D array which is updated at each accretion step
                if (differentiation == 0) then
                    q = heat(fALs, Altratio, EAl, LifeAl, fFes,Feratio,EFe,LifeFe, tac(count,nJ+1))
                    temp(nJ+1,nN) = ((2*bulkk(count,nN)*dt)/(rho(1)*bulkC*nN*(dr**2)))*(temp(nJ,nN+1)-temp(nJ,nN-1))+&
                    ((bulkk(count,nN)*dt)/(rho(1)*bulkC*(dr**2)))*(temp(nJ,nN+1)-2*temp(nJ,nN)+temp(nJ,nN-1))+&
                    temp(nJ,nN)+(dt/bulkC)*q
                    !Neumann boundary conditions
                    temp(nJ+1,1) = temp(nJ+1,2)
                endif
                if (differentiation==1) then
                    if (nN <= 200) then  
                        q = heat(0*fALs, Altratio, EAl, LifeAl, fFes,Feratio,EFe,LifeFe, tac(count,nJ+1))
                        temp(nJ+1,nN) = ((2*bulkk(count,nN)*dt)/(core_density*bulkC*nN*(dr**2)))*(temp(nJ,nN+1)-&
                        temp(nJ,nN-1))+((bulkk(count,nN)*&
                        dt)/(core_density*bulkC*(dr**2)))*(temp(nJ,nN+1)-2*temp(nJ,nN)+temp(nJ,nN-1))+temp(nJ,nN)+(dt/bulkC)*q
                    endif
                    if (nN>200)then
                        q = heat(fALs, Altratio, EAl, LifeAl, 0*fFes,Feratio,EFe,LifeFe, tac(count,nJ+1))
                        temp(nJ+1,nN) = ((2*bulkk(count,nN)*dt)/(mantle_density*bulkC*nN*(dr**2)))*(temp(nJ,nN+1)-&
                        temp(nJ,nN-1))+((bulkk(count,nN)*&
                        dt)/(mantle_density*bulkC*(dr**2)))*(temp(nJ,nN+1)-2*temp(nJ,nN)+temp(nJ,nN-1))+temp(nJ,nN)+(dt/bulkC)*q
                    endif
                    !Neumann boundary conditions
                    temp(nJ+1,1) = temp(nJ+1,2)
                    
                endif
            enddo

            !Computes residual heat of fusion values 
            do nN = 1,N

                do nW = 1,5

                    Hsil(nW,nN)=Hsil(nW,nN)-bulkC*(temp(nJ+1,nN)-M(1,nW))

                enddo

                Hmet(1,nN)=Hmet(1,nN)-bulkC*(temp(nJ+1,nN)-M(2,1))

                Hsulf(1,nN)=Hsulf(1,nN)-bulkC*(temp(nJ+1,nN)-M(3,1))

                do nW= 1,5

                    Hconj(nW,nN) = Hconj(nW,nN)-bulkC*(temp(nJ+1,nN)-M(4,nW))

                enddo

                !Following seciton uses an algorithm modified from Reynolds et al (1966)
                !Incorporates melting, algorithm is a subfunction called Renolds which is presented below
                !The output is a 1x2 matrix TH = [new temp, new residual heat of fusion]
                !Values from TH are then assigned to temp(nj+1,nN) and Hphase(w,n)

                !Silicates
              

                do nW = 1,5

                    th = reynolds(temp(nJ+1,nN),M(1,nw),Hsil(nW,nN),Hstart(1,nW),c(1),P(1))

                    temp(nJ+1,nN) = th(1)

                    Hmet(1,nN) = th(2)

                enddo

                !metals-only

                th = reynolds(temp(nJ+1,nN),M(2,1),Hmet(1,nN),Hstart(2,1),c(2),P(2))

                temp(nj+1,nN) = th(1)

                Hmet(1,nN) = th(2)

                !Sulfide-only

                th = reynolds(temp(nJ+1,nN),M(3,1),Hsulf(1,nN),Hstart(3,1),c(3),P(3))

                temp(nJ+1,nN) = th(1)

                Hsulf(1,nN) = th(2)

                !Conjoined grains

                do nW = 1,5

                    th = reynolds(temp(nJ+1,nN),M(4,nW),Hconj(nW,nN),Hstart(4,nW),c(4),P(4))

                    temp(nJ+1,nN) = th(1)

                    Hconj(nW,nN) = th(2)

                enddo   

                !Following section adjusts the value of thermal conductivity to account for decreasing pore space after partial silicate melting
                !If temperature is increasing 
                if (temp(nj+1,nN) > temp(nJ,nN)) then

                    !If accretion has finished
                    if (acc_con ==1) then    

                        !if the current space step (n) is within the regolith
                        if (nN > N - Reg) then 
                        !Do nothing

                        else 
                            !Check to see what ranking is required here
                            if(temp(nJ+1,nN) > M(1,5)) then
                                !If the current nN is NOT in the regolith
                                !If T exceeds silicate solidus
                                !Set specific heat capacity for this n to k(3)
                                if (nN > 50) then
                                    bulkk(50,:) = k(3) !Might need to swap ranking around
                                else 
                                    bulkk(nN,:) = k(3)
                                endif
                            !If T is below the silicate solidus
                            else

                                !Do nothing

                            endif

                        endif

                    !If accretion is continuing
                    !If T exceeds silicate solidus

                    else 

                        if (temp(nJ+1,nN) > M(1,5)) then
                            if (nN > 50) then
                                bulkk(50,:) = K(3)
                            else
                                bulkk(nN,:) = K(3)
                            endif

                        !If T is below the silicate solidus
                        else    

                        !DO nothing    

                        endif

                    endif
                !If temperature is not increasing

                else 
                    !DO nothing. this allows for a permanent thermal confuctivigy change after partial melting

                endif        

            enddo

        enddo
        !print*,'size of tT(1) is', SIZE(tac(count,:)), 'array size is', SIZE(tT(:,1))
        !print*,'size of tT(2) is', SIZE(temp)!, 'array size is', SIZE(tT(:,2:N+1))
            ! tT(1,1) = tac(count,1) - dt
            do i = 1, INT(SIZE(tac(count,:)))
                tT(i,1) = tac(count,i) !Might need to fix up
            enddo
        

        ! print*,'size of tT(i,2:N+1)', SIZE(tT(:,2:N)), 'array size', SIZE(temp(:,:))
     

       
        do i = 1,SIZE(temp(:,1))
            do j = 1, SIZE(temp(1,:))
                tT(i,j+1) = temp(i,j)
            enddo
        enddo

         
        allocate(thk(14,N))

        thk(1,:) = Temp(tlength,:)
  
        thk(2,:) = Hsil(1,:)
        thk(3,:) = Hsil(2,:)
        thk(4,:) = Hsil(3,:)
        thk(5,:) = Hsil(4,:)
        thk(6,:) = Hsil(5,:)

        thk(7,:) = Hmet(1,:)

        thk(8,:) = Hsulf(1,:)

        thk(9,:) = Hconj(1,:)
        thk(10,:) = Hconj(2,:)
        thk(11,:) = Hconj(3,:)
        thk(12,:) = Hconj(4,:)
        thk(13,:) = Hconj(5,:)
       
        do i = 1,SIZE(thk(14,:))
            thk(14,:) = bulkk(count,i)
        enddo

        print*,'dt =', dt 
        print*,' dr =', dr

        !Check to see if thk output is correct

        ! write(filename,"(a,i5.5,a)") 'thk',count,'.txt'
        ! print "(a)",' writing to '//trim(filename)
        ! open(newunit=iu,file=filename,status='replace',&
        ! action='write')
        ! write(iu,"(a)") '#  Thk'
        !     do i=1,SIZE(thk(:,1))
        !         write(iu,*) thk(i,:)!fmt='(10F15.2)') tT(i,:)
        !     enddo
        ! close(iu)

        ! !Check to see if tT output is correct
        ! write(filename,"(a,i5.5,a)") 'tt_',count,'.txt'
        ! print "(a)",' writing to '//trim(filename)
        ! open(newunit=iu,file=filename,status='replace',&
        ! action='write')
        ! write(iu,"(a)") '#  t, r'
        !     do i=1,SIZE(tT(:,1))
        !         write(iu,fmt='(10F15.2)') tT(i,:)
        !     enddo
        ! close(iu)
    end subroutine grad_a
    
end module grad
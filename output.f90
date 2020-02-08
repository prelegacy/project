module output
implicit none	
contains
    subroutine write_output(t,r,temp,dt)
        integer :: N,J,nJ,i,iu, outputfile,uu
        character(len=25) :: filename, other
        real, dimension(:), intent(inout) ::r,t,dt
        real, dimension(:,:),intent(out) :: temp
        N = SIZE(r)
        J = SIZE(t)

    print*,'How would you like the output?'
	print*,'(1) all values, plot in splash afterwards, there are ', J/100,' files to make'
	print*,'(2) max values at each radial value'
	print*,'(3) max values at each timestep'
	read*,outputfile
    select case(outputfile)
    case(1)
        open(newunit=uu,file='splash.titles',status='replace',&
            action='write')
        do nJ = 1,J,100
           if (nJ < 1e5) then
                write(uu, *) 't =',(nJ*dt(1))/1e6,'Myr'
                write(filename,"(a,i5.5,a)") 'output_0',nJ,'.dat'
           else if ((nJ > 1e5 ) .and. (nJ < 1e6)) then
                write(uu, *) 't =',(nJ*dt(1))/1e6,'Myr'
                write(filename,"(a,i6.5,a)") 'output_',nJ,'.dat'
            else if ((nJ > 1e6 ) .and. (nJ < 1e7)) then
                write(filename,"(a,i7.5,a)") 'output_',nJ,'.dat'
                write(uu, *) 't =',(nJ*dt(1))/1e6,'Myr'
            else if ((nJ > 1e7 ) .and. (nJ < 1e8)) then
                write(filename,"(a,i8.5,a)") 'output_',nJ,'.dat'
            else if ((nJ > 1e8 ) .and. (nJ < 1e9)) then
                write(filename,"(a,i9.5,a)") 'output_',nJ,'.dat'
            endif
            print "(a)",' writing to '//trim(filename)
            open(newunit=iu,file=filename,status='replace',&
            action='write')
            write(iu,"(a)") '#  r, temp'
                do i=1,N
                    write(iu,*) r(i),temp(Nj,i)
                enddo
                close(iu)
                write(filename,"(a,i5.5,a)") 'output_',nJ,'.dat'
        enddo

        
    case(2)       
        write(filename,"(a)")'maxoutput.dat'
        print "(a)",' writing to '//trim(filename)
        open(newunit=iu,file=filename,status='replace',&
        action='write')
        write(iu,"(a)") '#  r, temp'
        do i=1,N
            write(iu,*) r(i),maxval(temp(:,i))
        enddo
        close(iu)

    case(3)       
        write(filename,"(a)")'maxoutput.dat'
        print "(a)",' writing to '//trim(filename)
        open(newunit=iu,file=filename,status='replace',&
        action='write')
        write(iu,"(a)") '#  t, temp'
        do i=1,J
            write(iu,*) t(i),maxval(temp(i,:))
        enddo
        close(iu)

    ! call execute_command_line("python3 plot.py")
    end select    
           
           

        


        
    end subroutine write_output
end module output


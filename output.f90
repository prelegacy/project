module output
implicit none	
contains
    subroutine write_output(t,r,temp,dt,melting)
        integer :: N,J,nJ,i,iu, outputfile,uu,ans
        integer, intent(inout) :: melting
        character(len=25) :: filename, other
        character(*),parameter :: fileplace = "/home/ccnew1/project/output_files/"
        real, dimension(:), intent(inout) ::r,t,dt
        real, dimension(:,:),intent(out) :: temp
        N = SIZE(r)
        J = SIZE(t)

    print*,'How would you like the output?'
	print*,'(1) all values, plot in splash afterwards, there are ', J/100,' files to make'
	print*,'(2) max values at each radial value and timestep'
	! print*,'(3) max values at each timestep'
	read*,outputfile
    select case(outputfile)
    case(1)
        
        open(newunit=uu,file='splash.titles',status='replace',&
            action='write')
            if (melting == 1) then
                write(uu, *) 'Temp vs Radius plot for the regular H case'
            else if (melting == 2) then    
                write(uu, *) 'Temp vs Radius plot for the H case excluding metal sulfide melting'
            else if (melting == 3) then    
                write(uu, *) 'Temp vs Radius plot for the H case excluding conjoined melting'
            else if (melting == 4) then    
                write(uu, *) 'Temp vs Radius plot for the H case using alternate C'
            endif
        do nJ = 1,J,100
           if (nJ < 1e5) then
                write(filename,"(a,i5.5,a)") 'output_0',nJ,'.dat'
           else if ((nJ > 1e5 ) .and. (nJ < 1e6)) then
                write(filename,"(a,i6.5,a)") 'output_',nJ,'.dat'
            else if ((nJ > 1e6 ) .and. (nJ < 1e7)) then
                write(filename,"(a,i7.5,a)") 'output_',nJ,'.dat'
               
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
        print*,'basic info -----'
        print*,''
        print*,'Max T value = ',MAXVAL(temp)    
        print*,''
        print*,'files generated, would you like to run splash to view?'
        print*,'(1) Yes please!'
        print*,'(2) No thanks'
        read*,ans
        if (ans == 1) then
            call execute_command_line('splash *.dat')
        else
            !Do nothing
            print*,'program has ended'
        endif
       ! call execute_command_line("mv *.dat ~/project/output_files/")
        !call execute_command_line("cd ~/project/output_files/")
    case(2)       
        write(filename,"(a)")'rmaxoutput.dat'
        print "(a)",' writing to '//trim(filename)
        open(newunit=iu,file=filename,status='new',&
        action='write')
        write(iu,"(a)") '#  r, temp'
        do i=1,N
            write(iu,*) r(i),maxval(temp(:,i))
        enddo
        close(iu)

        write(filename,"(a)")'tmaxoutput.dat'
        print "(a)",' writing to '//trim(filename)
        open(newunit=iu,file=filename,status='new',&
        action='write')
        write(iu,"(a)") '#  t, temp'
        do i=1,J
            write(iu,*) t(i),maxval(temp(i,:))
        enddo
        close(iu)
        print*,'basic info -----'
        print*,''
        print*,'Max T value = ',MAXVAL(temp)    
        print*,''
        print*,'files generated, would you like to run python3 to view?'
        print*,'(1) Yes please!'
        print*,'(2) No thanks'
        read*,ans
        if (ans == 1) then
            call execute_command_line("python3 plot.py")
        else
            !Do nothing
            print*,'program has ended'
        endif
       
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

     call execute_command_line("python3 plot.py")
    end select    
           
           

        


        
    end subroutine write_output
end module output


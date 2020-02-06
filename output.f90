module output
implicit none	
contains
    subroutine write_output(t,r,temp)
        integer :: N,J,nJ,i,iu, outputfile
        character(len=25) :: filename, other
        real, dimension(:), intent(inout) ::r,t 
        real, dimension(:,:),intent(out) :: temp
        N = SIZE(r)
        J = SIZE(t)

    print*,'How would you like the output?'
	print*,'(1) all values, plot in splash afterwards'
	print*,'(2) max values at each radial value'
	print*,'(3) max values at each timestep'
	read*,outputfile
	select case(outputfile)
	case(1)
        do nJ = 1,J 
           
             write(filename,"(a,i5.5,a)") 'output_',nJ,'.dat'
            
            print "(a)",' writing to '//trim(filename)
            open(newunit=iu,file=filename,status='replace',&
            action='write')
            write(iu,"(a)") '#  t, r, temp'
                do i=1,N
                    write(iu,*) t(nJ),r(i),temp(Nj,i)
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


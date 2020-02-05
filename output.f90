module output
implicit none	
contains
    subroutine write_output(t,r,temp)
        integer :: N,J,nJ,i,iu, year,ii
        character(len=25) :: filename, other
        real, dimension(:), intent(inout) ::r,t 
        real, dimension(:,:),intent(out) :: temp
        N = SIZE(r)
        J = SIZE(t)

        
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


        
    end subroutine write_output
end module output


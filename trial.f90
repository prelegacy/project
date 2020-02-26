module trial
    implicit none
    
contains
    
    subroutine practice(tvals)

        real,allocatable,dimension(:),intent(in):: tvals
        integer ::nz, zval,i,nN,iu,nJ,ncounter
        character(len=25) :: filename
        print*,'shape of tac is', SHAPE(tvals)

        !Need to check out how the 





        write(filename,"(a)")'toutput.txt'
        print "(a)",' writing to '//trim(filename)
        open(newunit=iu,file=filename,status='replace',&
        action='write')
        write(iu,"(a)") '#  t'
        do i = 1,SIZE(tvals)
            write(iu,*)tvals(i)
        enddo
        ! do i=1,SIZE(tvals!SIZE(tac(1,:))
        !     ! print*,' i =',i, ' of ',SIZE(tT(1,:))
        !     write(iu,fmt='(50F15.2)') tvals(i,:) 
        ! enddo
        close(iu)


    end subroutine 
end module trial
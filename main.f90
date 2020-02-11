!A program that models the onion shell thermal model (Moskovitz & Gaidos 2011)

program main
use setup !(add in as required)
use functions
use heateq
use output
implicit none
	real, allocatable, dimension(:) :: k, rho, fal, ffe,c,a,r,t,dt,dr,P
	real, allocatable, dimension(:,:)::Hin,M,temp,Hsil,Hmet,Hsulf,Hconj,bulkk
	real :: al_ab, fe_ab,tau_al, tau_fe, E_al, E_fe, al, fe,trial,init,bdry,q,stab
	!integer :: n, ng, ntotal,nfile 
	!real :: t
	integer :: model, i,melting

	!original value for our output file
	!nfile = 0
	
	!initial t value
	! t=0
	print*, 'hello world'
	print*,'starting program'
	
	print*,'(1) Instantaneous Onion Shell model (Moskovitz & Gaidos 2011)?'
	print*,'(2)'
	print*,'(3)'
	read*,model
	select case(model)
	case(1)
		print*, 'starting setup'
		call setupinitial(k, rho, c, fal, ffe, al_ab, fe_ab, tau_al, tau_fe, E_al, E_fe,r,t,dt,dr,al,fe,P,init,bdry,Hin,M,melting)
		stab =stability(dt(1),dr(1))
		if (stab <=0.01) then
			print*,'setup complete'
			print*,'dr and dt values are ', dr(1),dt(1)
			print*,'with a stability value of ',stab
			print*, 'there will be ', SIZE(t),' timesteps'
			print*,'initializing heat equation'
			call heateqn(temp,r,t,dr,dt,init,bdry,Hin,Hsil,Hmet,Hsulf,Hconj,P,c,k,bulkk,q,fal, al, E_al, tau_al, ffe,fe,E_fe,tau_fe &
			,rho,M)
			print*,'heat eqn complete'
			print*,'initializing output'
			call write_output(t,r,temp,dt,melting)
		endif
		

	end select
	!case(1)
			! call setupinitial(pos, vel, mass, rho, h, pres, accel, n, ntotal,cs,ke)
			! print*,'starting derive'
			! call derivs(pos,vel,mass,rho,pres,accel,h,n,ng,cs)
			! print*,'writing output'
			! call write_output(pos, vel, mass, rho, h, pres, accel, n, ng,nfile,t,ke)
			! call tke_output(ke,vel,mass,t,n,ng)
			! call leap(pos,vel,mass,rho,pres,accel,h,cs,n,ng,ntotal,t,ke)
			! print*, 'we ended up with',ng,'ghosts'	

end program main


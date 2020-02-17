!A program that models the onion shell thermal model (Moskovitz & Gaidos 2011)

program main
use setup !(add in as required)
use functions
use heateq
use output
implicit none
	real, allocatable, dimension(:) :: k, rho, fal, ffe,c,a,r,t,dt,dr,P,Hstart_imp,rvals,tvals
	real, allocatable, dimension(:,:)::Hin,M,temp,Hsil,Hmet,Hsulf,Hconj,bulkk,Hstart,N,J,temps_time,rad,tac,delt,delx
	real :: al_ab, fe_ab,tau_al, tau_fe, E_al, E_fe, al, fe,trial,init,bdry,q,stab,final_rad,t_acc,t_dur,tfin
	!integer :: n, ng, ntotal,nfile 
	!real :: t
	integer :: model, i,melting, reg,Z,rstep_tot,tstep_dur,tstep_fin,tstep_tot

	!original value for our output file
	!nfile = 0
	
	!initial t value
	! t=0
	print*, 'hello world'
	print*,'starting program'
	
	print*,'(1) Instantaneous Accretion Onion Shell model (Moskovitz & Gaidos 2011)'
	print*,'(2) Gradual accretion model'
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
		
	case(2)
		call gradinitial(k,reg,rho,c,P,init,bdry,Hstart,Hstart_imp,M,Z,final_rad,rvals, rstep_tot,t_acc,t_dur,tfin,tvals &
		,tstep_dur,tstep_fin,N,J,tstep_tot,temps_time)
		call  heateqn_a(Z,rad,rvals,tac,tvals,tfin,tstep_fin,tstep_dur,delt,delx)
	end select
	
end program main


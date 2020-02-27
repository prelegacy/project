!A program that models the onion shell thermal model (Moskovitz & Gaidos 2011)

program main
use setup !(add in as required)
use functions
use heateq
use output
use trial
implicit none
	real, allocatable, dimension(:) :: k, rho, fal, ffe,c,a,r,t,dt,dr,P,Hstart_imp,rvals,tvals,tcounter,rcounter,delxx,deltt, &
	init_array
	! real(kind=8), allocatable,dimension(:)::
	real, allocatable, dimension(:,:)::Hin,M,temp,Hsil,Hmet,Hsulf,Hconj,bulkk,Hstart,N,J,temps_time,rad,tac,delt,delx,&
	tt,THK,tac_final
	real :: al_ab, fe_ab,tau_al, tau_fe, E_al, E_fe, al, fe,init,bdry,q,stab,final_rad,t_acc,t_dur,tfin,stab1, stab2
	!integer :: n, ng, ntotal,nfile 
	real :: number
	integer :: model, i,melting, reg,Z,rstep_tot,tstep_dur,tstep_fin,tstep_tot, acc_con

	!original value for our output file
	!nfile = 0
	
	!initial t value
	! t=0
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
        ,tstep_dur,tstep_fin,N,J,tstep_tot,temps_time,rad,tac,delt,delx,delxx,tcounter,rcounter,deltt,tac_final)
		!Determine the stability of the program


		stab1 = stability(real((tac(1,2)-tac(1,1))),1000.0)
        stab2 = stability(7000.0,1000.0)
		
        !Only continue on if stab1 and stab2 are less than 0.01
		if (stab1 < 0.01 .and. stab2 < 0.01) then
			print*, "Results will be stable with stability values of"
			print*, stab1, 'and ', stab2
			print*, "Starting heateq_a"
		else 
			print*,'results are not stable, try to improve your dt values'
			call exit()
		endif

		! call practice(tvals)
		
		call heateqn_a(k,Z,rad,reg,tac,deltt,delxx,temp,init,bdry,c,p,Hin,Hstart_imp,init_array,acc_con,rho,tT,temps_time,&
        bulkk,THK,m,Hstart,tstep_dur,tac_final,tstep_fin,tstep_tot)


		!(k,Z,rad,reg,tac,deltt,delxx,temp,init,bdry,c,p,Hin,Hstart_imp,init_array,acc_con,rho,tT,temps_time,&
        !bulkk,THK,m,Hstart,tstep_dur,tac_final,tstep_fin)
		!Create the write outputs
		
	end select
	
end program main


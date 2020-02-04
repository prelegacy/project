module functions
    implicit none
    contains    

    !Determines the stability of the finite difference method
    function stability(dt,dr) 
        integer, intent(in) :: dt,dr
        real:: stability
        
        stability = dt/(dr**2)
        
    end function stability

    !Determiens the heating term q
    function heat(fal, al, E_al, tau_al, ffe,fe,E_fe,tau_fe, t)
        real, intent(in):: fal, al, E_al, tau_al, ffe,fe,E_fe,tau_fe, t
        real:: heat

        heat = (fal*al*E_al/tau_al)*exp(-t/tau_al)+(ffe*fe*E_fe/tau_fe)*exp(-t/tau_fe)

    end function heat

    function reynolds(temp,M,H,Hin,bulkc,P)
        real,dimension(2):: reynolds
        real,intent(in) :: temp,M,Hin,bulkc,P
        real:: T2,T,H

        T = temp     

        if ( temp < M ) then
            if ( H>=Hin ) then
                H = Hin 
            else
                T2 = T+P*((Hin-H)/bulkc)
                if ( T2 < M ) then
                    T = T2
                    H = Hin
                else
                    T = M
                end if
            end if
        else
            if ( H <= 0 ) then
                H = 0
            else
                T2 = T - P*((H/bulkc))
                if ( T2 > M ) then
                    T = T2
                    H = 0
                else
                    T = M
                end if
            end if
        end if
    reynolds = (/T,H/)
    end function reynolds
end module functions
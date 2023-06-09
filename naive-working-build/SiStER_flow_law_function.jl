function SiStER_flow_law_function(type,pre,Ea,n,R,T,epsII,sII,PARAMS)
# [eta] = SiStER_flow_law_function(type,pre,Ea,n,R,T,epsII,sII)
#  gives generic form of ductile creep law 
#  to be used for both diffusion & dislocation creep in
#  get_ductile_rheology functions()
#  inputs:
#  type = "from_stress' | 'from_strain_rate" to set flow law expression
#  with usual parameters: prefactor [pre], activation energy [Ea, J/mol],
#  exponent [n], gas constant [R], temperature [T, in deg C], 
#  second invariant of strain rate [epsII, s^-1] | second invariant
#  of deviatoric stress [sII, Pa]

if type == "from_stress"
    
    eta = (pre.^-1).*sII.^(1 .-n).*exp.(Ea./(R.*(T.+273.15)));
    
    
elseif type == "from_strain_rate"
    
    eta = pre.^(-1 ./n).*epsII.^((1-n)./n).*exp((Ea)./(n.*R.*(T.+273.15)));
    
    
elseif type == "custom"
    
    disp("ERROR: CUSTOM VISCOSITY FUNCTION NOT CURRENTLY DEFINED")
    disp("this feature will be available in a future update.")
    eta = PARAMS.customviscofunction[phase,temperature,strain_rate,stress];
    
else

    disp("ERROR ? flow law is undefined.")

end

return eta

end
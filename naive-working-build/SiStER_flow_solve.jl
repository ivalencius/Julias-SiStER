# SiStER_flow_solve
# Performs inner solve of linear LS=R system as well as outer; iterative
# solution for non-linear dependence of L [viscosity] on S [vx,vy,P]
# Used to be named "run_Picard_iterations" but name changed by G.Ito 6/21/16
using LinearAlgebra

if PARAMS.BalanceStickyLayer==1
# BALANCE FLUXES ### JAO July 16; 2015
# RE-ADJUST BCs SO FLUX OF ROCK AND STICKY AIR MATERIAL ARE CONSERVED
# locate height of sticky layer - rock contact on the sides
    # bL=interp1(topo_x,topo_y,0);
    # bR=interp1(topo_x,topo_y,xsize);
    # bL = interpolate((topo_x,), topo_y, Gridded(Constant()))(0);
    # bR = interpolate((topo_x,), topo_y, Gridded(Constant()))(xsize);
    bL = linear_interpolation(topo_x, topo_y)(0);
    bR = linear_interpolation(topo_x, topo_y)(xsize);
    utop=BC.right[3]*(bL+bR)/xsize;
    ubot=BC.right[3]*(2*ysize-bL-bR)/xsize;
    BC.top[3]=utop;
    BC.bot[3]=-ubot;
end


ResL2=1; 

global pit = 1

for pittemp = 1:PARAMS.Npicard_max
    
    # SiStER_VEP_rheology cannot access global variable pittemp inside loop
    global pit = pittemp 
    # print(pit)
    
    if pit==1
        ResL2init=ResL2;
    end

    ## ---------------------------------------------------------------------------------
    # Compute visco-elasto-plastic viscosities
    #---------------------------------------------------------------------------------
    include("SiStER_VEP_rheology.jl");

    #---------------------------------------------------------------------------------
    # Assemble L & R matrices
    #---------------------------------------------------------------------------------
    include("SiStER_assemble_L_R.jl")
    L, R, Kc, Kb = SiStER_assemble_L_R(dx,dy,Zs.*etas,Zn.*etan,rho,BC,PARAMS,srhs_xx,srhs_xy); #G.Ito
    
    #---------------------------------------------------------------------------------
    # Residual:  L & R are from current solution S
    #---------------------------------------------------------------------------------
    if (@isdefined(S));
        Res=L*S-R;
        global ResL2=norm(Res,2)/norm(R,2);
    else
        global S=L\R; 
        print("First solve is linearized");
        Res=L*S-R;
        global ResL2=norm(Res,2)/norm(R,2);        
    end
    #---------------------------------------------------------------------------------
    # Solve for new solution S using Picard or approximate Newton | a
    # combination of the two 
    #---------------------------------------------------------------------------------  
    
    if (pit >= PARAMS.pitswitch);
        if pit==PARAMS.pitswitch; 
            print("switching from Picard to approx. Newton"); 
        end;
        beta=1;
        S=S-beta.*(L\Res);  # approximate Newton update, with L as approximation to Jacobian
        it_type="Newton: ";
    else
        S=L\R; # Picard update
        it_type="Picard: ";
    end
   
    include("SiStER_reshape_solver_output.jl")
    ptemp, vxtemp, vytemp = SiStER_reshape_solver_output(S,Kc,Nx,Ny);
    
    # Need to assign to the globals
    global p = ptemp
    global vx = vxtemp
    global vy = vytemp 

    ## ASSESS CONVERGENCE    
    if ResL2<PARAMS.conv_crit_ResL2 && pit >= PARAMS.Npicard_min
        print("\nFinal residual = "*string(ResL2))
        print("\n"*string(pit)*" iterations converged: L2 norm of residual dropped below "*string(PARAMS.conv_crit_ResL2));
        break;
	elseif pit==PARAMS.Npicard_max; #&& !(ResL2<PARAMS.conv_crit_ResL2 && pit >= PARAMS.Npicard_min)
        print("\nFinal residual = "*string(ResL2))
        print("\nWARNING! "*string(pit)*" Picard / approx. Newton iterations failed to converge within tolerance of "*string(PARAMS.conv_crit_ResL2));
    end
   
    ## get strain rate on nodes current solution
    include("SiStER_get_strain_rate.jl")
    EXXtemp, EXYtemp=SiStER_get_strain_rate(vx,vy,dx,dy,BC);

    # Need to assign to globals
    global EXX = EXXtemp
    global EXY = EXYtemp

    EXY_n=SiStER_interp_shear_to_normal_nodes(EXY);       
    EXX_s=SiStER_interp_normal_to_shear_nodes(EXX,dx,dy); 
    global epsII_n=sqrt.((EXX.^2)+(EXY_n.^2));
    global epsII_s=sqrt.((EXX_s.^2)+(EXY.^2));

    # helpful to visualize convergence
    # figure(1)
    # pcolor(X,Y,log10(etas))
    # set(gca,"ydir','reverse")
    # axis equal
    # caxis([18 25])
    # colorbar()
    # title(num2str(pit))
    # pause(.001)

    # RESIDUAL FOR INDIVIDUAL VARIABLES
    # [pres, vxres, vyres]=SiStER_reshape_solver_output[Res,Kc,Nx,Ny];
    # figure(1)
    # pcolor(X,Y,vxres)
    # set(gca,"ydir','reverse")
    # axis equal
    # colorbar()
    # title(num2str(pit))
    # pause(.001)

end
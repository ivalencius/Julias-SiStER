# SiStER_MAIN.m
#
# Simple Stokes solver with Exotic Rheologies
#
# Main routine doing initialization; time loop & outputs
#
#
# J.-A. Olive; B.Z. Klein; E. Mittelstaedt; M. Behn; G. Ito; S. Howell
# jaolive <at> ldeo.columbia.edu
# March 2011 - April 2017
using Statistics
using LinearAlgebra
using Plots
pyplot()

# All inclusions
include("SiStER_initialize_grid.jl")
include("SiStER_initialize_marker_positions.jl")
include("SiStER_locate_markers_in_grid.jl")
include("SiStER_initialize_marker_phases.jl")
include("SiStER_interp_shear_nodes_to_markers.jl")
include("SiStER_interp_phases_to_normal_nodes.jl")
include("SiStER_interp_phases_to_shear_nodes.jl")
include("SiStER_get_density.jl")
include("SiStER_get_cohesion.jl")
include("SiStER_get_friction.jl")
include("SiStER_interp_markers_to_shear_nodes.jl")
include("SiStER_get_elastic_moduli.jl")
include("SiStER_interp_markers_to_normal_nodes.jl")
include("SiStER_interp_normal_to_shear_nodes.jl")
include("SiStER_interp_shear_to_normal_nodes.jl")
include("SiStER_get_ductile_rheology_on_nodes_from_stresses.jl")
include("SiStER_assemble_L_R.jl")
include("SiStER_reshape_solver_output.jl")
include("SiStER_get_strain_rate.jl")
include("SiSter_interp_normal_nodes_to_markers.jl")
include("SiSter_interp_shear_nodes_to_markers.jl")
include("SiStER_set_timestep.jl")
include("SiStER_get_rotation_rate.jl")
include("SiStER_thermal_solver_sparse_CFD.jl")
include("SiStER_advect_markers.jl")
include("SiStER_locate_markers_in_grid.jl")
include("SiStER_patch_marker_holes.jl")
include("SiStER_topography_diffusion_solver.jl")

# INITIALIZATION
InpFil = "SiStER_Input_File_continental_rift.jl"
include(InpFil)

# include("SiStER_Help.jl")
# using .SiStER_Help
# names(SiStER_Help)

# Input File: loads parameter values; model geometry; boundary conditions
# if exist("running_from_SiStER_RUN','var")==0
#     clear() 
#     InpFil = input("Input file ? ','s");
# end
# run(InpFil)

# Make folder for figures
mkpath("figures/")

# construct grid & initialize marker / node arrays
println("---------------")
println("INITIALIZING")
println("---------------")

# include("SiStER_Initialize.jl");
# SiStER Initialize
# using Statistics

PARAMS.Nphase = Nphase; # for convenience

# construct staggered grids
# include("SiStER_initialize_grid.jl")
# @time = 0.000197 sec (49 allocations: 77.375 KiB)
(X,Y,x,y,xc,yc,dx,dy,Nx,Ny) = SiStER_initialize_grid(xsize,ysize,GRID);

# initialize marker arrays & positions
# @time = 33.94 sec (2.99 M allocations: 260.983 GiB, 30.25% fc time, 0.78% compilation time)
# New: @time = 0.148021 seconds (905.55 k allocations: 98.598 MiB, 4.31% gc time, 66.22% compilation time)
(xm, ym) = SiStER_initialize_marker_positions(xsize,ysize,dx,dy,Mquad);
print("**Number of markers: "*string(length(xm))*"**\n");

# locate markers with respect to grid()
# @time = 1.27 sec (4.35 M allocations: 367.611 MiB, 4,45% gc time, 93.55% compilation time)
(qd,icn,jcn) = SiStER_locate_markers_in_grid(xm,ym,x,y,dx,dy);
# icn = Int.(icn);
# jcn = Int.(jcn);

# assign marker phases
# @time = 0.192908 seconds (881.20 k allocations: 73.391 MiB, 4.83% gc time, 97.98% compilation time)
im = SiStER_initialize_marker_phases(Nphase,GEOM,xm,ym);
# im = im[1]; # returned as vector of a vector for some reason
# initialize marker plastic strain [to zero] & strain rate [to one]
ep=zeros(size(xm));
epNH=zeros(size(xm));
epsIIm=ones(size(xm));

# initialize marker stresses
sxxm=zeros(size(xm));
sxym=zeros(size(xm));

# initialize marker index [a unique number to identify & track each marker]
idm=1:length(xm);

# initialize temperature structure on nodes
T=PARAMS.a0.+(PARAMS.a1 .* Y) .+ (PARAMS.a2.*Y).^2 .+ (PARAMS.a3.*Y.^3);
T=T.+PARAMS.amp*sin.((2*pi).* X./PARAMS.lam);
if PARAMS.ynTreset==1 # reset T=T0 in top layer
    T[T.<PARAMS.T0].=PARAMS.T0;
end

# pass initial nodal T to markers
# @time = 0.504564 seconds (1.58 M allocations: 329.309 MiB, 6.43% gc time, 89.27% compilation time)
Tm = SiStER_interp_shear_nodes_to_markers(T,x,y,xm,ym,icn,jcn);
Tm = Tm';
Tm0 = copy(Tm);

# initialize nodal strain rate & other useful arrays
EXX=zeros(size(X));
EXY=zeros(size(X));
vx=zeros(size(X));
vy=zeros(size(X));
v=vx;
p=1e12.*ones(size(EXX));  #initialize to be high so plasticity doesnt activate at t=1, pit=1;
etan_new=zeros(Ny,Nx);
#-------------------------------------------------------------------------
# initialize dt_m small to keep things elastic & no plasticity at t=1; G.Ito
#-------------------------------------------------------------------------
# if (exist("dt_m','var")==0);
dt_m=1e2;
# end

# initialize marker chain to track base of layer 1 [sticky layer]
Ntopo=PARAMS.Ntopo_markers;
topo_x=LinRange(0, xsize, Ntopo);
topo_y=GEOM[1].bot.*ones(size(topo_x));
topo_marker_spacing=mean(diff(topo_x)); # initial mean spacing of topography markers

println("For continental drift: Nx = 106, NY=43 (from MATLAB solution)")
println("Mismatch due to round() rounding 4.5 to nearest even integer, not up to 5 like in MATLAB")
println("SiStER_initialize_grid.jl (line 40)")

# Save variables for test
# using JLD2
# @save "function_tests/test-vars.jld2" x y xm ym icn jcn im qd

### BEGIN TIME LOOP ############################################################

time=0;

# Start solve
for t = 1:Nt
# for t = 1:3

    timing = @elapsed begin
    println("STARTING ITERATION: " * string(t) * " out of " * string(Nt))
    
    # update time
    global time=time+dt_m;
    
    # Here we prepare nodal arrays to feed the Stokes solver
    ### Start of SiStER_material_props_on_nodes ################################

    # Get material properties on nodes from advected properties of markers
    # G.Ito 8/16 - replaces previous version; which relied on marker viscosity
    # for Picard iterations [J.-A.O.]
    # PHASE PROPORTIONS AT NORMAL AND SHEAR NODES. G.Ito 8/16

    # Two culprits of slowdown: SiStER_interp_phases_to_normal_nodes and SiStER_interp_phases_to_shear_nodes
    # @time = 1.407799 seconds (15.12 M allocations: 1.339 GiB, 11.47% gc time, 49.94% compilation time: 3% of which was recompilation)
    phase_n = SiStER_interp_phases_to_normal_nodes(xm,ym,icn,jcn,x,y,im, Nphase);
    # @time = 7.894677 seconds (36.21 M allocations: 1.877 GiB, 3.89% gc time, 81.22% compilation time)
    phase_s = SiStER_interp_phases_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,im, Nphase);

    # phase_n & _s is a Ny*Nx*Nphase array containing the proportion
    # of each phase at each node - this gets used in get_ductile_rheology
    # functions()

    # GET MARKER DENSITIES
    # @time = 0.088335 seconds (315.10 k allocations: 33.626 MiB, 95.79% compilation time)
    rhom = SiStER_get_density(im,Tm,MAT);

    # pass density to nodes
    # @time = 6.058524 seconds (30.60 M allocations: 1.600 GiB, 3.12% gc time, 79.76% compilation time)
    n2interp = SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,rhom);
    rho  = n2interp[1];


    # GET MARKER ELASTIC PROPERTIES  G.Ito 8/16
    # @time = 0.041856 seconds (35.91 k allocations: 10.643 MiB, 92.30% compilation time)
    Gm = SiStER_get_elastic_moduli(im,MAT);
    # pass shear modulus to nodes
    # @time = 1.617598 seconds (25.98 M allocations: 2.372 GiB, 20.31% gc time, 15.05% compilation time)
    n2interp = SiStER_interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,1 ./ Gm);
    Gn=1 ./ (n2interp[1]);
    # @time = 1.204658 seconds (25.46 M allocations: 1.303 GiB, 10.30% gc time)
    n2interp = SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,1 ./ Gm);
    Gs = 1 ./ (n2interp[1]);   

    # PROPERTIES FOR PLASTICITY  G.Ito 8/16
    # @time = 0.070012 seconds (179.41 k allocations: 24.294 MiB, 94.92% compilation time)
    cohes=SiStER_get_cohesion(im, ep, MAT); # cohesion depends on plastic strain
    n2interp = SiStER_interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,cohes);
    Cohes_n=n2interp[1];
    n2interp = SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,cohes);
    Cohes_s = n2interp[1];  

    # GET FRICTION BASED ON MARKERS J.A. Olive 4/17
    # @time = 0.052289 seconds (26.15 k allocations: 13.962 MiB, 92.95% compilation time)
    fric=SiStER_get_friction(im,ep,MAT); # friction depends on plastic strain
    n2interp = SiStER_interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,fric);
    Mu_n=n2interp[1];
    n2interp = SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,fric);
    Mu_s = n2interp[1]; 

    # ADVECTED strain rate invariant G.Ito 8/16
    n2interp = SiStER_interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,epsIIm);
    epsII_n=n2interp[1];
    n2interp = SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,epsIIm);
    epsII_s = n2interp[1];  


    # OLD STRESSES AND PRESSURES G.Ito 8/16
    n2interp = SiStER_interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,sxxm);
    sxxOLD = n2interp[1];
    sxxOLD_s=SiStER_interp_normal_to_shear_nodes(sxxOLD,dx,dy);

    n2interp = SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,sxym);
    sxyOLD = n2interp[1];
    sxyOLD_n = SiStER_interp_shear_to_normal_nodes(sxyOLD);

    #MIGHT WANT TO ADVECT PRESSURES [FOR SPEED?] G.Ito 8/16
    pold=p; 
    ps_old=SiStER_interp_normal_to_shear_nodes(p,dx,dy);

    EXYOLD=EXY;
    EXXOLD=EXX;
    EXX_sOLD=SiStER_interp_normal_to_shear_nodes(EXX,dx,dy);
    EXY_nOLD=SiStER_interp_shear_to_normal_nodes(EXY);

    #TEMPERATURE ARRAYS NEEDED FOR DUCTILE RHEOLOGY  G.Ito 8/16
    n2interp=SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,Tm);
    Ts=n2interp[1];
    n2interp=SiStER_interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,Tm);
    Tn=n2interp[1];

    ### End of SiStER_material_props_on_nodes ##################################

    ### SOLVE STOKES WITH NON-LINEAR RHEOLOGY HERE 
    ### Start of SiStER_flow_solve #############################################
    # SiStER_flow_solve
    # Performs inner solve of linear LS=R system as well as outer; iterative
    # solution for non-linear dependence of L [viscosity] on S [vx,vy,P]
    # Used to be named "run_Picard_iterations" but name changed by G.Ito 6/21/16
    if PARAMS.BalanceStickyLayer==1
    # BALANCE FLUXES ### JAO July 16; 2015
    # RE-ADJUST BCs SO FLUX OF ROCK AND STICKY AIR MATERIAL ARE CONSERVED
    # locate height of sticky layer - rock contact on the sides
        bL = linear_interpolation(topo_x, topo_y)(0);
        bR = linear_interpolation(topo_x, topo_y)(xsize);
        utop=BC.right[3]*(bL+bR)/xsize;
        ubot=BC.right[3]*(2*ysize-bL-bR)/xsize;
        BC.top[3]=utop;
        BC.bot[3]=-ubot;
    end


    ResL2=1; 

    # @time = 6.955675 seconds (8.18 M allocations: 6.318 GiB, 5.02% gc time, 37.41% compilation time)
    for pit = 1:PARAMS.Npicard_max
        
        if pit==1
            ResL2init=ResL2;
        end

        ## ---------------------------------------------------------------------------------
        # Compute visco-elasto-plastic viscosities
        #---------------------------------------------------------------------------------
        # SiStER_UPDATE_RHEOLOGY
        # calculates visco-elasto-plastic rheology terms on shear & normal nodes
        # G.Ito 8/20


        #--------------------------------------------------------------------------
        # plastic viscosity [Mohr-Coulomb]
        #--------------------------------------------------------------------------
        if PARAMS.YNPlas==1
            dp=p.-pold;
            ps=ps_old.+SiStER_interp_normal_to_shear_nodes(dp,dx,dy);  #<<<< Interpolation done here
            ps[ps.<0] .= 0;
            pnn = copy(p); 
            pnn[p.<0] .= 0;
            global yield_s=(Cohes_s+(Mu_s.*ps)).*cos.(atan.(Mu_s));
            global yield_n=(Cohes_n+(Mu_n.*pnn)).*cos.(atan.(Mu_n));
            if (PARAMS.YNElast==1) # elastic strain rate needs to be removed from total strain
                eta_plas_s=0.5.*yield_s./max.(epsII_s-(yield_s-sqrt.(sxxOLD_s.^2+sxyOLD.^2))./(2 .*Gs.*dt_m), minimum.(epsII_s).*1e-6);  
                eta_plas_n=0.5.*yield_n./max.(epsII_n-(yield_n-sqrt.(sxxOLD.^2+sxyOLD_n.^2))./(2. *Gn.*dt_m), minimum.(epsII_n).*1e-6); 
            else # elasticity is off; no need to remove elastic strain rate # JAO 04/17
                eta_plas_s=0.5.*yield_s./max(epsII_s,minimum(epsII_s)*1e-6);  
                eta_plas_n=0.5.*yield_n./max(epsII_n,minimum(epsII_n)*1e-6); 
            end
        end

        #-------------------------------------------------------------------------
        # Ductile viscosity
        #-------------------------------------------------------------------------


        # previous approach: based on strain rate 
        # problem = it should be only the viscous part of strain rate
        # Gerya [2010] p. 189 - best to base viscosity on stress
        #[etas_new]=SiStER_get_ductile_rheology[MAT,PARAMS,Ts,epsII_s,phase_s] 
        #[etan_new[2:end,2:end]]=SiStER_get_ductile_rheology[MAT,PARAMS,Tn[2:end,2:end],epsII_n[2:end,2:end],phase_n[2:end,2:end]];


        # VISCOSITY INDEXED ON STRESS
        if PARAMS.YNElast==1

            if t==1 & pit==1
                # at this stage we don't have any viscosity
                # but stresses are zero anyway
                sII_n=zeros(size(X));
                sII_s=zeros(size(X));
            else    
                # need sIIOLD_s to initialize bisection [done on shear nodes]
                sIIOLD_s=sqrt.(sxxOLD_s.^2+sxyOLD.^2);
                # BISECTION ALGORITHM FOR LOCAL ITERATIONS
                # calculates a deviatoric stress sII that can be fed to the creep laws
                # & is consistent with the current viscosity
                # this appears to be a necessary step to use a stress-based creep law
                # when elasticity is on [i.e., the strain rate is not entirely viscous]

                # STARTING VISCOSITY [TO BE ADJUSTED]
                visco=copy(etas_new);

                # INITIALIZE
                sIILOW=copy(sIIOLD_s);
                sIIHIGH=copy(sIIOLD_s);
                visco0=SiStER_get_ductile_rheology_on_nodes_from_stresses(MAT,PARAMS,Ts,sIIOLD_s,phase_s);
                if PARAMS.YNElast==1
                    Zs=Gs.*dt_m./(Gs.*dt_m+visco0);
                else
                    Zs=ones(Ny,Nx);
                end
                sIINEW=2 .*visco0.*epsII_s.*Zs+((1 .-Zs).*sIIOLD_s);
                sIIHIGH[sIINEW.>sIIHIGH]=sIINEW[sIINEW.>sIIHIGH];
                sIILOW[sIINEW.<sIILOW]=sIINEW[sIINEW.<sIILOW];

                viscoNEW=SiStER_get_ductile_rheology_on_nodes_from_stresses(MAT,PARAMS,Ts,sIINEW,phase_s);
                if PARAMS.YNElast==1
                    Zs=Gs.*dt_m./(Gs.*dt_m+viscoNEW);
                else
                    Zs=ones(Ny,Nx);
                end 
                sIINEW=2 .*viscoNEW.*epsII_s.*Zs+((1 .-Zs).*sIIOLD_s);
                sIIHIGH[sIINEW.>sIIHIGH]=sIINEW[sIINEW.>sIIHIGH];
                sIILOW[sIINEW.<sIILOW]=sIINEW[sIINEW.<sIILOW];

                # 10 iterations of a BISECTION ALGORITHM
                Nbisec=10;
                global sII = 0
                for i=1:Nbisec
                    
                    global sII=0.5 .*(sIIHIGH+sIILOW);
                    sIIPAST = copy(sII);
                    visco=SiStER_get_ductile_rheology_on_nodes_from_stresses(MAT,PARAMS,Ts,sII,phase_s);
                    if PARAMS.YNElast==1
                        Zs=Gs.*dt_m./(Gs.*dt_m+visco);
                    else
                        Zs=ones(Ny,Nx);
                    end
                    global sII=2 .*visco.*epsII_s.*Zs+((1 .-Zs).*sIIOLD_s);

                    sIILOW[sII.>sIIPAST]=sIIPAST[sII.>sIIPAST];
                    sIIHIGH[sII.<=sIIPAST]=sIIPAST[sII.<=sIIPAST];

                end

                global sII_s=copy(sII);
                global sII_n=SiStER_interp_shear_to_normal_nodes(sII_s);   
            end

            # & use the stresses to re-update viscosity
            global etas_new=SiStER_get_ductile_rheology_on_nodes_from_stresses(MAT,PARAMS,Ts,sII_s,phase_s);
            global etan_new[2:end,2:end]=SiStER_get_ductile_rheology_on_nodes_from_stresses(MAT,PARAMS,Tn[2:end,2:end],sII_n[2:end,2:end],phase_n[2:end,2:end,:]);
            
        else # if elasticity is off; the strain rate is entirely viscous; so we can use a strain rate-based viscosity law
            # this seems to yield an easier convergence than the bisection
            # algorithm in this case()
            global etas_new=SiStER_get_ductile_rheology_on_nodes_from_strain_rate(MAT,PARAMS,Ts,epsII_s,phase_s);
            global etan_new[2:end,2:end]=SiStER_get_ductile_rheology_on_nodes_from_strain_rate(MAT,PARAMS,Tn[2:end,2:end],epsII_n[2:end,2:end],phase_n[2:end,2:end,:]);
            
        end 

        if PARAMS.YNPlas==1 
            # identify yielding nodes to update ep
            global s_nodes_yield=findall(etas_new.>eta_plas_s);
            # n_nodes_yield=findall(etan_new.>eta_plas_n);
            # incorporate plastic viscosity into effective viscosity
            etan_new[2:end,2:end]=((1 ./eta_plas_n[2:end,2:end]) + (1 ./etan_new[2:end,2:end]) .+ (1 ./PARAMS.etamax)).^-1;
            etas_new=((1 ./eta_plas_s) + (1 ./etas_new) .+ (1 ./PARAMS.etamax)).^-1;
        else
            etan_new[2:end,2:end]=((1 ./etan_new[2:end,2:end]) .+ (1 ./PARAMS.etamax)).^-1;
            etas_new=((1 ./etas_new) .+ (1 ./PARAMS.etamax)).^-1;
        end

        global etan = copy(etan_new);
        global etas = copy(etas_new);
        etan[etan.<PARAMS.etamin].=PARAMS.etamin;
        etas[etas.<PARAMS.etamin].=PARAMS.etamin;


        #-------------------------------------------------------------------------
        # ELASTICITY TERMS [BASED ON LATEST VISCOSITY]
        #-------------------------------------------------------------------------
        if PARAMS.YNElast==0
            global Zs=ones(size(etas));
            global Zn=ones(size(etan));
        else
            global Zs=Gs*dt_m./(Gs*dt_m+etas);
            global Zn=Gn*dt_m./(Gn*dt_m+etan);
        end
        # right-hand side [1-Z]*sigmaOLD_ij
        srhs_xx=(1 .-Zn).*sxxOLD;
        srhs_xy=(1 .-Zs).*sxyOLD;

        #---------------------------------------------------------------------------------
        # Assemble L & R matrices
        #---------------------------------------------------------------------------------
        L, R, Kc, Kb = SiStER_assemble_L_R(dx,dy,Zs.*etas,Zn.*etan,rho,BC,PARAMS,srhs_xx,srhs_xy); #G.Ito
        
        #---------------------------------------------------------------------------------
        # Residual:  L & R are from current solution S
        #---------------------------------------------------------------------------------
        # if (@isdefined(S));
        #     Res=L*S-R;
        #     ResL2=norm(Res,2)/norm(R,2);
        if pit == 1
            global S=L\R; 
            println("First solve is linearized");
            Res=L*S-R;
            ResL2=norm(Res,2)/norm(R,2);    
        else
            Res=L*S-R;
            ResL2=norm(Res,2)/norm(R,2);    
        end
        #---------------------------------------------------------------------------------
        # Solve for new solution S using Picard or approximate Newton | a
        # combination of the two 
        #---------------------------------------------------------------------------------  
        
        if (pit >= PARAMS.pitswitch);
            if pit==PARAMS.pitswitch; 
                println("switching from Picard to approx. Newton"); 
            end;
            beta=1;
            S=S-beta.*(L\Res);  # approximate Newton update, with L as approximation to Jacobian
            it_type="Newton: ";
        else
            S=L\R; # Picard update
            it_type="Picard: ";
        end
    
        ptemp, vxtemp, vytemp = SiStER_reshape_solver_output(S,Kc,Nx,Ny);
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
    
    # end
        ## get strain rate on nodes current solutio
        EXXtemp, EXYtemp = SiStER_get_strain_rate(vx,vy,dx,dy,BC);
        global EXX = EXXtemp
        global EXY = EXYtemp

        global EXY_n=SiStER_interp_shear_to_normal_nodes(EXY);       
        global EXX_s=SiStER_interp_normal_to_shear_nodes(EXX,dx,dy); 
        epsII_n=sqrt.((EXX.^2)+(EXY_n.^2));
        epsII_s=sqrt.((EXX_s.^2)+(EXY.^2));
    end

    ### End of SiStER_flow_solve #################$#############################
    
    # GET STRAIN RATE FROM CURRENT SOLUTION
    global epsIIm = SiStER_interp_shear_nodes_to_markers(epsII_s,x,y,xm,ym,icn,jcn);
    
    # USE STRAIN RATE TO UPDATE STRESSES ON MARKERS
    ### Start of Sister_update_marker_stresses #################################
    # Updates stresses on markers for CURRENT solution.  Stress rotation
    # occurs after solutions are output
    # G.Ito 8/16

    # Compute STRESS Changes on nodes; interpolate to markers; & apply to marker stresses
    dsxx=(2*etan.*EXX-sxxOLD).*Zn;
    dsxy=(2*etas.*EXY-sxyOLD).*Zs;

    dsxxm=SiStER_interp_normal_nodes_to_markers(dsxx,xc,yc,xm,ym,icn,jcn);
    dsxym=SiStER_interp_shear_nodes_to_markers(dsxy,x,y,xm,ym,icn,jcn);
    global sxxm=sxxm+dsxxm';
    global sxym=sxym+dsxym';
    ### End of Sister_update_marker_stresses ###################################
    
    # BUILD UP PLASTIC STRAIN IN YIELDING AREAS IF PLASTICITY IS ACTIVATED
    # @time = 0.218199 seconds (515.86 k allocations: 272.884 MiB, 9.46% gc time, 74.46% compilation time)
    if (PARAMS.YNPlas==1)
        ### Start of SiStER_update_ep ##########################################
        # PLASTIC [EXCLUDING ELASTIC] STRAIN ACCUMULATION
        # Not sure whether to use the actual non-elastic strain (i.e., with
        # current stresses) | the approximate; in which stress is assumed to equal
        # yield stress
        # G.Ito 8/16; JAO 9/15 for non-healed ep; fixed by JAO 4/2017

        dep_s=zeros(size(epsII_s));
        dep_s[s_nodes_yield] = dt_m.*max.(
            epsII_s[s_nodes_yield]-(
                yield_s[s_nodes_yield]-sqrt.(sxxOLD_s[s_nodes_yield].^2+sxyOLD[s_nodes_yield].^2)
                )./
            (2 .*Gs[s_nodes_yield].*dt_m), 
            minimum(epsII_s)*1e-6);

        depm=SiStER_interp_shear_nodes_to_markers(reshape(dep_s,Ny,Nx),x,y,xm,ym,icn,jcn);
        global ep = (ep+depm')./(dt_m/PARAMS.tau_heal+1);
        global epNH=epNH+depm';
        ### End of SiStER_update_ep ############################################
    end
  
    # OUTPUT VARIABLES OF INTEREST [prior to rotation & advection] --> not implemented
    # if (t%dt_out==0 && dt_out>0) || t==1 || t==Nt # SAVING SELECTED OUTPUT
        # println("SAVING SELECTED VARIABLES TO OUTPUT FILE") 
        # filename=string(t);
        # [etam]=SiStER_interp_shear_nodes_to_markers[etas,x,y,xm,ym,icn,jcn]; # to visualize viscosity on markers
        # save(filename,"X','Y','vx','vy','p','time','xm','ym','etam','rhom','BC','etan','Tm','im','idm','epsIIm','sxxm','sxym','ep','epNH','icn','jcn','qd','topo_x','topo_y")
    # end
    
    # SET ADVECTION TIME STEP BASED ON CURRENT FLOW SOLUTION
    # @time = 0.056206 seconds (112.89 k allocations: 7.596 MiB, 99.92% compilation time)
    global dt_m = SiStER_set_timestep(dx,dy,vx,vy,PARAMS);

    # ROTATE ELASTIC STRESSES IN CURRENT FLOW FIELD
    # @time = 0.229503 seconds (683.16 k allocations: 311.985 MiB, 5.28% gc time, 72.22% compilation time)
    if (PARAMS.YNElast==1) 
        ### Start of SiStER_rotate_stresses ####################################
        # update elastic stresses on markers following a solve [but before advection]
        ROT = SiStER_get_rotation_rate(vx,vy,dx,dy,BC) 
        om = SiStER_interp_shear_nodes_to_markers(ROT,x,y,xm,ym,icn,jcn);

        # rotate markers
        α = (om*dt_m)';
        sxymtemp = (sxxm.*sin.(2 .*α)) + (sxym.*cos.(2 .*α));
        sxxm = (sxxm.*(cos.(α).^2 - sin.(α).^2)) - (sxym.*sin.(2 .*α));
        global sxym=sxymtemp;
        ### End of SiStER_rotate_stresses ######################################
    end
    
    # EVOLVE TEMPERATURE FIELD THROUGH DIFFUSION
    # @time = 1.501444 seconds (25.62 M allocations: 1.542 GiB, 14.99% gc time, 8.72% compilation time)
    if PARAMS.Tsolve==1
        ### Start of SiStER_thermal_update #####################################
        # SiStER THERMAL SOLVE
        # get previous temperature on nodes
        n2interp = SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,Tm);
        Told = n2interp[1];

        # enforce Dirichlet boundary conditions to avoid mismatch between markers
        # & nodes
        if BCtherm.top[1]==1
            Told[1,:].=BCtherm.top[2];
        end
        if BCtherm.bot[1]==1
            Told[Ny,:].=BCtherm.bot[2];
        end
        if BCtherm.left[1]==1
            Told[:,1].=BCtherm.left[2];
        end
        if BCtherm.right[1]==1
            Told[:,Nx].=BCtherm.right[2];
        end

        # GET VARIABLE DIFFUSIVITY AND CP
        if hasproperty(MAT[1], :cp) || hasproperty(MAT[1], :k) 
            cpfield=PARAMS.cpref*ones(size(T));
            kfield=PARAMS.kref*ones(size(T));
            rhofield=PARAMS.rhoref*ones(size(T));    
        else
            # include("SiStER_get_thermal_properties.jl")
            # km, cpm = SiStER_get_thermal_properties(im, MAT)
            # include("SiStER_get_density.jl")
            # rhom = SiStER_get_density(im,Tm,MAT)
            # n2interp = SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,km,cpm,rhom);
            # kfield=n2interp[1];
            # cpfield=n2interp[2];
            # rhofield=n2interp[3];
            println(" !!! NOT IMPLEMENTED YET !!!")
        end

        # THERMAL SOLVE
        # @time = 0.086090 seconds (61.87 k allocations: 13.556 MiB, 94.96% compilation time)
        global T = SiStER_thermal_solver_sparse_CFD(x,y,Told,rhofield,cpfield,kfield,dt_m,BCtherm,zeros(size(T)))[1];

        # temperature change
        dT=T.-Told;
        # enforce Dirichlet boundary conditions to avoid mismatch between markers
        # & nodes
        if BCtherm.top[1]==1
            dT[1,:].=0;
        end
        if BCtherm.bot[1]==1
            dT[Ny,:].=0;
        end
        if BCtherm.left[1]==1
            dT[:,1].=0;
        end
        if BCtherm.right[1]==1
            dT[:,Nx].=0;
        end

        # @time = 0.106679 seconds (86 allocations: 226.101 MiB, 63.78% gc time)
        global Tm = SiStER_interp_shear_nodes_to_markers(T,x,y,xm,ym,icn,jcn);

        if PARAMS.ynTreset==1 # reset T=T0 in top layer
            Tm[im.==1].=PARAMS.T0;
        end
        ### End of SiStER_thermal_update #######################################
    end

    # MARKER ADVECTION; REMOVAL; AND ADDITION #############################
    ### Start of SiStER_move_remove_and_reseed_markers #########################
    # @time = 4.062938 seconds (4.23 M allocations: 3.528 GiB, 2.80% gc time, 78.01% compilation time)
    xm_new, ym_new, _, _ = SiStER_advect_markers(x,y,xm,ym,dx,dy,dt_m,vx,vy);

    global xm = copy(xm_new);
    global ym = copy(ym_new);

    # eliminate markers that left domain
    Iin=vec((xm.<=xsize) .& (xm.>=0) .& (ym.>=0) .& (ym.<=ysize));

    #msg2="  markers removed: ";
    #msg=[msg2 num2str(length(xm)-length(Iin))];
    #disp(msg)
    global xm=xm[Iin];
    global ym=ym[Iin];
    global im=im[Iin];
    global ep=ep[Iin];
    global epNH=epNH[Iin];
    global Tm=Tm[Iin];
    global idm=idm[Iin];
    global sxxm=sxxm[Iin];
    global sxym=sxym[Iin];
    global epsIIm=epsIIm[Iin];

    # locate advected markers with respect to the eulerian grid()
    # @time = 0.085912 seconds (176 allocations: 105.934 MiB)
    quadtemp, icntemp, jcntemp = SiStER_locate_markers_in_grid(xm,ym,x,y,dx,dy);
    global quad = quadtemp
    global icn = icntemp
    global jcn = jcntemp

    # check for holes in the marker distribution; 
    # patch with new markers if necessary
    # those new markers immediately get assigned a value of phase [im], index 
    # (idm) & accumulated plastic strain [ep], i.e., the 2 variables that never get()
    # passed to nodes. 
    # @time = 0.422957 seconds (1.91 M allocations: 121.505 MiB, 8.41% gc time, 99.63% compilation time)
    xm, ym, im, Ifix, mp, ep, idm, Tm, sxxm, sxym, epNH, epsIIm =SiStER_patch_marker_holes(icn,jcn,quad,Nx,Ny,Mquad,Mquad_crit,xm,ym,x,y,dx,dy,im,ep,idm,Tm,sxxm,sxym,epNH, epsIIm);    

    # then they get assigned P; epsII & stresses from grid values

    if minimum(Ifix)>0
        

        xmFIX=xm[Ifix];
        ymFIX=ym[Ifix];
        
        # pass temperature; pressure; strain rate & stresses to the new
        # markers from their nodal values
        # locate new markers with respect to the eulerian grid()
        quadFIX, icnFIX, jcnFIX = SiStER_locate_markers_in_grid(xmFIX,ymFIX,x,y,dx,dy); 

        temp = SiStER_interp_normal_nodes_to_markers(p,xc,yc,xmFIX,ymFIX,icnFIX,jcnFIX);
        # pm = zeros()
        # pm[Ifix]=temp; # pressure
        pm = temp';
        
        
    end
        
        
    # locate all markers with respect to the eulerian grid()
    qdtemp, icntemp, jcntemp = SiStER_locate_markers_in_grid(xm,ym,x,y,dx,dy);
    global qd = qdtemp
    global icn = icntemp
    global jcn = jcntemp
    ### End of SiStER_move_remove_and_reseed_markers ###########################

    # Do the same for the markers outlining topography
    ### Start of SiStER_update_topography_markers ##############################
    # advect the marker chain that keeps track of topography 
    # in the current flow field
    # @time = 3.550728 seconds (4.76 M allocations: 260.666 MiB, 99.96% compilation time)
    topo_xtemp, topo_ytemp = SiStER_advect_markers(x,y,topo_x,topo_y,dx,dy,dt_m,vx,vy);
    global topo_x = topo_xtemp
    global topo_y = topo_ytemp

    # locate the interface between sticky layer & left / right edge
    if findfirst(x -> x < 0, topo_x) === nothing
        topoL = topo_y[1];
    else
        # topoL=interp1(topo_x,topo_y,0);
        topoL = LinearInterpolation(vec(topo_x), vec(topo_y), extrapolation_bc=Line())(0);
        # topoL = Int.(topoL)
    end

    if findfirst(x -> x > xsize, topo_x) === nothing
        topoR=topo_y[end];
    else
        # topoR=interp1(topo_x,topo_y,xsize);
        topoR = LinearInterpolation(vec(topo_x), vec(topo_y), extrapolation_bc=Line())(xsize);
        # topoR = Int.(topoR)
    end

    # eliminate topography markers that left domain; keep the first one out on both sides
    Iin=findall(x -> 0<x<xsize, topo_x);
    topo_x=topo_x[Iin];
    topo_y=topo_y[Iin];
    topo_x=push!(pushfirst!(topo_x, 0), xsize);
    topo_y=push!(pushfirst!(topo_y, topoL), topoR);

    # @time = 0.314203 seconds (751.34 k allocations: 52.225 MiB, 93.01% compilation time: 2% of which was recompilation)
    if PARAMS.YNSurfaceProcesses==1
        # ERODE TOPOGRAPHY
        topo_y=SiStER_topography_diffusion_solver(topo_x,topo_y,dt_m,PARAMS.topo_kappa);
        # RESET ROCK AND AIR [assumes topography is only interface between phase 1 & 2]
        # topomarkers=interp1(topo_x,topo_y,xm);
        topomarkers = linear_interpolation(vec(topo_x), vec(topo_y), extrapolation_bc=Line())(xm)
        im[(im.==1) .& (ym.>=topomarkers)] .= 2;
        im[(im.>=2) .& (ym.<topomarkers)] .= 1;
    end

    # if there has been too much stretching; regrid the surface topography
    # @time = 0.006591 seconds (422 allocations: 26.875 KiB, 99.53% compilation time)
    if maximum(diff(topo_x))>5*topo_marker_spacing || !issorted(topo_x)
        # surface regridding happens if somewhere 2 topo markers have been
        # stretched apart by more than 5 times the inital mean marker spacing
        # | if topo_x is no longer sorted due to compression.
        topo_xREGRID=range(0,xsize,length=Ntopo);
        # topo_yREGRID=interp1(topo_x,topo_y,topo_xREGRID[2:end-1]);
        topo_yREGRID = linear_interpolation(vec(topo_x), vec(topo_y), extrapolation_bc=Line())(topo_xREGRID[2:end-1])
        # topo_yREGRID=[topoL topo_yREGRID topoR];
        topo_y = push!(pushfirst!(topo_yREGRID, topoL), topoR);
        global topo_x=topo_xREGRID;
        global topo_y=topo_yREGRID;
        println("**REGRIDDING TOPOGRAPHY MARKERS**")
    end

    ### End of SiStER_update_topography_markers ################################
 
 
    #######################################################################
    
    # Create topography plot --> slow first time it is plotted
    plot(topo_x, topo_y)
    # scatter!(topo_x, topo_y, c="red")
    # plot!(xlims=(0,Inf), ylims=(0, Inf))
    plot!(size=(1000,180))
    yflip!(true)
    plot!(title = string(round(time/365.25/24/3600/1000, digits=3))*" kyrs.", xlabel = "x", ylabel = "y")
    savefig(pwd()*"/figures/"*string(t)*"-topography-plot.png")  

    # Create marker plot
    # fastscatter(xm(im>1)/1e3,ym(im>1) /1e3,log10(epsIIm(im>1)),'markersize’,2);
    # plotmask = im .> 1
    # plotmask = sample(axes(xm, 1), 10000)
    # scatter(xm[plotmask]./1000, ym[plotmask]./1000, marker_z=log10.(epsIIm[plotmask]), color=:heat, markershape=:circle, msc=:transparent, markersize=5)
    # yflip!(true)
    end

    println("TIME ELAPSED: "*string(timing)*" seconds")

    println("---------------")
    println("END OF ITERATION: "* string(t) *" out of "* string(Nt) *" - SIMULATION TIME: "* string(time/365.25/24/3600/1000) *" kyrs.")
    println("--------------------------------")
    println("--------------------------------")
end

print("FIN")
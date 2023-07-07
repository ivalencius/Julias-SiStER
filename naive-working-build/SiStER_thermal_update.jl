# SiStER THERMAL SOLVE

# get previous temperature on nodes
include("SiStER_interp_markers_to_shear_nodes.jl")
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
include("SiStER_thermal_solver_sparse_CFD.jl")
T = SiStER_thermal_solver_sparse_CFD(x,y,Told,rhofield,cpfield,kfield,dt_m,BCtherm,zeros(size(T)));


# temperature change
dT=T-Told;
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

Tm = SiStER_interp_shear_nodes_to_markers(T,x,y,xm,ym,icn,jcn);

if PARAMS.ynTreset==1 # reset T=T0 in top layer
     Tm[im.==1].=PARAMS.T0;
end
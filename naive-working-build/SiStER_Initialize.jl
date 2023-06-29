# SiStER Initialize
using Statistics

PARAMS.Nphase = Nphase; # for convenience

# construct staggered grids
include("SiStER_initialize_grid.jl")
(X,Y,x,y,xc,yc,dx,dy,Nx,Ny) = SiStER_initialize_grid(xsize,ysize,GRID);

# initialize marker arrays & positions
include("SiStER_initialize_marker_positions.jl")
(xm, ym) = SiStER_initialize_marker_positions(xsize,ysize,dx,dy,Mquad); # THIS IS VERY SLOW
print("**Number of markers: "*string(length(xm))*"**\n");

# locate markers with respect to grid()
include("SiStER_locate_markers_in_grid.jl")
(qd,icn,jcn) = SiStER_locate_markers_in_grid(xm,ym,x,y,dx,dy);
icn = Int.(icn);
jcn = Int.(jcn);

# assign marker phases
include("SiStER_initialize_marker_phases.jl")
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
include("SiStER_interp_shear_nodes_to_markers.jl")
Tm = SiStER_interp_shear_nodes_to_markers(T,x,y,xm,ym,icn,jcn);
Tm = Tm;
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
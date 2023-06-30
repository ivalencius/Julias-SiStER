# SiStER_update_topography_markers

# advect the marker chain that keeps track of topography 
# in the current flow field
include("SiStER_advect_markers.jl")
topo_x, topo_y = SiStER_advect_markers(x,y,topo_x,topo_y,dx,dy,dt_m,vx,vy);

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

if PARAMS.YNSurfaceProcesses==1
    # ERODE TOPOGRAPHY
    include("SiStER_topography_diffusion_solver.jl")
    topo_y=SiStER_topography_diffusion_solver(topo_x,topo_y,dt_m,PARAMS.topo_kappa);
    # RESET ROCK AND AIR [assumes topography is only interface between phase 1 & 2]
    # topomarkers=interp1(topo_x,topo_y,xm);
    topomarkers = linear_interpolation(vec(topo_x), vec(topo_y), extrapolation_bc=Line())(xm)
    im[(im.==1) .& (ym.>=topomarkers)] .= 2;
    im[(im.>=2) .& (ym.<topomarkers)] .= 1;
end

# if there has been too much stretching; regrid the surface topography
if maximum(diff(topo_x))>5*topo_marker_spacing || !issorted(topo_x)
    # surface regridding happens if somewhere 2 topo markers have been
    # stretched apart by more than 5 times the inital mean marker spacing
    # | if topo_x is no longer sorted due to compression.
    topo_xREGRID=range(0,xsize,length=Ntopo);
    # topo_yREGRID=interp1(topo_x,topo_y,topo_xREGRID[2:end-1]);
    topo_yREGRID = linear_interpolation(vec(topo_x), vec(topo_y), extrapolation_bc=Line())(topo_xREGRID[2:end-1])
    # topo_yREGRID=[topoL topo_yREGRID topoR];
    topo_y = push!(pushfirst!(topo_REGRIDy, topoL), topoR);
    topo_x=topo_xREGRID;
    topo_y=topo_yREGRID;
    println("**REGRIDDING TOPOGRAPHY MARKERS**")
end

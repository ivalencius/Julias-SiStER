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

# INITIALIZATION

# Input File: loads parameter values; model geometry; boundary conditions
# if exist("running_from_SiStER_RUN','var")==0
#     clear() 
#     InpFil = input("Input file ? ','s");
# end
# run(InpFil)
InpFil = "SiStER_Input_File_continental_rift.jl"
include(InpFil)

# Make folder for figures
mkpath("figures/")

# Redefined functions --> use as needed
# function meshgrid(x, y)
#     X = [i for i in x, _ in 1:length(y)]
#     Y = [j for _ in 1:length(x), j in y]
#     return X', Y'
# end

# Sourced from https://github.com/ChrisRackauckas/VectorizedRoutines.jl/blob/master/src/matlab.jl
# function accumarray2(subs, val, fun=sum, fillval=0; sz=maximum(subs,1), issparse=false)
#     counts = Dict()
#     for i = 1:size(subs,1)
#          counts[subs[i,:]]=[get(counts,subs[i,:],[]);val[i...]]
#     end
#     A = fillval*ones(sz...)
#     for j = keys(counts)
#          A[j...] = fun(counts[j])
#     end
#     issparse ? sparse(A) : A
#  end

# interp(x, v, xq) --> interpolate((x,), v, Gridded(Constant()))(xq)

# construct grid & initialize marker / node arrays
println("---------------")
println("INITIALIZING")
println("---------------")
include("SiStER_Initialize.jl");

println("For continental drift: Nx = 106, NY=43 (from MATLAB solution)")
println("Mismatch due to round() rounding 4.5 to nearest even integer, not up to 5 like in MATLAB")
println("SiStER_initialize_grid.jl (line 40)")

# BEGIN TIME LOOP #########################################################
global time=0; # Need to rename time variable??

# for t=1:Nt # time loop
using Plots
pyplot()

global t = 1

# Start solve
for temp = 1:Nt
# for temp = 1:10
    global t = temp
    println("STARTING ITERATION: " * string(t) * " out of " * string(Nt))
    
    # update time
    global time=time+dt_m;
    
    # Here we prepare nodal arrays to feed the Stokes solver 
    include("SiStER_material_props_on_nodes.jl");

    ### SOLVE STOKES WITH NON-LINEAR RHEOLOGY HERE 
    include("SiStER_flow_solve.jl");

    # using JLD2
    # # @save "running-vars.jld2"
    # @load "running-vars.jld2"
    
    # GET STRAIN RATE FROM CURRENT SOLUTION
    include("SiStER_interp_shear_nodes_to_markers.jl")
    global epsIIm = SiStER_interp_shear_nodes_to_markers(epsII_s,x,y,xm,ym,icn,jcn);
    
    # USE STRAIN RATE TO UPDATE STRESSES ON MARKERS
    include("SiStER_update_marker_stresses.jl");
    
    # BUILD UP PLASTIC STRAIN IN YIELDING AREAS IF PLASTICITY IS ACTIVATED
    if (PARAMS.YNPlas==1)
        include("SiStER_update_ep.jl");
    end
  
#     # OUTPUT VARIABLES OF INTEREST [prior to rotation & advection]
#     if (t%dt_out==0 && dt_out>0) || t==1 || t==Nt # SAVING SELECTED OUTPUT
#         disp("SAVING SELECTED VARIABLES TO OUTPUT FILE") 
#         filename=string(t);
#         [etam]=SiStER_interp_shear_nodes_to_markers[etas,x,y,xm,ym,icn,jcn]; # to visualize viscosity on markers
#         save(filename,"X','Y','vx','vy','p','time','xm','ym','etam','rhom','BC','etan','Tm','im','idm','epsIIm','sxxm','sxym','ep','epNH','icn','jcn','qd','topo_x','topo_y")
#     end
    
    # SET ADVECTION TIME STEP BASED ON CURRENT FLOW SOLUTION
    include("SiStER_set_timestep.jl")
    global dt_m = SiStER_set_timestep(dx,dy,vx,vy,PARAMS);

    # ROTATE ELASTIC STRESSES IN CURRENT FLOW FIELD
    if (PARAMS.YNElast==1) 
        include("SiStER_rotate_stresses.jl");
    end
    
    # EVOLVE TEMPERATURE FIELD THROUGH DIFFUSION
    if PARAMS.Tsolve==1
        include("SiStER_thermal_update.jl");
    end

    # MARKER ADVECTION; REMOVAL; AND ADDITION #############################
    include("SiStER_move_remove_and_reseed_markers.jl");
    # advect markers in current flow field
    # remove markers if necessary
    # add markers if necessary
    include("SiStER_update_topography_markers.jl");
    # here we do the same for the marker chain that keeps track of topography
    #######################################################################

    println("---------------")
    println("END OF ITERATION: "* string(t) *" out of "* string(Nt) *" - SIMULATION TIME: "* string(time/365.25/24/3600/1000) *" kyrs.")
    println("--------------------------------")
    println("--------------------------------")
    
    # Create topography plot
    plot(topo_x, topo_y)
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

print("FIN")

    
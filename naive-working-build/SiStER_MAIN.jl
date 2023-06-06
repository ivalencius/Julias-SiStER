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

# Redefined functions --> use as needed
# function meshgrid(x, y)
#     X = [i for i in x, _ in 1:length(y)]
#     Y = [j for _ in 1:length(x), j in y]
#     return X', Y'
# end

# construct grid & initialize marker / node arrays
print("---------------\n")
print("INITIALIZING\n")
print("---------------\n")
print("For continental drift: Nx = 106, NY=43 (from MATLAB solution)\n")
print("Mismatch due to round() rounding 4.5 to nearest even integer, not up to 5 like in MATLAB\n")
print("SiStER_initialize_grid.jl (line 40)\n")
include("SiStER_Initialize.jl")

using JLD2
@save "running-vars.jld2"
# @load "running-vars.jld2"

# BEGIN TIME LOOP #########################################################
time=0;

# for t=1:Nt # time loop
    
#     print("STARTING ITERATION: " * string(t) * " out of " * string(Nt)*"\n")
    
#     # update time
#     time=time+dt_m;
    
#     # Here we prepare nodal arrays to feed the Stokes solver 
#     SiStER_material_props_on_nodes

#     ### SOLVE STOKES WITH NON-LINEAR RHEOLOGY HERE 
#     SiStER_flow_solve
    
#     # GET STRAIN RATE FROM CURRENT SOLUTION
#     epsIIm=SiStER_interp_shear_nodes_to_markers[epsII_s,x,y,xm,ym,icn,jcn];
    
#     # USE STRAIN RATE TO UPDATE STRESSES ON MARKERS
#     SiStER_update_marker_stresses;
    
#     # BUILD UP PLASTIC STRAIN IN YIELDING AREAS IF PLASTICITY IS ACTIVATED
#     if (PARAMS.YNPlas==1) 
#         SiStER_update_ep;
#     end
  
#     # OUTPUT VARIABLES OF INTEREST [prior to rotation & advection]
#     if (t%dt_out==0 && dt_out>0) || t==1 || t==Nt # SAVING SELECTED OUTPUT
#         disp("SAVING SELECTED VARIABLES TO OUTPUT FILE") 
#         filename=string(t);
#         [etam]=SiStER_interp_shear_nodes_to_markers[etas,x,y,xm,ym,icn,jcn]; # to visualize viscosity on markers
#         save(filename,"X','Y','vx','vy','p','time','xm','ym','etam','rhom','BC','etan','Tm','im','idm','epsIIm','sxxm','sxym','ep','epNH','icn','jcn','qd','topo_x','topo_y")
#     end
    
#     # SET ADVECTION TIME STEP BASED ON CURRENT FLOW SOLUTION
#     [dt_m]=SiStER_set_timestep[dx,dy,vx,vy,PARAMS];

#     # ROTATE ELASTIC STRESSES IN CURRENT FLOW FIELD
#     if (PARAMS.YNElast==1) 
#         SiStER_rotate_stresses;
#     end
    
#     # EVOLVE TEMPERATURE FIELD THROUGH DIFFUSION
#     if PARAMS.Tsolve==1
#         SiStER_thermal_update;
#     end

#     # MARKER ADVECTION; REMOVAL; AND ADDITION #############################
#     SiStER_move_remove_and_reseed_markers;
#     # advect markers in current flow field
#     # remove markers if necessary
#     # add markers if necessary
#     SiStER_update_topography_markers
#     # here we do the same for the marker chain that keeps track of topography
#     #######################################################################

#     print("---------------")
#     print(["END OF ITERATION: ' num2str(t) ' out of ' num2str(Nt) ' - SIMULATION TIME: ' num2str(time/365.25/24/3600/1000) ' kyrs."])
#     print("--------------------------------")
#     print("--------------------------------")
    

# end

# print("FIN")

    
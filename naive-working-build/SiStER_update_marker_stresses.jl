# Updates stresses on markers for CURRENT solution.  Stress rotation
# occurs after solutions are output
# G.Ito 8/16

# Compute STRESS Changes on nodes; interpolate to markers; & apply to marker stresses
dsxx=(2*etan.*EXX-sxxOLD).*Zn;
dsxy=(2*etas.*EXY-sxyOLD).*Zs;

include("SiStER_interp_normal_nodes_to_markers.jl")
dsxxm=SiStER_interp_normal_nodes_to_markers(dsxx,xc,yc,xm,ym,icn,jcn);
include("SiStER_interp_shear_nodes_to_markers.jl")
dsxym=SiStER_interp_shear_nodes_to_markers(dsxy,x,y,xm,ym,icn,jcn);
sxxm=sxxm+dsxxm';
sxym=sxym+dsxym';
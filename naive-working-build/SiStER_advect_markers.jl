include("SiStER_locate_markers_in_grid.jl")
include("SiStER_get_marker_velocities.jl")

function SiStER_advect_markers(x,y,xm,ym,dx,dy,tstep,vx,vy)
# [xm_new,ym_new,vx_eff,vy_eff] = SiStER_advect_markers(x,y,xm,ym,dx,dy,tstep,vx,vy)
#
# Uses a Fourth-Order Runge-Kutta method to advect markers 
# in the current velocity field; over the advection time step tstep = dt_m
#
# First cut J.-A. Olive; 2011-2012 - modified by B.Z. Klein 2013-2014


# point A [particle pts]
XA = copy(xm);
YA = copy(ym);                          
qd, icn, jcn = SiStER_locate_markers_in_grid(XA,YA,x,y,dx,dy);
VxA, VyA = SiStER_get_marker_velocities(qd,icn,jcn,x,y,xm,ym,vx,vy,dx,dy);

# point B
XB = vec(XA .+ ((0.5*tstep).*VxA)');
YB = vec(YA .+ ((0.5*tstep).*VyA)');
qd,icn,jcn = SiStER_locate_markers_in_grid(XB,YB,x,y,dx,dy); 
VxB, VyB = SiStER_get_marker_velocities(qd,icn,jcn,x,y,xm,ym,vx,vy,dx,dy);

# point C
XC = vec(XA .+ ((0.5*tstep).*VxB)');
YC = vec(YA + ((0.5*tstep).*VyB)');
qd,icn,jcn = SiStER_locate_markers_in_grid(XC,YC,x,y,dx,dy);
VxC,VyC = SiStER_get_marker_velocities(qd,icn,jcn,x,y,xm,ym,vx,vy,dx,dy);

# point D
XD = vec(XA .+ (tstep.*VxC)');
YD = vec(YA .+ (tstep.*VyC)');
qd,icn,jcn = SiStER_locate_markers_in_grid(XD,YD,x,y,dx,dy); 
VxD,VyD = SiStER_get_marker_velocities(qd,icn,jcn,x,y,xm,ym,vx,vy,dx,dy);


# effective velocity
vx_eff = (1/6).*(VxA + (2 .* VxB) + (2 .* VxC) + VxD);
vy_eff = (1/6).*(VyA + (2 .* VyB) + (2 .* VyC) + VyD);


## Calculate new coordinates

PP1 = XA + (tstep.*vx_eff)';
PP2 = YA + (tstep*vy_eff)';

xm_new = PP1';
ym_new = PP2';

return xm_new, ym_new, vx_eff, vy_eff
end
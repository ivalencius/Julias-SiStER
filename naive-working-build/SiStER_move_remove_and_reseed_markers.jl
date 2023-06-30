# SiStER Update Markers
include("SiStER_advect_markers.jl")
xm_new, ym_new, _, _ = SiStER_advect_markers(x,y,xm,ym,dx,dy,dt_m,vx,vy);

xm = copy(xm_new);
ym = copy(ym_new);

# eliminate markers that left domain
Iin=vec((xm.<=xsize) .& (xm.>=0) .& (ym.>=0) .& (ym.<=ysize));

#msg2="  markers removed: ";
#msg=[msg2 num2str(length(xm)-length(Iin))];
#disp(msg)
xm=xm[Iin];
ym=ym[Iin];
im=im[Iin];
ep=ep[Iin];
epNH=epNH[Iin];
Tm=Tm[Iin];
idm=idm[Iin];
sxxm=sxxm[Iin];
sxym=sxym[Iin];
epsIIm=epsIIm[Iin];

# locate advected markers with respect to the eulerian grid()
include("SiStER_locate_markers_in_grid.jl")
quad,icn,jcn = SiStER_locate_markers_in_grid(xm,ym,x,y,dx,dy);
    
# check for holes in the marker distribution; 
# patch with new markers if necessary
# those new markers immediately get assigned a value of phase [im], index 
# (idm) & accumulated plastic strain [ep], i.e., the 2 variables that never get()
# passed to nodes. 

include("SiStER_patch_marker_holes.jl")
xm, ym, im, Ifix, mp, ep, idm, Tm, sxxm, sxym, epNH, epsIIm =SiStER_patch_marker_holes(icn,jcn,quad,Nx,Ny,Mquad,Mquad_crit,xm,ym,x,y,dx,dy,im,ep,idm,Tm,sxxm,sxym,epNH, epsIIm);    

# then they get assigned P; epsII & stresses from grid values

if minimum(Ifix)>0
    

    xmFIX=xm[Ifix];
    ymFIX=ym[Ifix];
    
    # pass temperature; pressure; strain rate & stresses to the new
    # markers from their nodal values
    # locate new markers with respect to the eulerian grid()
    quadFIX, icnFIX, jcnFIX = SiStER_locate_markers_in_grid(xmFIX,ymFIX,x,y,dx,dy); 
    
    include("SiStER_interp_normal_nodes_to_markers.jl")
    temp = SiStER_interp_normal_nodes_to_markers(p,xc,yc,xmFIX,ymFIX,icnFIX,jcnFIX);
    # pm = zeros()
    # pm[Ifix]=temp; # pressure
    pm = temp';
    
    
end
    
    
# locate all markers with respect to the eulerian grid()
qd,icn,jcn = SiStER_locate_markers_in_grid(xm,ym,x,y,dx,dy);
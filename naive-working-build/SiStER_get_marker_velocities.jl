using GroupSlices
include("SiStER_interp_grid_to_marker_vector.jl")

function SiStER_get_marker_velocities(quad,icn,jcn,x,y,xm,ym,vx,vy,dx,dy)
#
# [VX,VY] = SiStER_get_marker_velocities(quad,icn,jcn,x,y,xm,ym,vx,vy,dx,dy)
# interpolates velocities at nodes to velocities on markers
#
# First cut - J.A. Olive; 2011-2012; modified by B.Z. Klein 2013-2014

Nx = length(x);
Ny = length(y);
(m, _) = size(vx);

#boolean arrays for interpolation
quad1 = quad .== 1;
quad2 = quad .== 2;
quad3 = quad .== 3;
quad4 = quad .== 4;

top  = icn .== 1;
bot  = icn .== Ny-1;

left = jcn .== 1;
right= jcn .== Nx-1;

#linear index vector of icn & jcn for vx; vy
# IND = sub2ind(size(vx), icn, jcn);
s2i = LinearIndices(size(vx)); # reverse is CartesianIndices
IND = zeros(length(icn));
for ind ∈ 1:length(icn)
    IND[ind] = s2i[icn[ind], jcn[ind]]
end
IND = Int.(IND);

######################
## Vx Interpolation ##
######################

Vx_xnodes = [x[jcn]' x[jcn.+1]' x[jcn.+1]' x[jcn]']'; # 4 x num Mark

# Create a bunch of new arrays
Vx_ynodes = zeros(size(Vx_xnodes)) # Do I need to create a new array?
Vx_nodes = zeros(size(Vx_xnodes)) # Do I need to create a new array?
Vy_xnodes = zeros(size(Vx_xnodes)) # Do I need to create a new array?
Vy_nodes = zeros(size(Vx_xnodes)) # Do I need to create a new array?

## case 1, icn .== 1 [top]
c1 = vec(top);

Vx_ynodes[:,c1] =  cat(
    y[icn[c1]]+dy[icn[c1]]./2, 
    y[icn[c1]]+dy[icn[c1]]./2,
    y[icn[c1].+1]+dy[icn[c1].+1]./2,      
    y[icn[c1].+1]+dy[icn[c1].+1]./2, 
    dims = (2,2))';


#[y[icn[c1]]+dy[icn[c1]]./2 y[icn[c1]]+dy[icn[c1]]./2     y[icn[c1].+1]+dy[icn[c1].+1]./2      y[icn[c1].+1]+dy[icn[c1].+1]./2]';
Vx_nodes[:,c1]  =  cat(
    vx[IND[c1]],
    vx[IND[c1].+m],
    vx[IND[c1].+(1+m)],
    vx[IND[c1].+1],
    dims = (2,2))';

## case 2, icn .== Ny-1 [bot] && in quad 3 | 4
c2 = vec(bot); #& (quad3 | quad4);

Vx_ynodes[:,c2] = cat(
    y[icn[c2].-1]+dy[icn[c2].-1]./2,
    y[icn[c2].-1]+dy[icn[c2].-1]./2,
    y[icn[c2]]+dy[icn[c2]]./2,
    y[icn[c2]]+dy[icn[c2]]./2,
    dims = (2,2))';

Vx_nodes[:,c2] = cat(
    vx[IND[c2].-1],
    vx[IND[c2].+(m-1)],
    vx[IND[c2].+m],
    vx[IND[c2]],
    dims = (2,2))';

## case 3 icn ~top & ~bot; quad = 1 | 2
c3 = vec(.!top .& .!bot .& (quad1 .| quad2));

Vx_ynodes[:,c3] = cat(
    y[icn[c3]]-dy[icn[c3].-1]./2,
    y[icn[c3]]-dy[icn[c3].-1]./2,
    y[icn[c3]]+dy[icn[c3]]./2,
    y[icn[c3]]+dy[icn[c3]]./2,
    dims = (2,2))';
             
Vx_nodes[:,c3] = cat(
    vx[IND[c3].-1],
    vx[IND[c3].+(m-1)],
    vx[IND[c3].+m],
    vx[IND[c3]],
    dims = (2,2))';
            
## case 4 icn ~top & ~bot & quad = 3|4
c4 = vec(.!top .& .!bot .& (quad3 .| quad4));

Vx_ynodes[:,c4] = cat(
    y[icn[c4]]+dy[icn[c4]]./2,
    y[icn[c4]]+dy[icn[c4]]./2,
    y[icn[c4].+1]+dy[icn[c4].+1]./2,
    y[icn[c4].+1]+dy[icn[c4].+1]./2,
    dims = (2,2))';

Vx_nodes[:,c4] = cat(
    vx[IND[c4]],
    vx[IND[c4].+m],
    vx[IND[c4].+(1+m)],
    vx[IND[c4].+1],
    dims = (2,2))';

            
            
VX = SiStER_interp_grid_to_marker_vector(Vx_xnodes, Vx_ynodes, Vx_nodes, xm, ym);



#######################
## VY INTERPOLATION  ##
#######################

Vy_ynodes = [y[icn]' y[icn]' y[icn.+1]' y[icn.+1]']';

## case 1 jcn .== 1 [left]
c1 = vec(left);

Vy_xnodes[:,c1] = cat(
    x[jcn[c1]]+dx[jcn[c1]]./2,
    x[jcn[c1].+1]+dx[jcn[c1].+1]./2,
    x[jcn[c1].+1]+dx[jcn[c1].+1]./2,
    x[jcn[c1]]+dx[jcn[c1]]./2,
    dims = (2,2))';

Vy_nodes[:,c1] = cat(
    vy[IND[c1]],
    vy[IND[c1].+m],
    vy[IND[c1].+(1+m)],
    vy[IND[c1].+1],
    dims = (2,2))';

## case 2 jcn .== Nx-1 [right] & Quad = 2 | 3
c2 = vec(right);# & (quad2 | quad3);???????????????????????????

Vy_xnodes[:,c2] = cat(
    x[jcn[c2].-1]+dx[jcn[c2].-1]./2,
    x[jcn[c2]]+dx[jcn[c2]]./2,
    x[jcn[c2]]+dx[jcn[c2]]./2,
    x[jcn[c2].-1]+dx[jcn[c2].-1]./2,
    dims = (2,2))';
           
Vy_nodes[:,c2] = cat(
    vy[IND[c2].-m],
    vy[IND[c2]],
    vy[IND[c2].+1],
    vy[IND[c2].+(1-m)],
    dims = (2,2))';

 ## case 3 ~left & ~right & quad = 1 | 4
 c3 = vec(.!left .& .!right .& (quad1 .| quad4));
 
 Vy_xnodes[:,c3] = cat(
    x[jcn[c3]]-dx[jcn[c3].-1]./2,
    x[jcn[c3]]+dx[jcn[c3]]./2,
    x[jcn[c3]]+dx[jcn[c3]]./2,
    x[jcn[c3]]-dx[jcn[c3].-1]./2,
    dims = (2,2))';
              
 Vy_nodes[:,c3] = cat(
    vy[IND[c3].-m],
    vy[IND[c3]],
    vy[IND[c3].+1],
    vy[IND[c3].+(1-m)],
    dims = (2,2))';

## case 4 ~left & ~right & quad = 2 | 3
c4 = vec(.!left .& .!right .& (quad2 .| quad3));

Vy_xnodes[:,c4] = cat(
    x[jcn[c4]]+dx[jcn[c4]]./2,
    x[jcn[c4].+1]+dx[jcn[c4].+1]./2,
    x[jcn[c4].+1]+dx[jcn[c4].+1]./2,
    x[jcn[c4]]+dx[jcn[c4]]./2,
    dims = (2,2))';

Vy_nodes[:,c4] = cat(
    vy[IND[c4]],
    vy[IND[c4].+m],
    vy[IND[c4].+(m+1)],
    vy[IND[c4].+1],
    dims = (2,2))';
 

            
VY = SiStER_interp_grid_to_marker_vector(Vy_xnodes,Vy_ynodes,Vy_nodes, xm, ym);    

return VX, VY

end
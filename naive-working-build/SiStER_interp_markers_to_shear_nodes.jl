using Interpolations
using RecursiveArrayTools

function SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,quad,x,y,varargin)
# [n2interp] = SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,quad,x,y,varargin)
# interpolates marker properties to shear nodes
# First cut - J.A. Olive; March 2011
# Modified by E. Mittelstaedt; April 2011; to allow variable inputs.  
# Modified by B.Z. Klein; Spring 2014; for speedup
# Modified by B.Z. Klein, Summer 2014, for further speedup [vectorized]

function accumarray(subs, val, fun=sum, fillval=0; sz=maximum(subs,1), issparse=false)
    counts = Dict()
    for i = 1:size(subs,1)
         counts[subs[i,:]]=[get(counts,subs[i,:],[]);val[i...]]
    end
    A = fillval*ones(sz...)
    for j = keys(counts)
         A[j...] = fun(counts[j])
    end
    issparse ? sparse(A) : A
end

Nx=length(x);
Ny=length(y);
dx=diff(x);
dy=diff(y);

# MITTELSTAEDT - check for number of properties to interpolate
numV = size(varargin, 2);

# MITTELSTAEDT # establish interpolants matrices
# Converted from array of (1, numV) filled with (Ny, Nx) structs
# n2interp = repeat(zeros(Ny,Nx), outer = (1, 1, numV));
# Read: https://stackoverflow.com/questions/47818035/julia-three-dimensional-arrays-performance
# test n2interp = zeros(Ny, Nx, numV);
n2interp = VectorOfArray([zeros(Ny, Nx) for _ in 1:numV]);

# JCN = interp1(x, 1:length(x), xm, "nearest', 'extrap"); ## these are like the jcn & icn elsewhere, except the nodes are centered instead of upper left.
# ICN = interp1(y, 1:length(y), ym, "nearest', 'extrap"); ## this makes a lot of the indexing much simpler below.
JCN = interpolate((x,), 1:length(x), Gridded(Constant()))(xm);
ICN = interpolate((y,), 1:length(y), Gridded(Constant()))(ym);

## Interior Cells

center = (jcn.>1) .& (jcn.<Nx) .& (icn.>1) .& (icn.<Ny);
shiftLeft = (jcn.<Nx-1) .& (icn.>1) .& (icn.<Ny);
shiftUp = (jcn.>1) .& (jcn.<Nx) .& (icn.<Ny-1);
shiftBoth = (jcn.<Nx-1) .& (icn.<Ny-1);

xJCN = x[JCN];
yICN = y[ICN];

# Out of memory error -->  do in loop
# cell1 = center .& ((xm.-xJCN) .> 0) .& ((ym .- yICN) .> 0);  ## these are logical arrays that index the original quadrants
# cell2 = shiftLeft .& ((xm.-xJCN) .< 0) .& ((ym .- yICN) .> 0);
# cell3 = shiftBoth .& ((xm.-xJCN) .< 0) .& ((ym .- yICN) .< 0);
# cell4 = shiftUp .& ((xm.-xJCN) .> 0) .& ((ym .- yICN) .< 0);
cell1 = BitVector(undef, length(xm));
cell2 = BitVector(undef, length(xm));
cell3 = BitVector(undef, length(xm));
cell4 = BitVector(undef, length(xm));

for i ∈ 1:length(cell1)
    cell1[i] = center[i] & ((xm[i]-xJCN[i]) > 0) & ((ym[i] - yICN[i]) > 0);  ## these are logical arrays that index the original quadrants
    cell2[i] = shiftLeft[i] & ((xm[i]-xJCN[i]) < 0) & ((ym[i] - yICN[i]) > 0);
    cell3[i] = shiftBoth[i] & ((xm[i]-xJCN[i]) < 0) & ((ym[i] - yICN[i]) < 0);
    cell4[i] = shiftUp[i] & ((xm[i]-xJCN[i]) > 0) & ((ym[i] - yICN[i]) < 0);
end

### WEIGHTING [equal for now because that is what I'm running]

wc1 = 0.25;
wc2 = 0.25;
wc3 = 0.25;
wc4 = 0.25;

# cell 1 [i,j,1]

dxm = xm[cell1] .- x[JCN[cell1]];
dym = ym[cell1] .- y[ICN[cell1]];
ddx = dx[JCN[cell1]];
ddy = dy[ICN[cell1]];

wm1 = 1 .- (dxm.*dym + (ddx.-dxm).*dym + (ddy.-dym).*dxm)./(ddx.*ddy);
w1 = accumarray(hcat(ICN[cell1], JCN[cell1]), wm1, sz=(Ny, Nx));

# cell 2 [i, j-1, 2]

dxm = xm[cell2] .- x[JCN[cell2]];
dym = ym[cell2] .- y[ICN[cell2]];
ddx = dx[JCN[cell2].-1];
ddy = dy[ICN[cell2]];

wm2 = 1 .- (dxm.*dym + (ddx.-dxm).*dym + (ddy.-dym).*dxm)./(ddx.*ddy);
w2 = accumarray(hcat(ICN[cell2], JCN[cell2]), wm2, sz=(Ny, Nx));

# cell 3 [i-1, j-1, 3]

dxm = xm[cell3] .- x[JCN[cell3]];
dym = ym[cell3] .- y[ICN[cell3]];
ddx = dx[JCN[cell3].-1];
ddy = dy[ICN[cell3].-1];

wm3 = 1 .- (dxm.*dym + (ddx.-dxm).*dym + (ddy.-dym).*dxm)./(ddx.*ddy);
w3 = accumarray(hcat(ICN[cell3], JCN[cell3]), wm3, sz=(Ny, Nx));

# cell 4 [i-1, j, 4]

dxm = xm[cell4] .- x[JCN[cell4]];
dym = ym[cell4] .- y[ICN[cell4]];
ddx = dx[JCN[cell4]];
ddy = dy[ICN[cell4].-1];

wm4 = 1 .- (dxm.*dym + (ddx.-dxm).*dym + (ddy.-dym).*dxm)./(ddx.*ddy);
w4 = accumarray(hcat(ICN[cell4], JCN[cell4]), wm4, sz=(Ny, Nx));

#loop over material properties to interpolate
# THIS IS VERY SLOW --> using varargin[vn]*cell1[vn]
# vn = 50000 it yields 0.0025 for all
for vn = 1:numV
    # n2interp[vn] = (
    #     wc1.*accumarray(hcat(ICN[cell1], JCN[cell1]), (varargin[vn]*cell1[vn]).*wm1, sz=(Ny, Nx))./w1 +
    #     wc2.*accumarray(hcat(ICN[cell2], JCN[cell2]), (varargin[vn]*cell2[vn]).*wm2, sz=(Ny, Nx))./w2 +
    #     wc3.*accumarray(hcat(ICN[cell3], JCN[cell3]), (varargin[vn]*cell3[vn]).*wm3, sz=(Ny, Nx))./w3 +
    #     wc4.*accumarray(hcat(ICN[cell4], JCN[cell4]), (varargin[vn]*cell4[vn]).*wm4, sz=(Ny, Nx))./w4
    #     )./(wc1+wc2+wc4+wc4);
    n2interp[vn] = (
        wc1.*accumarray(hcat(ICN[cell1], JCN[cell1]), varargin[:, vn][cell1].*wm1, sz=(Ny, Nx))./w1 +
        wc2.*accumarray(hcat(ICN[cell2], JCN[cell2]), varargin[:, vn][cell2].*wm2, sz=(Ny, Nx))./w2 +
        wc3.*accumarray(hcat(ICN[cell3], JCN[cell3]), varargin[:, vn][cell3].*wm3, sz=(Ny, Nx))./w3 +
        wc4.*accumarray(hcat(ICN[cell4], JCN[cell4]), varargin[:, vn][cell4].*wm4, sz=(Ny, Nx))./w4
        )./(wc1+wc2+wc4+wc4);
end

## EDGES

### top edge

topEdge = (jcn.>1) .& (jcn.<Nx) .& (icn.==1);
shifted = (jcn.<Nx-1) .& (icn.==1);

# cell 1

cell1 = shifted .& (quad.==2);
cell1 = vec(cell1');

ddx = dx[JCN[cell1].-1];
ddy = dy[1];
dxm = xm[cell1] .- xJCN[cell1];
dym = ym[cell1] .- yICN[cell1];
wm1 = 1 .- (dxm.*dym + (ddx.-dxm).*dym + (ddy.-dym).*dxm)./(ddx.*ddy);
w1 = accumarray(hcat(ICN[cell1], JCN[cell1]), wm1, sz=(Ny, Nx));

# cell 2 

cell2 = topEdge .& (quad.==1);
cell2 = vec(cell2');

ddx = dx[JCN[cell2]];
ddy = dy[1];
dxm = xm[cell2] .- xJCN[cell2];
dym = ym[cell2] .- yICN[cell2];
wm2 = 1 .- (dxm.*dym + (ddx.-dxm).*dym + (ddy.-dym).*dxm)./(ddx.*ddy);
w2  = accumarray(hcat(ICN[cell2], JCN[cell2]), wm2, sz=(Ny,Nx));

#loop over material properties to interpolate

for vn = 1:numV
    temp = (
        wc1*accumarray(hcat(ICN[cell1], JCN[cell1]), varargin[:, vn][cell1].*wm1, sz=(1, Nx))./w1 +
        wc2*accumarray(hcat(ICN[cell2], JCN[cell2]), varargin[:, vn][cell2].*wm2, sz=(1, Nx))./w2
        )/(wc1+wc2);
    n2interp[vn][1,2:end] = temp[1, 2:end];
end

### bottom edge

bottomEdge = (jcn.>1) .& (jcn.<Nx) .& (icn.==Ny-1);
shifted =    (jcn.<Nx-1) .& (icn.==Ny-1);

# cell 1

cell1 = shifted .& (quad.==3);
cell1 = vec(cell1');

ddx = dx[JCN[cell1].-1];
ddy = dy[Ny-1];
dxm = xm[cell1] .- xJCN[cell1];
dym = ym[cell1] .- y[end-1];
wm1 = 1 .- (dxm.*dym + (ddx.-dxm).*dym + (ddy.-dym).*dxm)./(ddx.*ddy);
w1 = accumarray(hcat(ones(Int, sum(cell1)), JCN[cell1]), wm1, sz=(1, Nx));

# cell 2

cell2 = bottomEdge .& (quad.==4);
cell2 = vec(cell2');

ddx = dx[JCN[cell2]];
ddy = dy[Ny-1];
dxm = xm[cell2] .- xJCN[cell2];
dym = ym[cell2] .- y[end-1];
wm2 = 1 .- (dxm.*dym + (ddx.-dxm).*dym + (ddy.-dym).*dxm)./(ddx.*ddy);
w2  = accumarray(hcat(ones(Int, sum(cell2)), JCN[cell2]), wm2, sz=(1,Nx));

#loop over material properties to interpolate

for vn = 1:numV
    temp = (
        wc1*accumarray(hcat(ones(Int, sum(cell1)), JCN[cell1]), varargin[:, vn][cell1].*wm1, sz=(1,Nx))./w1 +
        wc2*accumarray(hcat(ones(Int, sum(cell2)), JCN[cell2]), varargin[:, vn][cell2].*wm2, sz=(1,Nx))./w2
        )/(wc1+wc2);
    n2interp[vn][Ny,2:end] = temp[1, 2:end];
end

### left edge

leftEdge = (jcn.==1) .& (icn.>1) .& (icn.<Ny);
shifted  = (jcn.==1) .& (icn.<Ny-1);

# cell 1

cell1 = shifted .& (quad.==4);
cell1 = vec(cell1');

ddx = dx[1];
ddy = dy[ICN[cell1].-1];
dxm = xm[cell1] .- x[1];
dym = ym[cell1] .- yICN[cell1];
wm1 = 1 .- (dxm.*dym + (ddx.-dxm).*dym + (ddy.-dym).*dxm)./(ddx.*ddy);
w1 = accumarray(hcat(ICN[cell1], ones(Int, sum(cell1))), wm1, sz=(Ny,1));

# cell 2

cell2 = leftEdge .& (quad.==1);
cell2 = vec(cell2');

ddx = dx[1];
ddy = dy[ICN[cell2]];
dxm = xm[cell2] .- x[1];
dym = ym[cell2] .- yICN[cell2];
wm2 = 1 .- (dxm.*dym + (ddx.-dxm).*dym + (ddy.-dym).*dxm)./(ddx.*ddy);
w2 = accumarray(hcat(ICN[cell2], ones(Int, sum(cell2))), wm2, sz=(Ny,1));

#loop over material properties to interpolate

for vn = 1:numV
    temp = (
        wc1*accumarray(hcat(ICN[cell1], ones(Int, sum(cell1))), varargin[:, vn][cell1].*wm1, sz=(Ny,1))./w1 +
        wc2*accumarray(hcat(ICN[cell2], ones(Int, sum(cell2))), varargin[:, vn][cell2].*wm2, sz=(Ny,1))./w2
        )/(wc1+wc2);
    n2interp[vn][2:end-1, 1] = temp[2:end-1];
end

### right edge

rightEdge = (jcn.==Nx-1) .& (icn.>1) .& (icn.<Ny);
shifted = (jcn.==Nx-1) .& (icn.<Ny-1);

# cell 1

cell1 = shifted .& (quad.==3);
cell1 = vec(cell1');

ddx = dx[Nx-1];
ddy = dy[ICN[cell1].-1];
dxm = xm[cell1] .- x[Nx-1];
dym = ym[cell1] .- yICN[cell1];
wm1 = 1 .- (dxm.*dym + (ddx.-dxm).*dym + (ddy.-dym).*dxm)./(ddx.*ddy);
w1 = accumarray(hcat(ICN[cell1], ones(Int, sum(cell1))), wm1, sz=(Ny,1));

# cell 2

cell2 = rightEdge .& (quad.==2);
cell2 = vec(cell2');

ddx = dx[Nx-1];
ddy = dy[ICN[cell2]];
dxm = xm[cell2] .- x[Nx-1];
dym = ym[cell2] .- yICN[cell2];
wm2 = 1 .- (dxm.*dym + (ddx.-dxm).*dym + (ddy.-dym).*dxm)./(ddx.*ddy);
w2 = accumarray(hcat(ICN[cell2], ones(Int, sum(cell2))), wm2, sz=(Ny,1));

#loop over material properties to interpolate

for vn = 1:numV
    temp = (
        wc1*accumarray(hcat(ICN[cell1], ones(Int, sum(cell1))), varargin[:, vn][cell1].*wm1, sz=(Ny,1))./w1 +
        wc2*accumarray(hcat(ICN[cell2], ones(Int, sum(cell2))), varargin[:, vn][cell2].*wm2, sz=(Ny,1))./w2
        )/(wc1+wc2);
    n2interp[vn][2:end-1, Nx] = temp[2:end-1];
end

## CORNERS

# upper left

upperLeft = (jcn.==1) .& (icn.==1) .& (quad.==1);
upperLeft = vec(upperLeft');

ddx = dx[1];
ddy = dy[1];
dxm = xm[upperLeft] .- x[1];
dym = ym[upperLeft] .- y[1];
wm  = 1 .- (dxm.*dym + (ddx.-dxm).*dym + (ddy.-dym).*dxm)/(ddx*ddy);
wco = sum(wm);

for vn = 1:numV
    # Check what this line does...
    # n2interp[vn][1,1] = sum((varargin[vn]*upperLeft[vn]).*wm)./wco;
    n2interp[vn][1,1] = sum(varargin[:, vn][upperLeft].*wm)./wco; # this should be what I want...
end

# upper right

upperRight = (icn.==1) .& (jcn.==Nx-1) .& (quad.==2);
upperRight = vec(upperRight');

ddx = dx[Nx-1];
ddy = dy[1];
dxm = xm[upperRight] .- x[Nx-1];
dym = ym[upperRight] .- y[1];
wm  = 1 .- (dxm.*dym + (ddx.-dxm).*dym + (ddy.-dym).*dxm)/(ddx*ddy);
wco = sum(wm);

for vn = 1:numV
    # n2interp[vn][1,Nx] = sum((varargin[vn]*upperRight[vn]).*wm)./wco;
    n2interp[vn][1,Nx] = sum(varargin[:, vn][upperRight].*wm)./wco;
end

# lower Right

lowerRight = (icn.==Ny-1) .& (jcn.==Nx-1) .& (quad.==3);
lowerRight = vec(lowerRight');

ddx = dx[Nx-1];
ddy = dy[Ny-1];
dxm = xm[lowerRight] .- x[Nx-1];
dym = ym[lowerRight] .- y[Ny-1];
wm  = 1 .- (dxm.*dym + (ddx.-dxm).*dym + (ddy.-dym).*dxm)/(ddx*ddy);
wco = sum(wm);

for vn = 1:numV
    n2interp[vn][Ny,Nx] = sum(varargin[:, vn][lowerRight].*wm)./wco;
end

# lower left

lowerLeft = (icn.==Ny-1) .& (jcn.==1) .& (quad.==4);
lowerLeft = vec(lowerLeft');

ddx = dx[1];
ddy = dy[Ny-1];
dxm = xm[lowerLeft] .- x[1];
dym = ym[lowerLeft] .- y[Ny-1];
wm  = 1 .- (dxm.*dym + (ddx.-dxm).*dym + (ddy.-dym).*dxm)/(ddx*ddy);
wco = sum(wm);

for vn = 1:numV
    n2interp[vn][Ny,1] = sum(varargin[:, vn][lowerLeft].*wm)./wco;
end


return n2interp
end
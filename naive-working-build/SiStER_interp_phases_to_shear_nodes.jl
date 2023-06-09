using Interpolations

function SiStER_interp_phases_to_shear_nodes(xm,ym,icn,jcn,quad,x,y,phases,maxPhases)
# [phasesS] = SiStER_interp_phases_to_shear_nodes(xm,ym,icn,jcn,quad,x,y,phases,maxPhases)
# B.Z. Klein, July 2017, an interp function specific to phase [to shear nodes], to enable 
# exact mixing of several phases

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


# JCN = interp1(x, 1:length(x), xm, "nearest', 'extrap"); ## these are like the jcn & icn elsewhere, except the nodes are centered instead of upper left.
# ICN = interp1(y, 1:length(y), ym, "nearest', 'extrap"); ## this makes a lot of the indexing much simpler below.
JCN = interpolate((x,), 1:length(x), Gridded(Constant()))(xm);
ICN = interpolate((y,), 1:length(y), Gridded(Constant()))(ym);

phasesS = zeros(Ny, Nx, maxPhases);

## Interior Cells

center = (jcn.>1) .& (jcn.<Nx) .& (icn.>1) .& (icn.<Ny);
shiftLeft = (jcn.<Nx-1) .& (icn.>1) .& (icn.<Ny);
shiftUp = (jcn.>1) .& (jcn.<Nx) .& (icn.<Ny-1);
shiftBoth = (jcn.<Nx-1) .& (icn.<Ny-1);

# These lines are SLOW --> replace with loop? sped up by reducing allocations
# cell1 = center .& (((xm.-x[JCN]) .> 0) .& ((ym .- y[ICN]) .> 0));  ## these are logical arrays that index the original quadrants
# cell2 = shiftLeft .& (((xm.-x[JCN]) .< 0) .& ((ym .- y[ICN]) .> 0));
# cell3 = shiftBoth .& (((xm.-x[JCN]) .< 0) .& ((ym .- y[ICN]) .< 0));
# cell4 = shiftUp .& (((xm.-x[JCN]) .> 0) .& ((ym .- y[ICN]) .< 0));

cell1 = BitArray(undef, size(center));
cell2 = BitArray(undef, size(shiftLeft));
cell3 = BitArray(undef, size(shiftBoth));
cell4 = BitArray(undef, size(shiftUp));

xJCN = x[JCN]
yICN = y[ICN]
for i ∈ 1:length(cell1)
    cell1[i] = center[i] & ((xm[i]-xJCN[i]) > 0) & ((ym[i]-yICN[i]) > 0);
    cell2[i] = shiftLeft[i] & ((xm[i]-xJCN[i]) < 0) & ((ym[i] .- yICN[i]) > 0);
    cell3[i] = shiftBoth[i] & ((xm[i]-xJCN[i]) < 0) & ((ym[i] .- yICN[i]) < 0);
    cell4[i] = shiftUp[i] & ((xm[i]-xJCN[i]) > 0) & ((ym[i] .- yICN[i]) < 0);
end

# Transpose because arrays are (x,) and convert from BitArray to BitVector
cell1 = vec(cell1');
cell2 = vec(cell2');
cell3 = vec(cell3');
cell4 = vec(cell4');

### WEIGHTING [equal for now because that is what I'm running]

wc1 = 0.25;
wc2 = 0.25;
wc3 = 0.25;
wc4 = 0.25;

# cell 1 [i,j,1]

dxm = xm[cell1] .- xJCN[cell1];
dym = ym[cell1] .- yICN[cell1];
ddx = dx[JCN[cell1]];
ddy = dy[ICN[cell1]];

wm1 = 1 .- (dxm.*dym + (ddx-dxm).*dym + (ddy-dym).*dxm)./(ddx.*ddy);
w1 = accumarray(hcat(ICN[cell1], JCN[cell1]), wm1, sz=(Ny, Nx));

# cell 2 [i, j-1, 2]

dxm = xm[cell2] - xJCN[cell2];
dym = ym[cell2] - yICN[cell2];
ddx = dx[JCN[cell2].-1];
ddy = dy[ICN[cell2]];

wm2 = 1 .- (dxm.*dym + (ddx-dxm).*dym + (ddy-dym).*dxm)./(ddx.*ddy);
w2 = accumarray(hcat(ICN[cell2], JCN[cell2]), wm2, sz=(Ny, Nx));

# cell 3 [i-1, j-1, 3]

dxm = xm[cell3] - xJCN[cell3];
dym = ym[cell3] - yICN[cell3];
ddx = dx[JCN[cell3].-1];
ddy = dy[ICN[cell3].-1];

wm3 = 1 .- (dxm.*dym + (ddx-dxm).*dym + (ddy-dym).*dxm)./(ddx.*ddy);
w3 = accumarray(hcat(ICN[cell3], JCN[cell3]), wm3, sz=(Ny, Nx));

# cell 4 [i-1, j, 4]

dxm = xm[cell4] - xJCN[cell4];
dym = ym[cell4] - yICN[cell4];
ddx = dx[JCN[cell4]];
ddy = dy[ICN[cell4].-1];

wm4 = 1 .- (dxm.*dym + (ddx-dxm).*dym + (ddy-dym).*dxm)./(ddx.*ddy);
w4 = accumarray(hcat(ICN[cell4], JCN[cell4]), wm4, sz=(Ny, Nx));

#loop over material properties to interpolate

for n = 1:maxPhases
    
    phaseMask = phases .== n;
    
    phasesS[:,:,n] = (
        wc1.*accumarray(hcat(ICN[cell1 .& phaseMask], JCN[cell1 .& phaseMask]), wm1[phaseMask[cell1]], sz=(Ny,Nx))./w1 + 
        wc2*accumarray(hcat(ICN[cell2 .& phaseMask], JCN[cell2 .& phaseMask]), wm2[phaseMask[cell2]], sz=(Ny, Nx))./w2 + 
        wc3*accumarray(hcat(ICN[cell3 .& phaseMask], JCN[cell3 .& phaseMask]), wm3[phaseMask[cell3]], sz=(Ny, Nx))./w3 + 
        wc4*accumarray(hcat(ICN[cell4 .& phaseMask], JCN[cell4 .& phaseMask]), wm4[phaseMask[cell4]], sz=(Ny, Nx))./w4
        )./(wc1+wc2+wc4+wc4);
end



## EDGES

### top edge

topEdge = (jcn.>1) .& (jcn.<Nx) .& (icn.==1);
shifted = (jcn.<Nx.-1) .& (icn.==1);

# cell 1

cell1 = shifted .& (quad.==2);
cell1 = vec(cell1');

ddx = dx[JCN[cell1].-1];
ddy = dy[1];
dxm = xm[cell1] .- xJCN[cell1];
dym = ym[cell1] .- yICN[cell1];
wm1 = 1 .- (dxm.*dym + (ddx-dxm).*dym + (ddy.-dym).*dxm)./(ddx.*ddy);
w1 = accumarray(hcat(ICN[cell1], JCN[cell1]), wm1, sz=(1, Nx));

# cell 2 

cell2 = topEdge .& (quad.==1);
cell2 = vec(cell2');

ddx = dx[JCN[cell2]];
ddy = dy[1];
dxm = xm[cell2] - xJCN[cell2];
dym = ym[cell2] - yICN[cell2];
wm2 = 1 .- (dxm.*dym + (ddx-dxm).*dym + (ddy.-dym).*dxm)./(ddx.*ddy);
w2  = accumarray(hcat(ICN[cell2], JCN[cell2]), wm2, sz=(1, Nx));

#loop over material properties to interpolate

for n = 1:maxPhases
    
    phaseMask = phases .== n;
    
    temp = (
        wc1*accumarray(hcat(ICN[cell1 .& phaseMask], JCN[cell1 .& phaseMask]), wm1[phaseMask[cell1]], sz=(1, Nx))./w1 +
        wc2*accumarray(hcat(ICN[cell2 .& phaseMask], JCN[cell2 .& phaseMask]), wm2[phaseMask[cell2]], sz=(1, Nx))./w2
        )/(wc1+wc2);
    phasesS[1,:,n] = temp;
end

### bottom edge

bottomEdge = (jcn.>1) .& (jcn.<Nx) .& (icn.==Ny.-1);
shifted = (jcn.<Nx.-1) .& (icn.==Ny.-1);

# cell 1

cell1 = shifted .& (quad.==3);
cell1 = vec(cell1');

ddx = dx[JCN[cell1].-1];
ddy = dy[Ny-1];
dxm = xm[cell1] - xJCN[cell1];
dym = ym[cell1] .- y[end-1];
wm1 = 1 .- (dxm.*dym + (ddx.-dxm).*dym + (ddy.-dym).*dxm)./(ddx.*ddy);
w1 = accumarray(hcat(ones(Int, sum(cell1)), JCN[cell1]), wm1, sz=(1, Nx));

# cell 2

cell2 = bottomEdge .& (quad.==4);
cell2 = vec(cell2');

ddx = dx[JCN[cell2]];
ddy = dy[Ny-1];
dxm = xm[cell2] .- x[JCN[cell2]];
dym = ym[cell2] .- y[end-1];
wm2 = 1 .- (dxm.*dym + (ddx.-dxm).*dym + (ddy.-dym).*dxm)./(ddx.*ddy);
w2  = accumarray(hcat(ones(Int, sum(cell2)), JCN[cell2]), wm2, sz=(1, Nx));

#loop over material properties to interpolate

for n = 1:maxPhases
    
    phaseMask = phases .== n;
    
    temp = (
        wc1*accumarray(hcat(ones(Int, sum(cell1 .& phaseMask)), JCN[cell1 .& phaseMask]), wm1[phaseMask[cell1]], sz=(1, Nx))./w1 +
        wc2*accumarray(hcat(ones(Int, sum(cell2 .& phaseMask)), JCN[cell2 .& phaseMask]), wm2[phaseMask[cell2]], sz=(1, Nx))./w2
        )/(wc1+wc2);
    phasesS[Ny,:,n] = temp;
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
w1 = accumarray(hcat(ICN[cell1], ones(Int, sum(cell1))), wm1, sz=(Ny, 1));

# cell 2

cell2 = leftEdge .& (quad.==1);
cell2 = vec(cell2');

ddx = dx[1];
ddy = dy[ICN[cell2]];
dxm = xm[cell2] .- x[1];
dym = ym[cell2] .- yICN[cell2];
wm2 = 1 .- (dxm.*dym + (ddx.-dxm).*dym + (ddy.-dym).*dxm)./(ddx.*ddy);
w2 = accumarray(hcat(ICN[cell2], ones(Int, sum(cell2))), wm2, sz=(Ny, 1));

#loop over material properties to interpolate

for n = 1:maxPhases
    
    phaseMask = phases .== n;
    
    temp = (
        wc1*accumarray(hcat(ICN[cell1 .& phaseMask], ones(Int, sum(cell1 .& phaseMask))), wm1[phaseMask[cell1]], sz=(Ny, 1))./w1 +
        wc2*accumarray(hcat(ICN[cell2 .& phaseMask], ones(Int, sum(cell2 .& phaseMask))), wm2[phaseMask[cell2]], sz=(Ny, 1))./w2
        )/(wc1+wc2);
    phasesS[:, 1, n] = temp;
end

### right edge

rightEdge = (jcn.==Nx-1) .& (icn.>1) .& (icn.<Ny);
shifted =   (jcn.==Nx-1) .& (icn.<Ny-1);

# cell 1

cell1 = shifted .& (quad.==3);
cell1 = vec(cell1');

ddx = dx[Nx-1];
ddy = dy[ICN[cell1].-1];
dxm = xm[cell1] .- x[Nx-1];
dym = ym[cell1] .- yICN[cell1];
wm1 = 1 .- (dxm.*dym + (ddx.-dxm).*dym + (ddy.-dym).*dxm)./(ddx.*ddy);
w1 = accumarray(hcat(ICN[cell1], ones(Int, sum(cell1))), wm1, sz=(Ny, 1));

# cell 2

cell2 = rightEdge .& (quad.==2);
cell2 = vec(cell2');

ddx = dx[Nx-1];
ddy = dy[ICN[cell2]];
dxm = xm[cell2] .- x[Nx-1];
dym = ym[cell2] .- yICN[cell2];
wm2 = 1 .- (dxm.*dym + (ddx.-dxm).*dym + (ddy.-dym).*dxm)./(ddx.*ddy);
w2 = accumarray(hcat(ICN[cell2], ones(Int, sum(cell2))), wm2, sz=(Ny, 1));

#loop over material properties to interpolate

for n = 1:maxPhases
    
    phaseMask = phases .== n;
    
    temp = (
        wc1*accumarray(hcat(ICN[cell1 .& phaseMask], ones(Int, sum(cell1 .& phaseMask))), wm1[phaseMask[cell1]], sz=(Ny, 1))./w1 +
        wc2*accumarray(hcat(ICN[cell2 .& phaseMask], ones(Int, sum(cell2 .& phaseMask))), wm2[phaseMask[cell2]], sz=(Ny, 1))./w2
        )/(wc1+wc2);
    phasesS[:,Nx, n] = temp;
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

for n = 1:maxPhases
    
    phaseMask = phases .== n;
    
    phasesS[1,1,n] = sum(wm[phaseMask[upperLeft]])./wco;
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

for n = 1:maxPhases
    
    phaseMask = phases .== n;
    
    phasesS[1,Nx,n] = sum(wm[phaseMask[upperRight]])./wco;
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

for n = 1:maxPhases
    
    phaseMask = phases .== n;
    
    phasesS[Ny,Nx,n] = sum(wm[phaseMask[lowerRight]])./wco;
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

for n = 1:maxPhases
    
    phaseMask = phases .== n;
    
    phasesS[Ny,1,n] = sum(wm[phaseMask[lowerLeft]])./wco;
end



return phasesS
end
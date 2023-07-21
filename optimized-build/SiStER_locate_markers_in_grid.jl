using Interpolations

function SiStER_locate_markers_in_grid(xm,ym,x,y,dx,dy)
# [quad,icn,jcn] = SiStER_locate_markers_in_grid(xm,ym,x,y,dx,dy)
# Tells a marker which cell (and which quadrant of that cell) it belongs to.
    # icn;jcn are the indexes of the upper-left shear node of the cell()
    # that a given marker is currently in
    # quad is the quadrant of the cell that contains the marker
    # (quad = 1 means bottom-right, then numbered clockwise)
    # sped up by B. Klein in Fall 2016 by using interp1 function

# print("% Free memory: "*string(Sys.free_memory()/Sys.total_memory())*"\n")
## Determine Location of Markers & Quadrant of Element 
M=length(xm);
icn = zeros(1,M);
jcn = zeros(1,M);

indX = 1:length(x);
indY = 1:length(y);

# Replacement --> use Julia 'broadcast'
# Ix = minimum.(abs.(broadcast(-, xm, x')));
# Iy = minimum.(abs.(broadcast(-, ym, y')));

# Need to deal with out of bounds indexes --> take closest in bounds value
Ix = extrapolate(interpolate((x,), indX, Gridded(Constant())), Flat())(xm);
Iy = extrapolate(interpolate((y,), indY, Gridded(Constant())), Flat())(ym);

Ix = Int.(Ix)
Iy = Int.(Iy)

# Ix = interp1(x, indX, xm, "nearest", "extrap"); # This is speedup over bsxfun
# Iy = interp1(y, indY, ym, "nearest", "extrap");

# [~, Ix] = min(abs(bsxfun[@minus, xm, x']));
# [~, Iy] = min(abs(bsxfun[@minus, ym, y']));

jcn[xm.>x[Ix]]  = Ix[xm.>x[Ix]];
jcn[xm.<=x[Ix]] = Ix[xm.<=x[Ix]].-1;
jcn = Int.(jcn);

icn[ym.>y[Iy]]  = Iy[ym.>y[Iy]];
icn[ym.<=y[Iy]] = Iy[ym.<=y[Iy]].-1;
icn = Int.(icn);

jcn[jcn.==0] .= 1;
jcn[jcn.>length(dx)] .= length(dx);

icn[icn.==0] .= 1;
icn[icn.>length(dy)] .= length(dy);

# print("% Free memory: "*string(Sys.free_memory()/Sys.total_memory())*"\n")

# out of memory error here --> Do in loop
# disx = abs.((xm.-x[jcn])./dx[jcn]);
# disy = abs.((ym-y[icn])./dy[icn]);
disx = zeros(size(jcn));
for ind ∈ 1:length(jcn)
    jcnVal = jcn[ind]
    disx[ind] = abs((xm[ind]-x[jcnVal])/dx[jcnVal])
end
disy = zeros(size(icn));
for ind ∈ 1:length(icn)
    icnVal = icn[ind]
    disy[ind] = abs((ym[ind]-y[icnVal])/dy[icnVal])
end

xRIGHT = disx .> 0.5;
yUP = disy .> 0.5;

quad = zeros(1,M); # quadrant 1 = bottom-right; numbered clockwise
quad[xRIGHT .& yUP] .= 3;
quad[xRIGHT .& .!yUP] .= 2;
quad[.!xRIGHT .& yUP] .= 4;
quad[.!xRIGHT .& .!yUP] .= 1;

return Int.(quad), Int.(icn), Int.(jcn)
end
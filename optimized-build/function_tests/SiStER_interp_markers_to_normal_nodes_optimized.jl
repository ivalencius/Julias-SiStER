using RecursiveArrayTools

function meshgrid(x, y)
    X = [i for i in x, _ in 1:length(y)]
    Y = [j for _ in 1:length(x), j in y]
    return X', Y'
end

function SiStER_interp_markers_to_normal_nodes_optimized(xm,ym,icn,jcn,x,y,varargin)
# [n2interp] = SiStER_interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,varargin)
# interpolates properties [in the order of input] from markers to normal nodes
#
# First cut - J.A. Olive; March 2011
# Modified by E. Mittelstaedt; April 2011 to allow multiple inputs  
# Modified by B.Z. Klein; Spring 2014 for speedup

Nx=length(x);
Ny=length(y);
dx=diff(x);
dy=diff(y);

# check for number of properties to interpolate
numV = size(varargin,2);
# n2interp[1:numV] = struct("data",zeros(Ny,Nx));
n2interp = VectorOfArray([zeros(Ny, Nx) for _ in 1:numV]);

# INDEX = sub2ind([Ny-1, Nx-1], icn, jcn);
s2i = LinearIndices(zeros(Ny-1, Nx-1)); # reverse is CartesianIndices
INDEX = zeros(length(icn));
for ind ∈ 1:length(icn)
    INDEX[ind] = s2i[icn[ind], jcn[ind]]
end
INDEX = Int.(INDEX);


# AcCell = bsxfun[@times, dy', dx];
AcCell = dx' .* dy

xN = x[1:Nx-1] .+ dx/2;
yN = y[1:Ny-1] .+ dy/2;
XN, YN = meshgrid(xN, yN);

AMvec = abs.((xm .- XN[INDEX]).*(ym .- YN[INDEX]));
WMvec = (AcCell[INDEX] .- AMvec)./AcCell[INDEX];

# From SiStER_material_props_on_nodes.jl line 43 --> this shouldnt have all 0's at bottom
# w = accumarray(hcat(icn', jcn'), WMvec, sz=(Ny-1, Nx-1));
w = zeros((Ny-1, Nx-1));
counter = 1
for (i, j) in zip(icn', jcn')
    val = WMvec[counter]
    w[i, j] += isnan(val) ? 0 : val
    counter+=1;
end

for vn = 1:numV
    VecData = varargin[:, vn].*WMvec;

    # accumarray(hcat(icn', jcn'), VecData, sz=(Ny-1, Nx-1))
    arr = zeros((Ny-1, Nx-1));
    counter = 1
    for (i, j) in zip(icn', jcn')
        val = VecData[counter]
        arr[i, j] += isnan(val) ? 0 : val
        counter+=1;
    end

    n2interp[vn][2:Ny,2:Nx] = arr./w;
end

return n2interp
end
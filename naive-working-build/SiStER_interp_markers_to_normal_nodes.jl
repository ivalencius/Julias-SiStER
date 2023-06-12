using RecursiveArrayTools
using GroupSlices

function SiStER_interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,varargin)
# [n2interp] = SiStER_interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,varargin)
# interpolates properties [in the order of input] from markers to normal nodes
#
# First cut - J.A. Olive; March 2011
# Modified by E. Mittelstaedt; April 2011 to allow multiple inputs  
# Modified by B.Z. Klein; Spring 2014 for speedup

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

function meshgrid(x, y)
    X = [i for i in x, _ in 1:length(y)]
    Y = [j for _ in 1:length(x), j in y]
    return X', Y'
end

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
w = accumarray(hcat(icn', jcn'), WMvec, sz=(Ny-1, Nx-1));

for vn = 1:numV
    VecData = varargin[:, vn].*WMvec;
    n2interp[vn][2:Ny,2:Nx] = accumarray(hcat(icn', jcn'), VecData, sz=(Ny-1, Nx-1))./w;
end



return n2interp
end
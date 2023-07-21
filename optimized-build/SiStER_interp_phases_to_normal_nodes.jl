using GroupSlices

function SiStER_interp_phases_to_normal_nodes(xm,ym,icn,jcn,x,y,phases, maxPhases)
# [phaseWeights] = SiStER_interp_phases_to_normal_nodes(xm,ym,icn,jcn,x,y,phases, maxPhases)
# B.Z. Klein, July 2017, an interp function specific to phase [for normal nodes], to enable 
# exact mixing of several phases

function meshgrid(x, y)
    X = [i for i in x, _ in 1:length(y)]
    Y = [j for _ in 1:length(x), j in y]
    return X', Y'
end

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

phaseWeights = zeros(Ny, Nx, maxPhases);

for n = 1:maxPhases
    
    phaseMask = phases.==n;
    phaseWeights[2:Ny, 2:Nx,n] = accumarray(
        hcat(icn[phaseMask'], jcn[phaseMask']),
        WMvec[phaseMask]', 
        sz=size(AcCell));
    
end
   
sumWeights = repeat(sum(phaseWeights, dims = 3), outer=[1, 1, maxPhases]);
phaseWeights = phaseWeights./sumWeights;


return phaseWeights
end
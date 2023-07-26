function meshgrid(x, y)
    X = [i for i in x, _ in 1:length(y)]
    Y = [j for _ in 1:length(x), j in y]
    return X', Y'
end

function SiStER_interp_phases_to_normal_nodes_optimized(xm,ym,icn,jcn,x,y,phases, maxPhases)
# [phaseWeights] = SiStER_interp_phases_to_normal_nodes(xm,ym,icn,jcn,x,y,phases, maxPhases)
# B.Z. Klein, July 2017, an interp function specific to phase [for normal nodes], to enable 
# exact mixing of several phases

# @btime begin
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

AcCell = dx' .* dy

xN = x[1:Nx-1] .+ dx/2;
yN = y[1:Ny-1] .+ dy/2;
XN, YN = meshgrid(xN, yN);

AMvec = abs.((xm .- XN[INDEX]).*(ym .- YN[INDEX]));
WMvec = (AcCell[INDEX] .- AMvec)./AcCell[INDEX];
# end # 11ms

# @btime begin # BIGGEST SLOWDOWN HERE
phaseWeights = zeros(Ny, Nx, maxPhases);

for n = 1:maxPhases
    
    phaseMask = phases.==n;
    tmp = zeros(size(AcCell));
    WMmasked = WMvec[phaseMask];
    counter = 1
    for (i, j) in zip(icn[phaseMask'], jcn[phaseMask'])
        val = WMmasked[counter]
        tmp[i, j] += isnan(val) ? 0 : val
        counter+=1;
    end
    phaseWeights[2:Ny, 2:Nx,n] = tmp;
end
# end # 115 ms

# @btime begin 
sumWeights = repeat(sum(phaseWeights, dims = 3), outer=[1, 1, maxPhases]);
phaseWeights = phaseWeights./sumWeights;
# end # 20 μs
return phaseWeights
end
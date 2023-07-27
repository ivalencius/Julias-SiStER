# Test the cost of reimplementation of accumarray in Julia

using BenchmarkTools
using GroupSlices
using JLD2
using Test

# StackOverflow solution
function stack_accumarray(subs, val, fun=sum, fillval=0; sz=maximum(subs,1), issparse=false)
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

# My tailored implementation
function SiStER_accumarray(index, data; sz=size(data), fun=sum)
    tmp = zeros(sz);
    counter = 1
    for (i, j) in zip(index[:,1], index[:,2])
        val = data[i]
        tmp[i, j] == isnan(val) ? tmp[i,j] : fun(sum)+tmp[i,j]
        counter+=1;
    end
    return tmp
end

function stack_test(phaseWeights, phases, maxPhases, icn, jcn, WMvec, AcCell, Nx, Ny)
    for n = 1:maxPhases
        phaseMask = phases.==n;
        phaseWeights[2:Ny, 2:Nx, n] = stack_accumarray(
            hcat(icn[phaseMask'], jcn[phaseMask']),
            WMvec[phaseMask]',sz=size(AcCell));
    end
    return phaseWeights
end

function homebrew_test1(phaseWeights, phases, maxPhases, icn, jcn, WMvec, AcCell, Nx, Ny)
    for n = 1:maxPhases
        phaseMask = phases.==n;
        phaseWeights[2:Ny, 2:Nx, n] = SiStER_accumarray(
            hcat(icn[phaseMask'], jcn[phaseMask']),
            WMvec[phaseMask], sz=size(AcCell)); 
    end
    return phaseWeights
end

function homebrew_test2(phaseWeights, phases, maxPhases, icn, jcn, WMvec, AcCell, Nx, Ny)
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
    return phaseWeights
end
# Test it on equivalent of SiStER_interp_phases_to_normal_nodes.jl
function test_implementation()
    @load "test-vars.jld2" xm ym icn jcn x y im qd
    phases = im
    maxPhases = 3
    Nx=length(x);
    Ny=length(y);
    dx=diff(x);
    dy=diff(y);
    
    function meshgrid(x, y)
        X = [i for i in x, _ in 1:length(y)]
        Y = [j for _ in 1:length(x), j in y]
        return X', Y'
    end

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
    
    @test stack_test(phaseWeights, phases, maxPhases, icn, jcn, WMvec, AcCell, Nx, Ny) ≈ homebrew_test1(phaseWeights, phases, maxPhases, icn, jcn, WMvec, AcCell, Nx, Ny)
    @test stack_test(phaseWeights, phases, maxPhases, icn, jcn, WMvec, AcCell, Nx, Ny) ≈ homebrew_test2(phaseWeights, phases, maxPhases, icn, jcn, WMvec, AcCell, Nx, Ny)
    
    println("Stackoverflow Test")
    @btime stack_test($phaseWeights, $phases, $maxPhases, $icn, $jcn, $WMvec, $AcCell, $Nx, $Ny)

    println("\nTailored (inside function)")
    @btime homebrew_test1($phaseWeights, $phases, $maxPhases, $icn, $jcn, $WMvec, $AcCell, $Nx, $Ny)
  
    println("\nTailored (non-function)")
    @btime homebrew_test2($phaseWeights, $phases, $maxPhases, $icn, $jcn, $WMvec, $AcCell, $Nx, $Ny)
    
end

test_implementation()
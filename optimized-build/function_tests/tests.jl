using Test
using BenchmarkTools
using JLD2

include("SiStER_interp_phases_to_normal_nodes.jl")
include("Sister_interp_phases_to_normal_nodes_optimized.jl")
function test_interp_phases_to_normal_nodes()
    println("\nSiStER_interp_phases_to_normal_nodes")
    @load "test-vars.jld2" xm ym icn jcn x y im qd
    phases = im
    maxPhases = 3
    soln = SiStER_interp_phases_to_normal_nodes(xm,ym,icn,jcn,x,y,phases,maxPhases);
    soln_opt = SiStER_interp_phases_to_normal_nodes_optimized(xm,ym,icn,jcn,x,y,phases,maxPhases);
    @test isapprox(soln, soln_opt, nans=true)
    println("Outputs Equivalent Within Error ✔")
    println("Timings")
    print("\tOriginal: \t")
    @btime soln = SiStER_interp_phases_to_normal_nodes($xm,$ym,$icn,$jcn,$x,$y,$phases,$maxPhases);
    print("\tOptimized: \t")
    @btime soln_opt = SiStER_interp_phases_to_normal_nodes_optimized($xm,$ym,$icn,$jcn,$x,$y,$phases,$maxPhases);
end

include("SiStER_interp_phases_to_shear_nodes.jl")
include("Sister_interp_phases_to_shear_nodes_optimized.jl")
function test_interp_phases_to_shear_nodes()
    println("\nSiStER_interp_phases_to_shear_nodes")
    @load "test-vars.jld2" xm ym icn jcn x y im qd
    phases = im
    quad = qd
    maxPhases = 3
    soln = SiStER_interp_phases_to_shear_nodes(xm,ym,icn,jcn,quad,x,y,phases,maxPhases);
    soln_opt = SiStER_interp_phases_to_shear_nodes_optimized(xm,ym,icn,jcn,quad,x,y,phases,maxPhases);
    @test isapprox(soln, soln_opt, nans=true)
    println("Outputs Equivalent Within Error ✔")
    println("Timings")
    print("\tOriginal: \t")
    @btime soln = SiStER_interp_phases_to_shear_nodes($xm,$ym,$icn,$jcn,$quad,$x,$y,$phases,$maxPhases);
    print("\tOptimized: \t")
    @btime soln_opt = SiStER_interp_phases_to_shear_nodes_optimized($xm,$ym,$icn,$jcn,$quad,$x,$y,$phases,$maxPhases);
end

include("SiStER_interp_markers_to_normal_nodes.jl")
include("Sister_interp_markers_to_normal_nodes_optimized.jl")
function test_interp_markers_to_normal_nodes()
    println("\nSiStER_interp_markers_to_normal_nodes")
    @load "test-vars.jld2" xm ym icn jcn x y im qd
    varargin = rand(Float64, length(xm)).*100
    soln = SiStER_interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,varargin);
    soln_opt = SiStER_interp_markers_to_normal_nodes_optimized(xm,ym,icn,jcn,x,y,varargin);
    @test isapprox(soln, soln_opt, nans=true)
    println("Outputs Equivalent Within Error ✔")
    println("Timings")
    print("\tOriginal: \t")
    @btime soln = SiStER_interp_markers_to_normal_nodes($xm,$ym,$icn,$jcn,$x,$y,$varargin);
    print("\tOptimized: \t")
    @btime soln_opt = SiStER_interp_markers_to_normal_nodes_optimized($xm,$ym,$icn,$jcn,$x,$y,$varargin);
end

include("SiStER_interp_markers_to_shear_nodes.jl")
include("Sister_interp_markers_to_shear_nodes_optimized.jl")
function test_interp_markers_to_shear_nodes()
    println("\nSiStER_interp_markers_to_shear_nodes")
    @load "test-vars.jld2" xm ym icn jcn x y im qd
    quad = qd
    varargin = rand(Float64, length(xm)).*100
    soln = SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,quad,x,y,varargin);
    soln_opt = SiStER_interp_markers_to_shear_nodes_optimized(xm,ym,icn,jcn,quad,x,y,varargin);
    @test isapprox(soln, soln_opt, nans=true)
    println("Outputs Equivalent Within Error ✔")
    println("Timings")
    print("\tOriginal: \t")
    @btime soln = SiStER_interp_markers_to_shear_nodes($xm,$ym,$icn,$jcn,$quad,$x,$y,$varargin);
    print("\tOptimized: \t")
    @btime soln_opt = SiStER_interp_markers_to_shear_nodes_optimized($xm,$ym,$icn,$jcn,$quad,$x,$y,$varargin);
end

########################### Toggle tests on/off ################################
test_interp_phases_to_normal_nodes()
test_interp_phases_to_shear_nodes()
test_interp_markers_to_normal_nodes()
test_interp_markers_to_shear_nodes()

# Replacement code for accumarray
# tmp = zeros(size(AcCell));
# WMmasked = WMvec[phaseMask];
# counter = 1
# for (i, j) in zip(icn[phaseMask'], jcn[phaseMask'])
#     val = WMmasked[counter]
#     tmp[i, j] += isnan(val) ? 0 : val
#     counter+=1;
# end
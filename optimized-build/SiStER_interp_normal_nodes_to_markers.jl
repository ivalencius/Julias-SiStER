include("SiStER_interp_grid_to_marker_vector.jl")
using GroupSlices

function SiStER_interp_normal_nodes_to_markers(varN,xc,yc,xm,ym,icn,jcn)
# [varm]=SiStER_interp_normal_nodes_to_markers(varN,xc,yc,xm,ym,icn,jcn)
# interpolates properties from normal nodes to markers
# J.-A. Olive, 2011 [First cut]
# B.Z. Klein, 2014 [speedup]

Nx=length(xc)+1;
Ny=length(yc)+1;
(m, _) = size(varN);

ic = copy(icn);
jc = copy(jcn);

jc[xm' .> xc[jc]] = jc[xm' .> xc[jc]].+1;
jc[jc.<2] .= 2;
jc[jc.>(Nx-1)] .= Nx-1;

ic[ym'.>yc[ic]] = ic[ym'.>yc[ic]].+1;
ic[ic.<2] .= 2;
ic[ic.>(Ny-1)] .= Ny-1;

xNnodes = [xc[jc.-1]' xc[jc]' xc[jc]' xc[jc.-1]']';
yNnodes = [yc[ic.-1]' yc[ic.-1]' yc[ic]' yc[ic]']';

# IND = sub2ind(size(varN), ic, jc);
s2i = LinearIndices(size(varN)); # reverse is CartesianIndices
IND = zeros(length(ic));
for ind ∈ 1:length(ic)
    IND[ind] = s2i[ic[ind], jc[ind]]
end
IND = Int.(IND);
VARnodes = [varN[IND] varN[IND.+m] varN[IND.+(1+m)] varN[IND.+1]]';

varm = SiStER_interp_grid_to_marker_vector(xNnodes,yNnodes,VARnodes,xm,ym);


return varm

end
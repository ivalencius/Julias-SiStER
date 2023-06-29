include("SiStER_interp_grid_to_marker_vector.jl")
function SiStER_interp_shear_nodes_to_markers(varS,x,y,xm,ym,icn,jcn)
# [varm]=SiStER_interp_shear_nodes_to_markers(varS,x,y,xm,ym,icn,jcn)
# interpolates properties from shear nodes to markers

(m, _) = size(varS) 

xnodes = vcat(x[jcn], x[jcn.+1], x[jcn.+1], x[jcn]);
ynodes = vcat(y[icn], y[icn], y[icn.+1], y[icn.+1]);

# INDEX = sub2ind(size(varS), icn, jcn);
s2i = LinearIndices(size(varS)); # reverse is CartesianIndices
# foo = hcat(icn, jcn)
# INDEX = hcat(s2i[icn]', s2i[jcn]');
INDEX = zeros(length(icn));
for ind ∈ 1:length(icn)
    INDEX[ind] = s2i[icn[ind], jcn[ind]]
end
INDEX = Int.(INDEX);
VARnodes = [varS[INDEX] varS[INDEX.+m] varS[INDEX.+(m+1)] varS[INDEX.+1]]';

varm = SiStER_interp_grid_to_marker_vector(xnodes,ynodes,VARnodes,xm,ym);


return varm
end


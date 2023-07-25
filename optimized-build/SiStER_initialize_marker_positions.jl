include("SiStER_seed_markers_uniformly.jl")

function SiStER_initialize_marker_positions(xsize,ysize,dx,dy,Mquad)
# [xm, ym] = SiStER_initialize_marker_positions(xsize,ysize,dx,dy,Mquad)
# assigns markers their inital position [coordinates]
# markers are seeded by cell quadrant; to make sure there is enough of them
# in each quadrant 



# smallest quadrant sizes
mdx=minimum(dx)/2;
mdy=minimum(dy)/2;

# creating a regular grid with step = that of the smallest quadrant
xx=0:mdx:xsize;
yy=0:mdy:ysize;

nxx=length(xx);
nyy=length(yy);

xm = zeros((nxx-1)*(nyy-1)*Mquad);
ym = zeros((nxx-1)*(nyy-1)*Mquad);
midx=1;

# @time 0.139 sec (2.45 M allocations: 110.030 MiB, 50.30% gc time)
# Using vcat method: 112.716897 seconds (2.25 M allocations: 782.566 GiB, 37.28% gc time, 0.00% compilation time)
for i=1:nyy-1
    for j=1:nxx-1
        # @time 0.0649 sec (136.41 k allocations, 9.205 MiB, 99.93% compilation time)
        (xmtemp, ymtemp)=SiStER_seed_markers_uniformly(xx[j],yy[i],mdx,mdy,Mquad);
        xm[midx:midx+Mquad-1]=xmtemp;
        ym[midx:midx+Mquad-1]=ymtemp;
        midx=midx+Mquad;
        # xmtemp and ymtemp just fill in xm and ym - Ilan
        # xm = vcat(xm, xmtemp)
        # ym = vcat(ym, ymtemp)        
    end
end

return [xm, ym]
end
using Random
function SiStER_seed_markers_uniformly(x,y,dx,dy,N)
# [xm, ym]=SiStER_seed_markers_uniformly(x,y,dx,dy,N)
# seeds N markers in the cell | quadrant whose upper-left node 
# has coordinates [x,y] & width, height [dx,dy]
function meshgrid(x, y)
    X = [i for i in x, _ in 1:length(y)]
    Y = [j for _ in 1:length(x), j in y]
    return X', Y'
end

# create a subgrid

fact=0.4; # randomization factor()

nx=ceil(Int, sqrt(N*dx/dy));
ny=ceil(Int, N/nx);
actual_N=nx*ny;

ddx=range(x,x+dx,length=nx+1);
ddy=range(y,y+dy,length=ny+1);
(DX, DY)=meshgrid(ddx,ddy);

xm=zeros(1,actual_N);
ym=zeros(1,actual_N);

k=0;
for i=1:ny
    for j=1:nx
        
        k=k+1;
        xsub=DX[i,j];
        ysub=DY[i,j];
        dxsub=dx/nx;
        dysub=dy/ny;
       
        # randomize marker position around sub-cell center
        xm[k]=xsub+(dxsub/2)*(1+fact*2*(rand()-0.5));
        ym[k]=ysub+(dysub/2)*(1+fact*2*(rand()-0.5));

    end
end


if actual_N>N # only keep N random markers out of actual_N markers

    idx=randperm(actual_N);
    idx=idx[1:N];
    xm=xm[idx];
    ym=ym[idx];
    
end

return [xm, ym]
end
using RecursiveArrayTools
using StatsBase
include("SiStER_seed_markers_uniformly.jl")

function SiStER_patch_marker_holes(icn,jcn,quad,Nx,Ny,Mquad,Mquad_crit,xm,ym,x,y,dx,dy,im,ep,idm,Tm,sxxm,sxym,epNH,epsIIm)
# function [xm, ym, im, Ifix, mp, ep, idm, Tm, sxxm, sxym, epNH, epsIIm]=SiStER_patch_marker_holes(icn,jcn,quad,Nx,Ny,Mquad,Mquad_crit,xm,ym,x,y,dx,dy,im,ep,idm,Tm,sxxm,sxym,epNH,epsIIm)
#
# seeds new markers in all quadrants where marker density has fallen below
# the threshold 
# AND THEN assigns those new markers a marker index; a marker phase & a plastic
# strain based on the properties of the markers surrounding them
# (these particular properties are never passed on the grid, e.g. plastic strain...)
#
# J.-A. Olive; & B.Z. Klein; 2012-2014

M=length(xm);
md_crit=copy(Mquad_crit);
mark_per_quad=copy(Mquad);
mdx=minimum(dx)/2;
mdy=minimum(dy)/2;

############### LOOK FOR EMPTY [no more markers] QUADRANTS
# mp = accumarray(hcat(icn, jcn, Int.(quad)), 1, sz=(Ny-1, Nx-1, 4));
# Redo accumarray
mp = zeros(Ny-1, Nx-1, 4);
for (i, j, q) ∈ zip(icn, jcn, quad)
    mp[i, j, q] += 1
end

empty = findall(mp.==0);

display_message_if_empty=0;
if display_message_if_empty==1
    if(!isempty(empty))
        # [iEmpty,jEmpty,kEmpty] = ind2sub(size(mp), empty);
        # s2i = CartesianIndices(size(mp)); # reverse is CartesianIndices
        # IND = zeros(length(icn));
        # for ind ∈ 1:length(icn)
        #     IND[ind] = s2i[icn[ind], jcn[ind], quad[ind]]
        # end
        # IND = Int.(IND);
        # iEmpty, jEmpty, kEmpty = getindex.(empty, [1 2 3])
        temp_inds = getindex.(empty, [1 2 3]);
        iEmpty = temp_inds[:, 1];
        jEmpty = temp_inds[:, 2];
        kEmpty = temp_inds[:, 3];
        for n = 1:min(length(iEmpty), 100)
            println("WARNING ! Empty quadrant number "*string(kEmpty[n])*" in cell i = "*string(iEmpty[n])*" j = "*string(jEmpty[n]))
        end
    end
end
    
    
##### LOCATE QUADRANTS WHERE MARKER DENSITY IS BELOW THRESHOLD
mpCrInd = findall(mp.<md_crit);
# iicr, jjcr, qqcr = ind2sub(size(mp), mpCrInd);
if !isempty(mpCrInd)
    temp_inds = getindex.(mpCrInd, [1 2 3]);
    iicr = temp_inds[:, 1]';
    jjcr = temp_inds[:, 2]';
    qqcr = temp_inds[:, 3]';
else
    iicr, jjcr, qqcr = [], [], [];
end
    
# NEED TO SEED NEW MARKERS IN THOSE QUADRANTS 
# SO THAT WE ARE BACK TO THE INITIAL MARKER DENSITY IN THOSE QUADRANTS
    
xrsd=[] 
yrsd=[] 
im_fix=[]; # marker phase
ep_fix=[]; # plastic strain
epNH_fix=[]; # non-healed plastic strain
te_fix=[]; # temperature
sxx_fix=[]; # stress
sxy_fix=[]; # stress
sr_fix=[]; # strain rate
   
if !isempty(iicr) # if there are critical quadrants
         
    for c ∈ 1:length(iicr) # go through all critical quadrants

        # the upper-left node of the corresponding cell is()
        icell=iicr[c];
        jcell=jjcr[c];
        qcell=qqcr[c];
        # the quadrant area in that cell is() 
        qsize=0.25*dx[jcell]*dy[icell];
        # if the smallest quadrant in the domain [area mdx*mdy] has
        # mark_per_quad markers; then this critical quadrant needs
        Nfix=ceil(qsize/(mdx*mdy))*mark_per_quad;
        Nfix=max(Nfix-mp[iicr[c],jjcr[c],qqcr[c]],1);
        Nfix = Int.(Nfix);
        # find the coordinates of the upper-left corner of the quadrant
        if qcell==1 || qcell==4
            xcorn=x[jcell];
        else
            xcorn=x[jcell]+(dx[jcell]/2);
        end

        if qcell==1 || qcell==2
            ycorn=y[icell];
        else
            ycorn=y[icell]+(dy[icell]/2);
        end
            
        # draw random marker location
        xmrr, ymrr=SiStER_seed_markers_uniformly(xcorn,ycorn,dx[jcell]./2 , dy[icell]./2 ,Nfix)
            
        xrsd=append!(xrsd, xmrr);
        yrsd=append!(yrsd, ymrr);
            
        #### NOW THAT THE NEW MARKERS ARE SEEDED;
        ##### ASSIGN PARAMETERS THAT ARE NEVER STORED ON THE EULERIAN GRID
        
        # the value is assigned based on the average; | max. value of the
        # markers that remain in the corresponding CELL [since there's no grid value to interpolate from]
                    
        # CAREFUL THIS CANNOT WORK IF WE END UP WITH AN EMPTY CELL
        # if that was to happen; let's just draw "im" randomly
            
                       
        if isempty((ep[vec((icn.==icell) .& (jcn.==jcell))]))
            println("WARNING ! - EMPTY CELL - SOMETHING IS VERY WRONG...")    
            im_fix=1+floor(rand(1,Nfix)*maximum(im)); # random phase number
            ep_fix=zeros(1,Nfix);
            epNH_fix=zeros(1,Nfix);
            te_fix=zeros(1,Nfix);
            sxx_fix=zeros(1,Nfix);
            sxy_fix=zeros(1,Nfix);
            sr_fix=zeros(1,Nfix);
        else
            # assign the average phase of the markers that are left in the cell()
            phase_fix=round(mode(im[vec((icn.==icell) .& (jcn.==jcell))]));
            # assign the greatest plastic strain of the markers that are left in the cell()
            strain_fix=maximum((ep[vec((icn.==icell) .& (jcn.==jcell))]));
            strainNH_fix=maximum((epNH[vec((icn.==icell) .& (jcn.==jcell))]));
            # assign the average temperature of the markers that are left in
            # the cell()
            temp_fix=mean((Tm[vec((icn.==icell) .& (jcn.==jcell))]));
            # reassign mean stress / strain rate of markers left in cell()
            stress_xx_fix=mean((sxxm[vec((icn.==icell) .& (jcn.==jcell))]));
            stress_xy_fix=mean((sxym[vec((icn.==icell) .& (jcn.==jcell))]));
            strainrate_fix=mean((epsIIm[vec((icn.==icell) .& (jcn.==jcell))]));
        end

        im_fix=append!(im_fix, phase_fix.*ones(Nfix))
        ep_fix=append!(ep_fix, strain_fix.*ones(1,Nfix))
        epNH_fix=append!(epNH_fix, strainNH_fix.*ones(1,Nfix))
        te_fix=append!(te_fix, temp_fix.*ones(1,Nfix))
        sxx_fix=append!(sxx_fix, stress_xx_fix.*ones(1,Nfix))
        sxy_fix=append!(sxy_fix, stress_xy_fix.*ones(1,Nfix)) 
        sr_fix=append!(sr_fix, strainrate_fix.*ones(1,Nfix));
    
    end


       
    # NOW ASSIGN PROPERTIES TO THOSE MARKERS  
    Npatch=length(xrsd);
    index_fix=maximum(idm)+1:maximum(idm)+Npatch;

    Ifix=M+1:1:M+Npatch; # total number of markers added to fix holes in critical quadrants
    # xm[Ifix]=xrsd;
    # ym[Ifix]=yrsd;
    # im[Ifix]=im_fix;
    # ep[Ifix]=ep_fix;
    # epNH[Ifix]=epNH_fix;
    # idm[Ifix]=index_fix;
    # Tm[Ifix]=te_fix;
    # sxxm[Ifix]=sxx_fix;
    # sxym[Ifix]=sxy_fix;
    # epsIIm[Ifix]=sr_fix;
    xm = append!(xm, xrsd)
    ym = append!(ym, yrsd)
    im = append!(im, im_fix)
    ep = append!(ep, ep_fix)
    epNH = append!(epNH, epNH_fix)
    idm = append!(idm, index_fix)
    Tm = append!(Tm, te_fix)
    sxxm = append!(sxxm, sxx_fix)
    sxym = append!(sxym, sxy_fix)
    epsIIm = append!(epsIIm, sr_fix)

    # uncomment to display number of added markers
    #@sprintf("\n%d%s%d%s\n", length(Ifix), ' markers added in ', length(iicr), " cell quadrants.")
    
    else
            
    Ifix=0;
           
end

return xm, ym, im, Ifix, mp, ep, idm, Tm, sxxm, sxym, epNH, epsIIm

end
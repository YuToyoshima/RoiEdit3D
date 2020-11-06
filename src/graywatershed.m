function ret = graywatershed(im,pos,thr)
%%% This function segments image using pre-defined seed (pos).
%%% This function uses pixels larger than the threshold (thr).
%%% 
%%% [numdim,numseed] = size(pos); 
%%% 
%%% Yu Toyoshima, 2013/12/05


% sort image pixel
[imsort,idxsort] = sort(im(:),'descend');
idxuse = idxsort(imsort>thr);


% calc differential index
% assuming medium neighbor (8-neighbor in 2D and 26-neighbor in 3D)
numelblock = cumprod([1,size(im)]);
numdim = ndims(im);
diffarr = [-1,0,1]';
for p=2:numdim
    tmparr = zeros([ones(1,p-1),3]);
    tmparr(1) = -numelblock(p);
    tmparr(3) =  numelblock(p);
    diffarr = bsxfun(@plus,diffarr,tmparr);
end
idxdiff = diffarr(diffarr~=0);


% initialize return matrix
[~,numpos] = size(pos);
idxpos = numelblock(1:numdim)*bsxfun(@minus,pos,1)+1;
[~,idxpossort] = sort(im(idxpos),'ascend');
ret = zeros(size(im));
for p=1:numpos
    ret(idxpos(idxpossort(p))) = p;
end


% flood-fill
numvoxel = numel(im);
numpixeluse = numel(idxuse);
numlocalpeak = 0;
for p=1:numpixeluse
    if ret(idxuse(p))~=0; continue; end;
    
    idxsurround = idxuse(p) + idxdiff;
    retsurround = ret(idxsurround(idxsurround>0 & idxsurround<=numvoxel));
    numpeakassigned = unique(retsurround(retsurround~=0));
    
    switch numel(numpeakassigned)
        case 0 % any surround voxels are not assigned;
            % local peak is found; update peak number and assign
            numlocalpeak = numlocalpeak-1;
            ret(idxuse(p)) = numlocalpeak; % set new peak number
        case 1 % surrond voxels are assigned to the single peak;
            % assign this voxel to the surround voxel
            ret(idxuse(p)) = numpeakassigned;
        otherwise % surround voxels are already assigned to multiple peaks;
            % assign this voxel as a border(=0).
            
            % if the part of the the surround voxels were assigned to one
            % of the local peaks and if another part were assigned to a
            % pre-defined peak, re-assign the former (and the current voxel)
            % to the latter.
            if any(numpeakassigned<0) && sum(numpeakassigned>0)==1
                ret(idxuse(p)) = max(numpeakassigned);
                tmpv = numpeakassigned(numpeakassigned<0);
                for q=1:numel(tmpv)
                    ret(ret==tmpv(q)) = max(numpeakassigned);
                end
            end
    end
end

end
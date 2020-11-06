function [K,H,k1,k2] = curvature(im,idx)
%%% obtain gaussian curvature from 3D image.
%%% 
%%% Input:
%%%  im  ... Grayscale intensities of the image.
%%%  idx ... (If needed) indices of voxels that the curvatures are obatined.
%%% 
%%% Output: cuvatures at idx of input image 
%%%  K   ... Gaussian Curvature.
%%%  H   ... Mean Curvature.
%%%  k1  ... larger Principal Curvature.
%%%  k2  ... smaller Principal Curvature.
%%%
%%% Ref1: http://en.wikipedia.org/wiki/Gaussian_curvature
%%% Ref2: 
%%%  "Computing the Differential Characteristics of Isointensity
%%%  Surfaces", Thirion J., Gourdon A., 1995, Computer Vision and Image
%%%  Understanding.
%%%

numdim = ndims(im);

switch numdim
    case 2 % 2D image
        [pos(1,:),pos(2,:)] = ind2sub(size(im),idx);
        fx = gradient_pos(im,pos,1);
        fy = gradient_pos(im,pos,2);
        fxx = gradient2_pos(im,pos,1,1);
        fxy = gradient2_pos(im,pos,2,1);
        fyy = gradient2_pos(im,pos,2,2);
        
        fx2 = fx.^2;
        fy2 = fy.^2;
        
        K = (2*fx.*fy.*fxy - fx2.*fyy - fy2.*fxx) ./ (fx2+fy2).^(3/2);
        H  = K;
        k1 = K;
        k2 = K;
    
    case 3 % 3D image

        [pos(1,:),pos(2,:),pos(3,:)] = ind2sub(size(im),idx);
        fx = gradient_pos(im,pos,1);
        fy = gradient_pos(im,pos,2);
        fz = gradient_pos(im,pos,3);
        fxx = gradient2_pos(im,pos,1,1);
        fxy = gradient2_pos(im,pos,2,1);
        fxz = gradient2_pos(im,pos,3,1);
        fyy = gradient2_pos(im,pos,2,2);
        fyz = gradient2_pos(im,pos,3,2);
        fzz = gradient2_pos(im,pos,3,3);
        
        
        %%% obtain gaussian curvature using second fundamental form (modified)
        
        fx2 = fx.^2;
        fy2 = fy.^2;
        fz2 = fz.^2;
        
        Lmod = fz.*(fxx.*fz-2.*fx.*fxz)     +    fx2.*fzz;
        Nmod = fz.*(fyy.*fz-2.*fy.*fyz)     +    fy2.*fzz;
        Mmod = fz.*(fxy.*fz-fx.*fyz-fxz.*fy)+ fx.*fy.*fzz;
        
        Amod = (fx2+fy2+fz2); % definition changed from older version
        
        K = (Lmod.*Nmod-Mmod.^2) ./ (fz2.*Amod.^2);
        
        if nargout==1; return; end;
        
        
        %%% obtain mean curvature
        Emod = fx2 + fz2;
        Fmod = fx.*fy;
        Gmod = fy2 + fz2;
        
        H = (Fmod.*Mmod - Emod.*Nmod - Gmod.*Lmod) ./ (2*sqrt(Amod).*Amod.*fz2);
        
        if nargout==2; return; end;
        
        %%% obatain principal curvature
        k1 = H + sqrt(H.^2 - K); % principal curvature (larger one)
        k2 = K./k1;              % principal curvature (smaller one)
        
end

end


function ret = gradient2_pos(im,pos,dim1,dim2)
%%% assuming pos = [numdim,numpos];

% get indices of specified points
siz = size(im);
% numelblock = cumprod([1,siz])';
% idx = (pos'-1)*numelblock(1:end-1)+1;
% d = numelblock(dim1);

% switch to forward/backward diffeerence on edge point
flag_minus = pos(dim1,:)' <=        1 ; % x-1 < 1
flag_plus  = pos(dim1,:)' >= siz(dim1); % x+1 > end
idx_normal = find(~flag_minus & ~flag_plus);

d_plus  = zeros(size(pos));
d_minus = zeros(size(pos));
d_plus (dim1,~flag_plus)  = 1;
d_minus(dim1,~flag_minus) = 1;

% calc derivatives and correct central difference
% ret = im(idx + d*~flag_plus) - im(idx - d*~flag_minus);
ret = gradient_pos(im,pos+d_plus,dim2) - gradient_pos(im,pos-d_minus,dim2);
ret(idx_normal) = ret(idx_normal)*0.5;

end


function ret = gradient_pos(im,pos,dim)
%%% assuming pos = [numdim,numpos];

% get indices of specified points
siz = size(im);
numelblock = cumprod([1,siz])';
idx = (pos'-1)*numelblock(1:end-1)+1;
d = numelblock(dim);

% switch to forward/backward diffeerence on edge point
flag_minus = pos(dim,:)' <=       1 ; % x-1 <= 1
flag_plus  = pos(dim,:)' >= siz(dim); % x+1 >= end
idx_normal = find(~flag_minus & ~flag_plus);

% calc derivatives and correct central difference
ret = im(idx + d*~flag_plus) - im(idx - d*~flag_minus);
ret(idx_normal) = ret(idx_normal)*0.5;

end
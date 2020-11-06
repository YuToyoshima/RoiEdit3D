function peak_detection_14_3(varargin)
%%% peak detection for images of h20p::nls4-mCherry
%%%
%%% Usage:
%%% Please specify path of tiff image(s).
%%% Option file located with program or specified directory will be loaded
%%% and default options will be overwritten.
%%%
%%% Output:
%%%     *_aligned.tif : the image after subpixel alignment
%%%     *_blurred.tif : the image after blurred
%%%     *.mat : data including parameters of fitted elipsoids for nucleus
%%%
%%%
%%% Algorithms and third-party scripts used in this script:
%%%     Miji (FIJI interface for matlab)
%%%     subpixel alignment (dftregistration)
%%%     prefilter (for denoising)
%%%     subtract background
%%%     gaussian blur
%%%     thresholding
%%%     peak detection; get local maxima using image dilation
%%%     grayscale watershed (seeded)
%%%     clump splitting algorithm based on negative curvature
%%%
%%%
%%% This function uses Fiji (another ImageJ distribution).
%%% Please comipile by mcc as following;
%%% mcc('-m','peak_detection_14_3.m','-a',[fileparts(fileparts(which('Miji')))])
%%%
%%% subpixel allignment comes from:
%%% http://www.mathworks.com/matlabcentral/fileexchange/18401-efficient-subpixel-image-registration-by-cross-correlation
%%%
%%% minmaxfilt comes from:
%%% http://www.mathworks.com/matlabcentral/fileexchange/24705-min-max-filter
%%%
%%% marker-controlled watershed comes from:
%%% https://github.com/ijpb/MorphoLibJ/releases
%%%
%%% Written by Yu Toyoshima, 2018/7/18
%%%  originated from peak_detection_140517_13.m
%%%


optDefault = getDefaultOption(); % load default option
[pathImage,optIn] = parseInput(varargin);
numImages = numel(pathImage);

%%% load global option file
pathOptGlobal = [mfilename('fullpath'),'_param.txt'];
[optGlobal,flagOptGlobalFound] = loadOption(pathOptGlobal,optDefault);

%%% main loop
for s=1:numImages
    [imdir,imname,imext] = fileparts(pathImage{s}); %#ok<ASGLU>
    
    %%% load local option file
    pathOptLocal = fullfile(imdir,[mfilename,'_param.txt']);
    [opt,flagOptLocalFound] = loadOption(pathOptLocal,optGlobal);
    if ~isempty(optIn)
        opt = optIn;
    elseif ~flagOptGlobalFound && ~flagOptLocalFound
        disp('Use default option ...');
    end
    
    runMiji(opt.flag_debug); % run Miji
    
    se_peak      = elipsoid_se(opt.radius_dilation);
    se_curvature = elipsoid_se(opt.curvature_distance);
    
    
    %%% read image
    disp(['Read image header: ',pathImage{s}]); ht = tic;
    [~,imp,mfivs] = getImageMiji(pathImage{s},1,1,1);
    myluts = getLuts(imp); %#ok<NASGU>
    fprintf('%g (sec) \n',toc(ht));    
    
    
    %%% subpixel alignment
    disp('subpixel alignment: '); ht = tic;    
    cum2 = calc_alignment(pathImage{s},opt);
    mfivs.setParallelShift(cum2);
    mfivs.setAlign(true);    
    fprintf('subpixel alignment total : %g (sec) \n',toc(ht));        
    
    
    %%% filtering
    fprintf('filtering aligned image. \n'); ht = tic;    
    if opt.flag_reslice_z
        z_ratio = opt.scaling(3)/opt.scaling(1); 
        mfivs.zscale = z_ratio;
        mfivs.setInterpZ(true); % added from v14_3
    end
    impBlurred = setImageMiji(getImageMiji(imp,[],opt.channel_use,opt.frame_use),'blurred');
    ij.IJ.run(impBlurred, opt.prefilter_method, opt.prefilter_option);
    ij.IJ.run(impBlurred,'Subtract Background...',...
        ['rolling=',num2str(opt.radius_background),' stack']);
    ij.IJ.run(impBlurred,opt.blur_method,opt.blur_option);
    imBlurred = double(getImageMiji(impBlurred,[],[]));    
    clear imp mfivs impBlurred;
    fprintf('filterling alignd image total: %g (sec) \n',toc(ht));    
    
    
    %%% thresholding
    fprintf('thresholding (%s) ... ',opt.threshold_method); ht = tic;
    minint = min(imBlurred(:));
    maxint = max(imBlurred(:));
    histogram = int32(hist(imBlurred(:),256));
    hat = ij.process.AutoThresholder();
    threshold_idx = hat.getThreshold(opt.threshold_method,histogram);
    threshold_value = minint + (maxint-minint)/255*(threshold_idx+1);
    clear hat;    
    fprintf('%g (sec) \n',toc(ht));
        
    
    %%% peak detection
    fprintf('peak detection ... '); ht = tic;
    idxOverThr = find(imBlurred>=threshold_value);
    imDilated = minmaxfilt(imBlurred,opt.radius_dilation*2-1,'max','same'); % fast but cube only
    flag_peak = (imDilated(idxOverThr)==imBlurred(idxOverThr));
    pos = zeros(3,sum(flag_peak(:)));
    [pos(1,:),pos(2,:),pos(3,:)] = ind2sub(size(imBlurred),idxOverThr(flag_peak));    
    if ~opt.flag_debug; clear imDilated flag_peak; end    
    fprintf('%g (sec) \n',toc(ht));
    
    
    %%% gray scale watershed
    fprintf('gray scale watershed (ImageJ) ... '); ht = tic;
    impBlurredInv = setImageMiji(-imBlurred,'im_blurred_t0_inv');
    im_peaks = zeros(size(imBlurred));
    for p=1:size(pos,2)
        im_peaks(pos(1,p),pos(2,p),pos(3,p))=p;
    end
    imp_peaks = setImageMiji(im_peaks,'im_peaks');
    im_thr = zeros(size(imBlurred));
    im_thr(idxOverThr) = 1;
    imp_thr = setImageMiji(im_thr,'im_thr');
    hmcwt3d = inra.ijpb.watershed.MarkerControlledWatershedTransform3D(...
        impBlurredInv,imp_peaks,imp_thr,26);
    hmcwt3d.setVerbose(false);
    imp_gwsj = hmcwt3d.applyWithPriorityQueueAndDams();
    im_gws_tmp = getImageMiji(imp_gwsj);    
    clear hmcwt3d;    
    fprintf('%g (sec) \n',toc(ht));
    
    
    %%% remove too-small object
    fprintf('remove too-small object ... '); ht = tic;
    
    hist_gws_tmp = hist(im_gws_tmp(:),1:max(im_gws_tmp(:)));
    idx_large = find(hist_gws_tmp>=opt.min_voxel_remove);
    im_gws = zeros(size(im_gws_tmp));
    idx_gws = cell(1,numel(idx_large));
    for p=1:numel(idx_large)
        idx_gws(p) = {find(im_gws_tmp==idx_large(p))};
        im_gws(idx_gws{p})= p;
    end
    pos = pos(:,idx_large);
    
    im_border = false(size(imBlurred));
    im_border(idxOverThr) = im_gws(idxOverThr)==0;
    numpos = size(pos,2);
    
    if opt.flag_debug
        writePointROI(pos,[imdir,filesep,imname,'_peak_RoiSet.zip']); 
        imp_watershed = setImageMiji(im_gws,'watershed');
        ij.IJ.saveAsTiff(imp_watershed,[imdir,filesep,imname,'_watershed.tif']);
        imp_watershed.close();
        imp_border = setImageMiji(im_border,'border');
        ij.IJ.saveAsTiff(imp_border,[imdir,filesep,imname,'_border.tif']);
        imp_border.close();
    else
        clear idx_voxel;
    end
    
    fprintf('%g (sec) \n',toc(ht));
    
    
    %%% obtain threshold for each sub image
    fprintf('obtain threshold for each sub image ... '); ht = tic;
    
    im_log = -inf(size(imBlurred));
    im_log(idxOverThr) = log10(imBlurred(idxOverThr));
    
    im_localthr = false(size(im_log));
    for p=1:size(pos,2)
        idx_sub = idx_gws{p};
        tmpim_log = im_log(idx_sub);
        tmpim_sc = (tmpim_log-min(tmpim_log)) / (max(tmpim_log)-min(tmpim_log));
        thr = graythresh(tmpim_sc);
        idx_sub2 = idx_sub(tmpim_sc>thr);
        im_localthr(idx_sub2) = true;
    end
    
    if opt.flag_debug
        imp_log = setImageMiji(im_log,'im_log'); 
        ij.IJ.saveAsTiff(imp_log,[imdir,filesep,imname,'_log.tif']);
        imp_log.close();
        imp_localthr = setImageMiji(im_localthr,'im_localthr');
        ij.IJ.saveAsTiff(imp_localthr,[imdir,filesep,imname,'_localthr.tif']);
        imp_localthr.close();
    else
    end
    
    clear im_log idx_overthr;
    fprintf('%g (sec) \n',toc(ht));
    
    
    %%% calc principal curvature of data
    fprintf('calc principal curvature of data ... '); ht = tic;
    
    idx_localthr = find(im_localthr);
    [~,~,~,k2_data] = curvature(imBlurred,idx_localthr);
    im_k2_data = zeros(size(imBlurred));
    im_k2_data(idx_localthr) = k2_data;
    
    %%% remove isolated points
    im_k2_data_gt0 = im_k2_data>0;
    im_k2_data_gt0_holes = imfill(im_k2_data_gt0,'holes') ~= im_k2_data_gt0;
    im_k2_data(im_k2_data_gt0_holes) = 0;
    
    if opt.flag_debug
        imp_k2_data = setImageMiji(im_k2_data<0,'im_k2_data<0'); 
        ij.IJ.saveAsTiff(imp_k2_data,[imdir,filesep,imname,'_k2_data_lt0.tif']);
        imp_k2_data.close();
    else
        clear k2_data idx_localthr im_k2_data_gt0 im_k2_data_gt0_holes;
    end    
    
    fprintf('%g (sec) \n',toc(ht));
    
    
    %%% detecting under-segmentation using principal curvature
    fprintf('detecting under-segmentation using principal curvature ... ');
    ht = tic;
    
    im_border_dist = imdilate(im_border,se_curvature);
    im_k2lt0 = im_k2_data<0 & im_border_dist==0;
    flag_underseg = hist(im_gws(im_k2lt0),1:numpos)>opt.curvature_numvoxels;
    idx_underseg = find(flag_underseg);
    
    if opt.flag_debug
        imp_border_dist = setImageMiji(im_border_dist,'im_border_dist'); 
        ij.IJ.saveAsTiff(imp_border_dist,[imdir,filesep,imname,'_border_dist.tif']);
        imp_border_dist.close();
        imp_k2lt0 = setImageMiji(im_k2lt0,'im_k2lt0');
        ij.IJ.saveAsTiff(imp_k2lt0,[imdir,filesep,imname,'_k2lt0.tif']);
        imp_k2lt0.close();
    else
    end
    clear im_border im_border_dist im_k2_data;        
    fprintf('%g (sec) \n',toc(ht));
    
    
    %%% repair under-segmented objects
    fprintf('repair under-segmented objects ... '); ht = tic;
    
    pos2 = [];
    pos2derived = [];
    im_repair = im_gws.*(~im_k2lt0 & im_localthr);
    ret2 = zeros(size(im_gws));
    for p=1:numpos
        if any(p==idx_underseg)
            im_tmp = (im_repair==p);
            [x,y,z] = ind2sub(size(im_tmp),find(im_tmp));
            im_tmp_sub = im_tmp(min(x):max(x),min(y):max(y),min(z):max(z));
            
            im_dist = bwdist(~im_tmp_sub);
            im_dist_noise = im_dist.*reshape(linspace(1,1+1e-3,numel(im_dist)),size(im_dist));
            im_peak2 = (imdilate(im_dist_noise,se_peak)==im_dist_noise & im_dist_noise>0);
            [xp,yp,zp] = ind2sub(size(im_peak2),find(im_peak2));
            
            %%% remove too-near peaks
            if numel(xp)>1
                disttable = squareform(pdist([xp,yp,zp]));
                flag_too_near = false(size(xp));
                for q=1:numel(xp)
                    flag_too_near(q) = any(disttable(q,q+1:end)<=2);
                end
                xp = xp(~flag_too_near);
                yp = yp(~flag_too_near);
                zp = zp(~flag_too_near);
            end
            
            if size(im_dist_noise,3)>1
                rettmp = graywatershed(im_dist_noise,[xp,yp,zp]',1);
            else
                rettmp = graywatershed(im_dist_noise,[xp,yp]',1);
            end
            rettmp(rettmp>0) = rettmp(rettmp>0)+size(pos2,2);
            ret2(min(x):max(x),min(y):max(y),min(z):max(z)) = ...
                ret2(min(x):max(x),min(y):max(y),min(z):max(z)) + rettmp;
            pos2 = cat(2,pos2,[xp(:)'+min(x)-1;yp(:)'+min(y)-1;zp(:)'+min(z)-1]);
            pos2derived = cat(2,pos2derived,p*ones(1,numel(xp)));
        else
            pos2 = cat(2,pos2,pos(:,p));
            ret2(im_repair==p) = size(pos2,2);
            pos2derived = cat(2,pos2derived,p);
        end
        
    end
    
    if opt.flag_debug
        writePointROI(pos2,[imdir,filesep,imname,'_peak2_RoiSet.zip']); 
        imp_repair = setImageMiji(im_repair,'im_repair');
        ij.IJ.saveAsTiff(imp_repair,[imdir,filesep,imname,'_repair.tif']);
        imp_repair.close();
    else
    end
    clear ret im_k2lt0 im_localthr im_repair im_tmp;    
    fprintf('%g (sec) \n',toc(ht));
    
    
    %%% fit elipsoid with EM-like optimization
    fprintf('fit elipsoid with EM-like optimizaiton ... \n');
    ht = tic;
    
%     if opt.flag_reslice_z
%         %%% reset z scaling
%         imBlurred = double(im_blurred_t0_orig); 
%         pos2(1,:) = round((pos2(1,:)-1)/opt.reslice_scale_factor + 1);
%         pos2(2,:) = round((pos2(2,:)-1)/opt.reslice_scale_factor + 1);
%         pos2(3,:) = round((pos2(3,:)-1)/(z_ratio*opt.reslice_scale_factor) + 1);
%     end
    
    numpos2 = size(pos2,2);
    params_em_in = zeros(10,numpos2);
    params_em_in(2:4,:) = pos2;
    params_em_in(  5,:) = opt.sigma_default(1);
    params_em_in(  8,:) = opt.sigma_default(2);
    params_em_in( 10,:) = opt.sigma_default(3);
    
    emopt.thrint = threshold_value;
    emopt.tol = opt.reltol;
    emopt.thrdist = opt.thrdist;
    emopt.m1 = 3;
    emopt.m2 = 7;
    
    params_em = elipsoid_em_13(params_em_in,imBlurred,emopt);
    
    idx_repaired = find(hist(pos2derived,1:max(pos2derived))>1);
    flag_repaired = any(bsxfun(@eq,pos2derived,idx_repaired'));
    if numel(idx_repaired)==1
        flag_repaired = (pos2derived==idx_repaired);
    end
    
    params_em_before_removefp = params_em; %#ok<NASGU>
    flag_repaired_before_removefp = flag_repaired; %#ok<NASGU>
    
    
    %%% remove false positives based on relative distances
    flag_removefp = false(1,numpos2);
    tmpvec = 1:numpos2;
    for q=1:opt.removefp_numloop
        disp([' remove false positives, loop ',num2str(q)]);
        
        %%% find false positives based on minimum values of scaled distance
        scaled_distance = calc_scaled_distance(params_em);
        scaled_distance(find(eye(numpos2))) = inf; %#ok<FNDSB>
        [min_scdist_val,min_scdist_idx] = min(scaled_distance);
        list_fp = find(min_scdist_val<opt.removefp_min_scdist);
        for p=1:numel(list_fp)
            fp_pair = [list_fp(p),min_scdist_idx(list_fp(p))];
            if    min_scdist_idx(fp_pair(2)) == fp_pair(1) ...
                    && all(~flag_removefp(fp_pair))
                [~,minid] = min(params_em(1,fp_pair));
                flag_removefp(fp_pair(minid)) = true;
            end
        end
        
        %%% find false positives based on minimum values of distance
        raw_disance = squareform(pdist(bsxfun(@times,params_em(2:4,:),opt.scaling')'));
        raw_disance(find(eye(numpos2))) = inf; %#ok<FNDSB>
        [min_dist_val,min_dist_idx] = min(raw_disance);
        list_fp2 = find(min_dist_val<opt.removefp_min_dist);
        for p=1:numel(list_fp2)
            fp_pair = [list_fp2(p),min_dist_idx(list_fp2(p))];
            if    min_dist_idx(fp_pair(2)) == fp_pair(1) ...
                    && all(~flag_removefp(fp_pair))
                [~,minid] = min(params_em(1,fp_pair));
                flag_removefp(fp_pair(minid)) = true;
            end
        end
        
        flag_removefp = flag_removefp | any(~isfinite(params_em),1); % remove illegal roi
        
        %%% if there are no false positives, exit the loop
        if all(~flag_removefp)
            disp('false positives not found; exit the loop.');
            break;
        elseif all(flag_removefp)
            disp('all positives seem to be false; please check the condition.');
            break;
        end
        
        %%% remove false positives and refine
        params_em = params_em(:,~flag_removefp);
        flag_repaired = flag_repaired(~flag_removefp);
        tmpvec = tmpvec(~flag_removefp);
        
        numpos2 = size(params_em,2);
        
        flag_removefp = flag_removefp(~flag_removefp);
        params_em = elipsoid_em_13(params_em);
        
    end
    
    if opt.flag_debug
        writePointROI(params_em(2:4,:),[imdir,filesep,imname,'_peak_em_RoiSet.zip']); 
    else
        clear imBlurred;
    end
    
    fprintf('%g (sec) \n',toc(ht));
    
    
    %%%%% evaluation of Out of Bounds
    
    if opt.flag_reslice_z
        rx = round(params_em(2,:)); 
        ry = round(params_em(3,:));
        rz = round((params_em(4,:)-1)*z_ratio)+1;
    else
        rx = round(params_em(2,:));
        ry = round(params_em(3,:));
        rz = round(params_em(4,:));
    end
    flag_oob =   rx<1 | rx>size(ret2,1) | isnan(rx)...
        | ry<1 | ry>size(ret2,2) | isnan(ry)...
        | rz<1 | rz>size(ret2,3) | isnan(rz);
    
    flag_oob(~flag_oob) = (tmpvec(~flag_oob) ~= ...
        ret2(sub2ind(size(ret2),rx(~flag_oob),...
        ry(~flag_oob),...
        rz(~flag_oob)))); %#ok<NASGU>
    
    
    %%% reset z scaling
%     if opt.flag_reslice_z        
%         params_em(4,:)  = (params_em(4,:)-1)/z_ratio + 1; % zc
%         params_em(7,:)  = params_em(6,:)/z_ratio; % S_13
%         params_em(9,:)  = params_em(9,:)/z_ratio; % S_23
%         params_em(10,:) = params_em(10,:)/(z_ratio*z_ratio); %#ok<NASGU> % S_33
%     end
    
    %%% set variables
    flag_virtual_stack = opt.flag_virtual_stack; %#ok<NASGU>
    idx_frame = opt.frame_use(1,ones(1,numpos2)); %#ok<NASGU>
    scaling = opt.scaling; %#ok<NASGU>
    sigma_default = opt.sigma_default;  %#ok<NASGU>
    flagAlign = true; %#ok<NASGU>
    flagInterpZ = opt.flag_reslice_z;  %#ok<NASGU>
    detopt = opt; %#ok<NASGU>
    
    %%% save mat file
    fprintf('save mat file... ');
    tic;
    
    clear('imp*');
    save([imdir,filesep,imname,'.mat'],...
        'params_em','flag_oob','flag_repaired','emopt','myluts','scaling','sigma_default',...
        'idx_frame','imname','cum2','flag_virtual_stack','detopt','imext','flagAlign','flagInterpZ');
    
    if opt.flag_debug
        save([imdir,filesep,imname,'_debug.mat']); 
    end
    
    fprintf('%g (sec) \n',toc);
    
    disp(' ');
    
end

end



function opt = getDefaultOption()

%%% parameters
opt.flag_align_centroid = false; % whether to use centroid method (faster) than cross-corelation method (precise) in subpixel alignment (only for t-frame)
opt.flag_align_subtract = false; % whether to subtract median value from image in subpixel alignment (sensitive to weak signals)
opt.align_filter_radius = 2; % specify filter radius for alignment; if set to 0, the filtering step will be ignored.
opt.align_max_move = [10,10]; % [x,y] pixels of maximum difference between adjuscent z-slices; Detected movement larger than this value will be ignored.
opt.channel_use = 1; % specify channel for detect objects (for multichannel image)
opt.frame_use = 1; % specify time frame for detect objects
opt.prefilter_method = 'Median 3D...';
opt.prefilter_option = 'x=2 y=2 z=2';
opt.radius_background = 50; % radius of rolling-ball for background subtraction
opt.blur_method = 'Gaussian Blur 3D...';
opt.blur_option = 'x=2 y=2 z=2';
opt.threshold_method = 'Mean'; % for CCD/CMOS images, try Huang, Li, Triangle, etc.
opt.radius_dilation = [5,5,5]; % [x,y,z] (pixels) of dilation radius
opt.min_voxel_remove = 64; % to remove too-small objects
opt.curvature_distance = [5,5,5]; % a negative curvature voxel is ignored when the distance from segmentation border is smaller than this [x,y,z] value
opt.curvature_numvoxels = 50; % a segmented area is marked as under-segmented when the number of negative curvature voxels is larger than this value
opt.sigma_default = [10,10,10]; % default value of [x,y,z] (pixels) of sigma of elipsoid (multivariate gaussian distribution)
opt.thrdist = 4; % specify option parameters for elipsoid fitting
opt.reltol = 5e-8; % specify option parameters for elipsoid fitting
opt.scaling = [0.240,0.240,0.252]; % [x,y,z] um/pixel.
opt.removefp_numloop = 0; % maximum number of loops for remove false positives.
opt.removefp_min_dist = 1.5; % [um] if minimum distance is smaller than this threshold, either roi will be removed
opt.removefp_min_scdist = 1; % if minimum relative distance is smaller than this threshold, either roi will be removed
opt.flag_virtual_stack = true; %%% DEPRECATED; default is TRUE. % if true, memory requirement will be reduced by loading image data through Virtual Stack of ImageJ.
opt.flag_reslice_z = false; % if true, z-slices are interporated (from peak detection to repairing under-segmentation) so that the scaling of z and x is same.
% opt.reslice_scale_factor = 2; %%% DEPRECATED % rescale; valid if flag_reslice_z is true.
% opt.reslice_prefilter_method = 'Median 3D...'; %%% DEPRECATED % valid if flag_reslice_z is true.
% opt.reslice_prefilter_option = 'x=2 y=2 z=2'; %%% DEPRECATED % valid if flag_reslice_z is true.
% opt.reslice_blur_method = 'Gaussian Blur 3D...'; %%% DEPRECATED % valid if flag_reslice_z is true.
% opt.reslice_blur_option = 'x=2 y=2 z=2'; %%% DEPRECATED % valid if flag_reslice_z is true.

% debug mode
% flag_debug = true;
opt.flag_debug = false;

end


function [opt,flagFound] = loadOption(pathOpt,opt)
%%% load option file
flagFound = false;
disp(['Search option file:',pathOpt]);
if numel(dir(pathOpt))==1 % an option file found
    flagFound = true;
    disp('Load option file...')
    fid = fopen(pathOpt,'rt');
    tline = fgetl(fid);
    while ischar(tline)
        try
            eval(['opt.',tline]);
        catch
        end
        tline = fgetl(fid);
    end
    fclose(fid);
else
    disp('Option file not found');
end
end


function runMiji(flag)
%%% run MIJ
if ~exist('MIJ','class') || numel(ij.IJ.getInstance())==0
    strDir = pwd();
    Miji(flag);
    cd(strDir);
end
end


function [pathImage,opt] = parseInput(input)

%%% parse input and get path of image files
% state = warning('off','MATLAB:MKDIR:DirectoryExists'); %#ok<NASGU>
opt = [];
if numel(input)==0
    [imgnames,imgdir] = uigetfile('*.tif;*.mif;*.pif;*.dcv,*.sif','MultiSelect','on');
    if ~iscell(imgnames) % imgnames is string
        if imgnames==0
            return;
        end
        imgnames = {imgnames};
    end
    pathImage = cell(numel(imgnames),1);
    for p=1:numel(imgnames)
        pathImage{p} = fullfile(imgdir,imgnames{p});
    end
else
    flagStruct = cellfun(@isstruct,input);
    if any(flagStruct)
        opt = input{flagStruct};
    end
    pathImage = input(~flagStruct);
end
end


function se = elipsoid_se(radius)
%%% calc structure element for dilation
numdim = numel(radius);
dsum = 0;
for p=1:numdim
    r = ceil(radius(p));
    d = ((-r:r)/radius(p)).^2;
    pidx = 1:numdim;
    pidx(p) = 1;
    pidx(1) = p;
    dsum = bsxfun(@plus,dsum,permute(d(:),pidx));
end
se = dsum<=1;

end


function myluts = getLuts(imp)
%%% save luts for display
luts = imp.getLuts();
if isempty(luts) % single channel image
    luts = imp.getProcessor().getLut();
end
for p=1:numel(luts)
    tmplut = luts(p).getBytes();
    myluts(p).r = typecast(tmplut(  1:256),'uint8');  %#ok<AGROW>
    myluts(p).g = typecast(tmplut(257:512),'uint8');  %#ok<AGROW>
    myluts(p).b = typecast(tmplut(513:768),'uint8');  %#ok<AGROW>
    myluts(p).min = luts(p).min;  %#ok<AGROW>
    myluts(p).max = luts(p).max;  %#ok<AGROW>
end
end


function [cum2,imp,mfivs] = calc_alignment(pathImage,opt)
mini_scale = 4;
flag_use_mini_z = false; % if true, fast but less precise alignment for z

[~,imp,mfivs] = getImageMiji(pathImage,1,1,1);
width  = imp.getWidth();
height = imp.getHeight();
depth  = imp.getNSlices();
numframe = imp.getNFrames();

xm = 1:(width/mini_scale/2);
ym = 1:(height/mini_scale/2);
numxm = numel(xm);
numym = numel(ym);
cmat = false(width,height);
cmat(xm            ,[ym,ym-numym+height]) = true;
cmat(xm-numxm+width,[ym,ym-numym+height]) = true;


% calc alignment with z-slice
output = zeros(4,depth-1,numframe);
filtermatrix = fspecial('average');
fiji_directory = fileparts(fileparts(which('Miji')));
if depth>2
    % align with z-slice
    % get image data with padding
    pad_width  = round( width*0.1);
    pad_height = round(height*0.1);
    % im_pad = zeros(width+pad_width*2,height+pad_height*2);
    range_width  =  (1:width) + pad_width;
    range_height = (1:height) + pad_height;
    
    % fast fourier transform and subpixel registration
    fprintf('\n');
    parfor r=1:numframe
        %     for r=1:numframe
        ht = tic;
        fprintf('calc movement in frame # %d ... ',r);
        
        % for noisy image
        Miji_mod(false,fiji_directory); % setup Miji for parfor environment
        [~,himp2,mfivs2] = getImageMiji(pathImage,1,1,1);
        mfivs2.setMappingSize(100); % set MappedByteBuffer size as 100 image
        mfivs2.clearAndResizeCache(100); % Cache 100 image
        
        if opt.align_filter_radius~=0 %#ok<PFBNS>
            mfivs2.setMedian(true);
            mfivs2.setBlur(true);
        end
        
        tmpim = double(getImageMiji(himp2,1,opt.channel_use,r));
        im_pad = zeros(width+pad_width*2,height+pad_height*2);
        im_pad(range_width,range_height) ...
            = tmpim - opt.flag_align_subtract*mode(tmpim(:));
        if opt.flag_align_subtract
            im_fft_old = fft2(imfilter(im_pad,filtermatrix));
        else
            im_fft_old = fft2(im_pad);
        end
        if flag_use_mini_z
            im_fft_old = reshape(im_fft_old(cmat),[numxm*2,numym*2]); %#ok<UNRCH>
        end
        for p=1:depth-1
            tmpim = double(getImageMiji(himp2,p+1,opt.channel_use,r));
            im_pad(range_width,range_height) ...
                = tmpim - opt.flag_align_subtract*mode(tmpim(:));
            if opt.flag_align_subtract
                im_fft_new = fft2(imfilter(im_pad,filtermatrix));
            else
                im_fft_new = fft2(im_pad);
            end
            if flag_use_mini_z
                im_fft_new = reshape(im_fft_new(cmat),[numxm*2,numym*2]); %#ok<UNRCH>
            end
            [output(:,p,r)] = dftregistration(im_fft_old,im_fft_new,100);
            im_fft_old = im_fft_new;
        end
        fprintf('%g (sec) \n',toc(ht));
    end
    clear('im_pad','im_fft_old','im_fft_new');
end

% alignment
if flag_use_mini_z
    output(3:4,:,:) = output(3:4,:,:)*mini_scale; %#ok<UNRCH>
end

% for noisy image
flag_max_move = any(bsxfun(@ge,abs(output(3:4,:,:)),opt.align_max_move(:)));
output(flag_max_move(ones(1,4),:,:)) = 0;

% temporary aligned in z but currently not for t
cum = cat(2,zeros(2,1,numframe),cumsum(output(3:4,:,:),2));
mfivs.setParallelShift(cum);
mfivs.setAlign(true);
mmpfivs = MyMaxProjFileInfoVirtualStack(mfivs);
impProj = mmpfivs.getImage();

if opt.flag_align_centroid % image alignment by centroid
%     im = getImageMiji(imp,[],opt.channel_use,1);
%     im_proj = double(max(im,[],3));
    im_proj = double(getImageMiji(impProj,1,opt.channel_use,1));
    sum_proj = sum(im_proj(:));
    [width_proj,height_proj] = size(im_proj);
    gx = (1:width_proj);
    gy = (1:height_proj)';
    centroid = zeros(2,numframe);
    centroid(1,1) = sum(gx*im_proj) / sum_proj;
    centroid(2,1) = sum(im_proj*gy) / sum_proj;
    for r=2:numframe
        ht = tic;
        fprintf('calc movement of frame # %d from former frame ... ',r);
%         im = getImageMiji(imp,[],opt.channel_use,1);
%         im_proj = double(max(im,[],3));
        im_proj = double(getImageMiji(impProj,1,opt.channel_use,r));
        sum_proj = sum(im_proj(:));
        centroid(1,r) = sum(gx*im_proj) / sum_proj;
        centroid(2,r) = sum(im_proj*gy) / sum_proj;
        fprintf('%g (sec) \n',toc(ht));
    end
    output2([3,4],:,:) = -diff(centroid,[],2);
    cumshift = bsxfun(@plus,cum,cat(3,zeros(2,1,1),cumsum(output2(3:4,:,:),3)));
    
else % image alignment by minimum spanning tree
    
    %%% make mini fft image
    disp('FFT and make low resolution images...');
    tic;
    im_fft_mini = zeros(numxm*2*numym*2,numframe);
    for ct=1:numframe
%         im = getImageMiji(imp,[],opt.channel_use,ct);
%         im_proj = double(max(im,[],3)); % assume z-stack is already aligned
        im_proj = double(getImageMiji(impProj,1,opt.channel_use,ct));
        im_fft_tmp = fft2(im_proj);
        im_fft_mini(:,ct) = im_fft_tmp(cmat);
    end
    im_fft_mini = reshape(im_fft_mini,[numxm*2,numym*2,numframe]);
    toc;
    
    %%% dft registration by mini images
    disp('Obtain image similarity by dftregistration using low resolution images...')
    tic;
    dftres_mini = zeros(numframe,numframe,4);
    parfor ct1 = 2:numframe
        % for ct1 = 2:numframe
        ht = tic;
        dftres_tmp = zeros(numframe,4);
        for ct2 = 1:ct1-1
            dftres_tmp(ct2,:) = ...
                dftregistration(im_fft_mini(:,:,ct1),im_fft_mini(:,:,ct2)); %#ok<PFBNS>
        end
        dftres_mini(ct1,:,:) = dftres_tmp;
        disp(['Asymmetric loop of dftregistration for frame #',num2str(ct1),...
            ' ...  ',num2str(toc(ht)),' (sec)']);
    end
    disttable = dftres_mini(:,:,1) + dftres_mini(:,:,1)' + diag(inf*ones(1,numframe));
    toc;
    
    %%% calc minimum spanning tree (MST)
    %%% MST minimizes sum of edge weights when the tree covers all the nodes
    disp('Building minimum spanning tree based on similarity...');
    tic;
    g = graph(disttable);
    g.Nodes.Name = cellfun(@num2str,num2cell(1:numframe)','UniformOutput',false);
    t = minspantree(g);
    idx_ecc = 1;
    edgetonew = t.bfsearch(idx_ecc,'edgetonew'); % edges in table format
    numedges = size(edgetonew,1);
    % figure;
    % plot(t,'NodeLabel',t.Nodes.Name,'Layout','layered','Source',num2str(idx_ecc));
    toc;
    
    %%% determine shift amount by 2nd dft registration
    disp('Determine shift amount by 2nd dftregistration using high resolution images...');
    tic;
    dftres = zeros(numedges,4);
    parfor c=1:numedges
        dftres(c,:) = ...
            dftfun(pathImage,opt.channel_use,edgetonew(c,:),cum,fiji_directory); %#ok<PFBNS>
    end
    cumshift = zeros(numframe,2);
    for c=1:numedges
        cumshift(edgetonew(c,2),:) = cumshift(edgetonew(c,1),:)+dftres(c,3:4);
    end
    toc;
end

cum2 = bsxfun(@plus,permute(cumshift,[2,3,1]),cum);

end


function dftres = dftfun(pathImage,channel_use,tvec,cum,fiji_directory)
disp(['treefit (fitfunc): ',num2str(tvec(1)),' to ',num2str(tvec(2))]);
Miji_mod(false,fiji_directory); % setup Miji for parfor environment

% [~,himp_orig,mfivs_orig] = getImageMiji(pathImage,[],[],1);
[~,~,mfivs] = getImageMiji(pathImage,[],[],1);
mfivs.setParallelShift(cum);
mfivs.setAlign(true);
mfivs.setMappingSize(100); % set MappedByteBuffer size as 100 image
mfivs.clearAndResizeCache(100); % Cache 100 image
mmpfivs = MyMaxProjFileInfoVirtualStack(mfivs);
impProj = mmpfivs.getImage();
% im = getImageMiji(himp_orig,[],channel_use,tvec);
% im_proj = double(max(im,[],3)); % assume aligned
im_proj = double(getImageMiji(impProj,1,channel_use,tvec));

im_fft = fft2(im_proj);
dftres = dftregistration(im_fft(:,:,:,:,1),im_fft(:,:,:,:,2),100);
end




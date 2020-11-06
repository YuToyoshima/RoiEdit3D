function [im,imp,mmfivs] = getImageMiji(imp,zvec,cvec,tvec)
%%% get data from ImageJ(Miji)'s image
%%% 
%%% This function can be used like as MIJ.getImage()
%%%

%%% run MIJ
% if ~exist('MIJ','class') || numel(ij.IJ.getInstance())==0
if ~exist('MIJ','class')
    strDir = pwd();
    Miji(false);
    cd(strDir);
end


%%% when filepath specified, open the image using MIJ
if ischar(imp) % filepath specified
    if exist('zvec','var') || exist('cvec','var') || exist('tvec','var')
        %%% if zvec | cvec | tvec is given, open the image with virtual stack        
        
        %%% read by self class; cannot handle formal BigTiff.
        %%% ImageJ saves large Tiff File (>4GB) as own way.
        %%% Inherent classes (and inherited classes) can handle it.
        % imp = MyFileInfoVirtualStack(imp,false).getImage();
        mmfivs = MyMultiFileInfoVirtualStack(imp,false);
        imp = mmfivs.getImage(); 

        %%% read by bio-formats plugin; can handle BigTiff but too slow.
%         loci.common.DebugTools.enableLogging('OFF');
%         hio = loci.plugins.in.ImporterOptions();
%         hio.setId(imp);
%         hio.setVirtual(true);
%         imps = loci.plugins.BF.openImagePlus(hio);
%         imp = imps(1);
    else
        imp = ij.IJ.openImage(imp);
    end
end

width  = imp.getWidth();
height = imp.getHeight();
depth  = imp.getNSlices();
bitdepth = imp.getBitDepth();
numcolor = imp.getNChannels();
numframe = imp.getNFrames();
hstack = imp.getStack();

if ~exist('zvec','var') || isempty(zvec); zvec = 1:depth; end;
if ~exist('cvec','var') || isempty(cvec); cvec = 1:numcolor; end;
if ~exist('tvec','var') || isempty(tvec); tvec = 1:numframe; end;


switch bitdepth
    case 8
        strType = 'uint8';
    case 16
        strType = 'uint16';
    case 32
        strType = 'single';
end

% get image data
% im = zeros(width*height,depth,numcolor,numframe,strType);
% for q=1:numcolor
%     for p=1:depth
%         for r=1:numframe
%             idxstack = imp.getStackIndex(q,p,r);
%             im(:,p,q,r) = typecast(hstack.getPixels(idxstack),strType);
%         end
%     end
% end
% 
% im = reshape(im,[width,height,depth,numcolor,numframe]);

% get image data
im = zeros(width*height,numel(zvec),numel(cvec),numel(tvec),strType);
for q=1:numel(cvec)
    for p=1:numel(zvec)
        for r=1:numel(tvec)
            idxstack = imp.getStackIndex(cvec(q),zvec(p),tvec(r));
            im(:,p,q,r) = typecast(hstack.getPixels(idxstack),strType);
        end
    end
end

im = reshape(im,[width,height,numel(zvec),numel(cvec),numel(tvec)]);



end
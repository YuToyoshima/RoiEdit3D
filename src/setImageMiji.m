function imp = setImageMiji(im,imp)
%%% set data to ImageJ(Miji)'s image
%%% 
%%% This function can be used like as MIJ.createImage()
%%%

%%% run MIJ
% if ~exist('MIJ','class') || numel(ij.IJ.getInstance())==0
if ~exist('MIJ','class')
    strDir = pwd();
    Miji(false);
    cd(strDir);
end

%%% when title specified, create the image
if ischar(imp) % title specified
    [width,height,depth,numcolor,numframe] = size(im);
    strType = class(im);
    
    switch strType
        case 'logical'
            strType = 'uint8';
            bitdepth = 8;
        case 'uint8'
            bitdepth = 8;
        case 'uint16'
            bitdepth = 16;
        case 'single'
            bitdepth = 32;
        case 'double' % loss of precision
            strType = 'single';
            bitdepth = 32; 
    end
    
    imp = ij.IJ.createHyperStack(imp,width,height,numcolor,depth,numframe,bitdepth);
    
else
    width  = imp.getWidth();
    height = imp.getHeight();
    depth  = imp.getNSlices();
    numcolor = imp.getNChannels();
    numframe = imp.getNFrames();
    bitdepth = imp.getBitDepth();
    
    switch bitdepth
        case 8
            strType = 'uint8';
        case 16
            strType = 'uint16';
        case 32
            strType = 'single';
    end
    
end


hstack = imp.getStack();

% set image data
im = reshape(im,[width*height,depth,numcolor,numframe]);
for q=1:numcolor
    for p=1:depth
        for r=1:numframe
            idxstack = imp.getStackIndex(q,p,r);
            hstack.setPixels(cast(im(:,p,q,r),strType),idxstack);
        end
    end
end
imp.setStack(hstack);

end
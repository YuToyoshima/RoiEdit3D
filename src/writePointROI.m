function writePointROI(pos,strpath)
%%% write ROI data as ImageJ's RoiSet.zip
%%%
%%%


%%% run MIJ
if ~exist('MIJ','class') || numel(ij.IJ.getInstance())==0
    strDir = pwd();
    Miji(false);
    cd(strDir);
end

hrm = ij.plugin.frame.RoiManager.getInstance();
if isempty(hrm)
    hrm = ij.plugin.frame.RoiManager();
end
hrm.setVisible(false);
if hrm.getCount()~=0 % RoiManager already contains some Rois.
    hrm.runCommand('Deselect');
    hrm.runCommand('Delete');
end

numpoint = size(pos,2);
for p=1:numpoint
    tmproi = ij.gui.PointRoi(pos(1,p),pos(2,p));
    tmproi.setPosition(pos(3,p));
    hrm.addRoi(tmproi);
    hrm.select(p-1);
    hrm.runCommand('Rename',sprintf('%04d-%04d-%04d',round(pos(3,p)),...
                                                     round(pos(2,p)),...
                                                     round(pos(1,p))));
end
hrm.runCommand('Deselect');
hrm.runCommand('Save',strpath);
hrm.runCommand('Delete');

end
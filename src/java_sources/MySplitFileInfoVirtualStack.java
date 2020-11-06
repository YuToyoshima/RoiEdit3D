import ij.IJ;

/***
 * This plugin handles multi-page TIFF files under MyMultiFileInforVirtualStack.
 ***/
public class MySplitFileInfoVirtualStack extends MyFileInfoVirtualStack {
    
    public String[] strPaths;
    public MyFileInfoVirtualStack mfivsIn;
    public int nFile = 1;
    public int nImagesInFile = 0;
    public int idxFile = 0;
    
    
    public MySplitFileInfoVirtualStack(String[] paths, boolean show, int nz, int nt) {
//     public MySplitFileInfoVirtualStack(String[] paths, booleasn show) {
        super(paths[0],false); // call super-class constructor
        super.setMappingSize(1); // reduce memory consumpution, may not required
        strPaths = paths;
        nFile = strPaths.length;
        openNewMfivs(0);
        nImagesInFile = mfivsIn.getSize();
        myimp = mfivsIn.getImage();
        nxorig = myimp.getWidth();
        nyorig = myimp.getHeight();
        ncorig = 1;
        nzorig = nz;
        ntorig = nt;
        nImages = nzorig * ntorig;
    }
    
    
    public void openNewMfivs(int idx) {
        idxFile = idx;
        mfivsIn = new MyFileInfoVirtualStack(strPaths[idxFile],false);
        mfivsIn.clearAndResizeCache(1);
    }
    
    
    
    
    
    @Override
    public Object getSubImage(int slice) {
//         System.out.println("getSubImage, nImagesInFile="+nImagesInFile);
        if (nImagesInFile == 0) { return super.getSubImage(slice); } // in init
        int idxFileTarget  = (slice-1) / nImagesInFile; // file index is 0-based
        int idxSliceTarget = (slice-1) % nImagesInFile + 1; // slice index is 1-based
        if (flagDebug) { System.out.println("getSubImage, slice="+slice+",idxFileTarget="+idxFileTarget+", idxSliceTarget="+idxSliceTarget); }
        if (idxFileTarget!=idxFile) { openNewMfivs(idxFileTarget); }
        return mfivsIn.imc.getUnchecked((double)idxSliceTarget);
    }
    
    
    @Override
    public Object getSubImage(int slice, int[] range) {
        if (nImagesInFile == 0) { return super.getSubImage(slice); } // in init
        int idxFileTarget  = (slice-1) / nImagesInFile; // file index is 0-based
        int idxSliceTarget = (slice-1) % nImagesInFile + 1; // slice index is 1-based
        if (idxFileTarget!=idxFile) { openNewMfivs(idxFileTarget); }
        return mfivsIn.getSubImage(idxSliceTarget, range);
    }
       
       
    @Override
    public void setMappingSize(int msize) {
        mfivsIn.setMappingSize(msize);
    }
    
    @Override
    public void resetMappingSize() {
        mfivsIn.resetMappingSize();
    }
    
    @Override
    public String getPath() { // for write mif file
        StringBuffer buf = new StringBuffer();
        buf.append(strPaths[0]);
        for (int i=1; i<strPaths.length; i++) {
            buf.append(",");
            buf.append(strPaths[i]);
        }
        return buf.toString();
    }
    
    @Override
    public void deleteSlice(int n) {
        IJ.error("Virtual Stack","Cannot delete slice !");
        return;
    }
    
    
}

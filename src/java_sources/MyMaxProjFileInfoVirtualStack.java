import ij.ImagePlus;
import ij.VirtualStack;
import ij.CompositeImage;
import ij.io.FileInfo;
import ij.process.ImageProcessor;
import ij.process.ByteProcessor;
import ij.process.ShortProcessor;
import ij.process.FloatProcessor;
import java.io.File;
import com.google.common.cache.*;

public class MyMaxProjFileInfoVirtualStack extends VirtualStack {
    
    ImagePlus myimp;
    public MyFileInfoVirtualStack mfivs;
    public MyMultiFileInfoVirtualStack mmfivs;
    public boolean flagUseMMFIVS = false;
    public boolean flagDebug = false;
    public int nxorig = 1;
    public int nyorig = 1;
    public int nzorig = 1;
    public int ncorig = 1;
    public int ntorig = 1;
    public int nImages;
    public FileInfo fi;
    
    /* On-memory cache for image */
    public LoadingCache<Integer,Object> imc = CacheBuilder.newBuilder()
    .maximumSize(1000)
    .build( new CacheLoader<Integer, Object>() {
        public Object load(Integer key) { return calcImage(key); }
    });
    
    public MyMaxProjFileInfoVirtualStack(MyFileInfoVirtualStack in) {
        flagUseMMFIVS = false;
        mfivs = in;
        init();
    }
    
    public MyMaxProjFileInfoVirtualStack(MyMultiFileInfoVirtualStack in) {
        flagUseMMFIVS = true;
        mmfivs = in;
        init();
    }
    
    public void init() {
        int[] tmpar = getOriginalDimensionsInner();
        nxorig = tmpar[0];
        nyorig = tmpar[1];
        nzorig = 1;
        ncorig = tmpar[3];
        ntorig = tmpar[4];
        nImages = ncorig * ntorig;
        
        File file = new File(getPathInner());
        ImagePlus tmpimp = getImageInner();
        
        fi = new FileInfo();
        fi.width = nxorig;
        fi.height = nyorig;
        fi.offset = 0;
        fi.fileType = tmpimp.getFileInfo().fileType;
        fi.fileName = "MaxProj_"+file.getName();
        fi.directory = file.getParent()+file.separator;
        fi.nImages = 1;        
        
        myimp = new ImagePlus(fi.fileName, this);
        myimp.setFileInfo(fi);
        myimp.setDimensions(ncorig,nzorig,ntorig);
        if (tmpimp instanceof CompositeImage) {
            myimp = new CompositeImage(myimp, CompositeImage.COMPOSITE);
            ((CompositeImage)(myimp)).copyLuts(tmpimp);
        }
    }
    
    public ImagePlus getImage() { return myimp; }
    
    public ImagePlus getImageInner() {
        if (flagUseMMFIVS) { return mmfivs.getImage(); }
        else { return mfivs.getImage(); }
    }
    
    public int[] getOriginalDimensions() {
        return new int[] {nxorig,nyorig,nzorig,ncorig,ntorig};
    }
    
    public int[] getOriginalDimensionsCZT() {
        return new int[] {nxorig,nyorig,ncorig,nzorig,ntorig};
    }
    
    public int[] getOriginalDimensionsInner() {
        if (flagUseMMFIVS) { return mmfivs.getOriginalDimensions(); }
        else { return mfivs.getOriginalDimensions(); }
    }
    
    public String getPathInner() {
        if (flagUseMMFIVS) { return mmfivs.getPath(); }
        else { return mfivs.getPath(); }
    }
    
    /** Returns an ImageProcessor for the specified image,
     * were 1<=n<=nImages. Returns null if the stack is empty.
     */
    public ImageProcessor getProcessor(int n) {
        if (n<1 || n>nImages) {
            throw new IllegalArgumentException("Argument out of range: "+n);
        }
        if (flagDebug) { System.out.println("getProcessor, n="+n); }
//         MyFileOpener mfo = new MyFileOpener(fi);
//         ImagePlus imp = mfo.open(false,imc.getUnchecked(n+rotation));
//         return imp.getProcessor();
        int[] dims = getImageInner().getDimensions();
        int w  = dims[0];
        int h  = dims[1];
        switch (getBitDepth()) {
            case 8:
                return new ByteProcessor(w,h,(byte[])imc.getUnchecked(n));
            case 16:
                ImageProcessor ip =new ShortProcessor(w,h);
                ip.setPixels(imc.getUnchecked(n));
                return ip;
            case 32:
            default:
                return new FloatProcessor(w,h,(float[])imc.getUnchecked(n));
        }
        
    }
    
    public Object calcImage(int key) {
        int cc = (key-1)%ncorig+1; // 1-based
        int ct = (key-cc)/ncorig+1; // 1-based
        
        Object ret;        
        int[] dims = getImageInner().getDimensions();
        int w  = dims[0];
        int h  = dims[1];
        int nz = dims[3];
        switch (getBitDepth()) {
            case 8: ret = new byte[w*h]; break;
            case 16: ret = new short[w*h]; break;
            case 32: default: ret = new float[w*h]; break;
        }
        for (int cz=1; cz<=nz; cz++) {
            int que = ((ct-1)*nz+(cz-1))*ncorig+cc;
            if (flagDebug) { System.out.println("cc="+cc+", cz="+cz+", ct="+ct+", que="+que); }
            ret = this.max(ret,getProcessorInner(que));
        }
        return ret;
    }
    
    public ImageProcessor getProcessorInner(int key) {
        if (flagUseMMFIVS) { return mmfivs.getProcessor(key);}
        else { return mfivs.getProcessor(key);}
    }
    
    public Object max(Object ar, ImageProcessor ip) {
        switch (ip.getBitDepth()) {
            case 8: return max((byte[])(ar),(byte[])(ip.getPixels()));
            case 16: return max((short[])(ar),(short[])(ip.getPixels()));
            case 32: default: return max((float[])(ar),(float[])(ip.getPixels()));
        }
    }
    
    public byte[] max(byte[] ar1, byte[] ar2) {
        for (int i=0; i<ar1.length; i++) { ar1[i] = (byte)(Math.max(ar1[i],ar2[i])); }
        return ar1;
    }
    
    public short[] max(short[] ar1, short[] ar2) {
        for (int i=0; i<ar1.length; i++) { ar1[i] = (short)(Math.max(ar1[i],ar2[i])); }
        return ar1;
    }
    
    public float[] max(float[] ar1, float[] ar2) {
        for (int i=0; i<ar1.length; i++) { ar1[i] = Math.max(ar1[i],ar2[i]); }
        return ar1;
    }
    
    public int getBitDepth() {
        if (flagUseMMFIVS) { return mmfivs.getBitDepth(); }
        else { return mfivs.getBitDepth(); }
    }
    
    public int getSize() { return nImages; }
    public String getSliceLabel(int n) { return null; }
    
    public void clearAndResizeCache(long maxnum){
        imc = CacheBuilder.newBuilder()
        .maximumSize(maxnum)
        .build( new CacheLoader<Integer, Object>() {
            public Object load(Integer key) { return calcImage(key); }
        });
    }
    
}
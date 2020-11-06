import ij.*;
import ij.plugin.filter.BackgroundSubtracter;
import ij.process.*;
import ij.gui.*;
import ij.io.*;
import ij.measure.Calibration;
import java.awt.*;
import java.io.*;
import java.nio.ByteBuffer;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.util.Properties;
import java.util.ArrayList;
import java.util.List;
import com.google.common.cache.*;

/** This plugin opens a multi-page TIFF file as a virtual stack.
 *  It implements the File/Import/TIFF Virtual Stack command. */
public class MyFileInfoVirtualStack extends VirtualStack {
    public FileInfo[] info;
    public int nImages;
    public RandomAccessFile raf;
    // public List<ByteBuffer> mappings = new ArrayList<ByteBuffer>();
    public MappedByteBuffer mbb;
    public FileChannel fch;
    public long fileSize, offsetFileInfo;
    public int idx_mbb;
    public static final int MAPPING_SIZE_MAX = 1 << 30; // MappedByteBuffer cannot handle over 2GB
    public long mapping_size;
    public long sliceBytes;
    public long sliceBytesWithGap;
    ImagePlus myimp;
    
    public String strPath;
    public int interpolationMethod = ij.process.ImageProcessor.NONE;
    public boolean flagAlign = false;
    public boolean flagFlipX = false;
    public boolean flagFlipY = false;
    public boolean flagFlipZ = false;
    public boolean flagLn = false;
    public boolean flagInterpZ = false;
    public boolean flagSubBG = false;
    public boolean flagBlur = false;
    public boolean flagMedian = false;
    public boolean flagDebug = false;
    public double[][][] parallelShift = null;
    public double rotation = 0; // 0<rotation<0.5
    public double zscale = 1;
    public int nxorig = 1;
    public int nyorig = 1;
    public int nzorig = 1;
    public int ncorig = 1;
    public int ntorig = 1;
    public double radiusSubBG = 50;
    public BackgroundSubtracter bs = new BackgroundSubtracter();
    
    public VirtualStack vsoir = null;
    public OirDecoder od = null;
    
    /* On-memory cache for image
     * n-th image of r[rad] rotated hyperstack can be obtained by calling imc.get(K).
     * K is a key for cache and encoded by n + 0.5*r', where r' satisfies 0<=r'<=1.
     * Such r' can be obtained by r' = r/(2*pi) for 0<r<2*pi, or
     * mod(mod(r,2*pi)+2*pi,2*pi)/(4*pi) for any r of real value.
     * If r' equals to 0, the cache returns both Aligned and Fliped image.
     * If r' does not equals to 0, the cache returns rotated images of that for r=0.
     */
    public LoadingCache<Double,Object> imc = CacheBuilder.newBuilder()
    .maximumSize(500)
    .build( new CacheLoader<Double, Object>() {
        public Object load(Double key) { return calcImage(key); }
    });
    
    public LoadingCache<Double,ArrayList> grid = CacheBuilder.newBuilder()
    .maximumSize(100)
    .build( new CacheLoader<Double,ArrayList>() {
        public ArrayList load(Double radian) { return getGridRotated(radian); }
    });
    
    /* Constructs a FileInfoVirtualStack from a FileInfo object. */
    public MyFileInfoVirtualStack(FileInfo fi) {
        info = new FileInfo[1];
        info[0] = fi;
        myOpen(true);
    }
    
    /* Constructs a FileInfoVirtualStack from a FileInfo
        object and displays it if 'show' is true. */
    public MyFileInfoVirtualStack(FileInfo fi, boolean show) {
        info = new FileInfo[1];
        info[0] = fi;
        myOpen(show);
    }
    
    public MyFileInfoVirtualStack(FileInfo[] info, boolean show) {
        this.info = info;
        myOpen(show);
    }
    
    public MyFileInfoVirtualStack(String path, boolean show) {
        strPath = path;
        File file = new File(path);
        try {
            raf = new RandomAccessFile(file,"r");
        } catch (FileNotFoundException e) {
            IJ.error("Virtual Stack","File Not Found");
        }
        
        String strDir = file.getParent()+file.separator;
        String strName = file.getName();
        try {
            vsoir = null;
            od = null;
            if(strName.endsWith(".oir")) {
                od = new OirDecoder(strDir,strName);
                myOpenOir(show);
                return;
            } else {
                if (strName.endsWith(".pif")) {
                    PifDecoder pd = new PifDecoder(strDir,strName);
                    info = pd.getPifInfo();
                } else if(strName.endsWith(".sif")) {
                    SifDecoder sd = new SifDecoder(strDir,strName);
                    info = sd.getSifInfo();
                } else if(strName.endsWith(".dcv")) {
                    DcvDecoder dd = new DcvDecoder(strDir,strName);
                    info = dd.getDcvInfo();
                } else {
                    TiffDecoder td = new TiffDecoder(strDir,strName);
                    info = td.getTiffInfo();
                }
            }
        }
        catch (IOException e) {
            String msg = e.getMessage();
            if (msg==null||msg.equals("")) msg = ""+e;
            IJ.error("try-catch in TiffDecoder", msg);
            return;
        }
        if (info==null || info.length==0) {
            IJ.error("Virtual Stack", "This does not appear to be a correct image stack");
            return;
        }
        if (IJ.debugMode)
            IJ.log(info[0].debugInfo);
        myOpen(show);
    }
    
    void myOpen(boolean show) {
        FileInfo fi = info[0];
        int n = fi.nImages;
        sliceBytes = fi.width*fi.height*fi.getBytesPerPixel();
        sliceBytesWithGap = sliceBytes + fi.gapBetweenImages;
        if (info.length==1 && n>1) {
            info = new FileInfo[n];
            for (int i=0; i<n; i++) {
                info[i] = (FileInfo)fi.clone();
                info[i].nImages = 1;
                info[i].longOffset = fi.getOffset() + i*sliceBytesWithGap;
            }
        }
        //nImages = n;
        nImages = info.length;
        
        // cannot handle compressed image
        if (fi.compression!=FileInfo.COMPRESSION_NONE
        && fi.compression!=FileInfo.COMPRESSION_UNKNOWN) {
            System.err.println("Error: Cannot handle compressed image. Exit.");
            return;
        }
        
        /* setup MappedByteBuffer */
        fch = raf.getChannel();
        try {
            fileSize = raf.length();
        } catch (IOException e) {
            IJ.error("Virtual Stack","Error occurred in setup RandomAccessFile");
        }
        mapping_size = (MAPPING_SIZE_MAX/sliceBytesWithGap)*sliceBytesWithGap;
        offsetFileInfo = fi.getOffset();
        updateMappedByteBuffer(0);
        
        FileOpener fo = new FileOpener(info[0]);
        ImagePlus imp = fo.open(false);
        if (nImages==1 && fi.fileType==FileInfo.RGB48) {
            if (show) imp.show();
            myimp = imp;
            return;
        }
        Properties props = fo.decodeDescriptionString(fi);
        
        setBitDepth(imp.getBitDepth()); // for avoiding error occured in the next line
        
        ImagePlus imp2 = new ImagePlus(fi.fileName, this);
        imp2.setFileInfo(fi);
        if (imp!=null && props!=null) {
            //setBitDepth(imp.getBitDepth());
            imp2.setCalibration(imp.getCalibration());
            imp2.setOverlay(imp.getOverlay());
            if (fi.info!=null)
                imp2.setProperty("Info", fi.info);
            int channels = getInt(props,"channels");
            int slices = getInt(props,"slices");
            int frames = getInt(props,"frames");
            if (channels*slices*frames==nImages) {
                imp2.setDimensions(channels, slices, frames);
                if (getBoolean(props, "hyperstack"))
                    imp2.setOpenAsHyperStack(true);
            }
            nxorig = imp2.getWidth();
            nyorig = imp2.getHeight();
            ncorig = channels;
            nzorig = slices;
            ntorig = frames;
            zscale = info[0].pixelDepth / info[0].pixelHeight;
            if (channels>1 && fi.description!=null) {
                //int mode = IJ.COMPOSITE;
                int mode = CompositeImage.COMPOSITE;
                if (fi.description.indexOf("mode=color")!=-1)
                    //mode = IJ.COLOR;
                    mode = CompositeImage.COLOR;
                else if (fi.description.indexOf("mode=gray")!=-1)
                    //mode = IJ.GRAYSCALE;
                    mode = CompositeImage.GRAYSCALE;
                imp2 = new CompositeImage(imp2, mode);
            }
        }
        if (show) imp2.show();
        myimp = imp2;
    }
    
    
    public void myOpenOir(boolean show){
        vsoir = od.getVirtualStack();
        nImages = vsoir.getSize();
        ImagePlus imp = od.getImage();
        FileInfo tmpfi = imp.getFileInfo();
        info = new FileInfo[nImages];
        for (int i=0; i<nImages; i++) {
            info[i] = (FileInfo)tmpfi.clone();
            info[i].nImages = 1;
        }
        setBitDepth(imp.getBitDepth()); // may not be required
        int[] dims = imp.getDimensions();
        nxorig = dims[0];
        nyorig = dims[1];
        ncorig = dims[2];
        nzorig = dims[3];
        ntorig = dims[4];
        zscale = info[0].pixelDepth / info[0].pixelHeight;
        imp.setStack(this);
        ((CompositeImage)imp).setMode(CompositeImage.COMPOSITE);
        if (show) imp.show();
        myimp = imp;
    }
    
    
    public ImagePlus getImage() {
        return myimp;
    }
    
    public String getPath() { return strPath; }
    
    int getInt(Properties props, String key) {
        Double n = getNumber(props, key);
        return n!=null?(int)n.doubleValue():1;
    }
    
    Double getNumber(Properties props, String key) {
        String s = props.getProperty(key);
        if (s!=null) {
            try {
                return Double.valueOf(s);
            } catch (NumberFormatException e) {}
        }
        return null;
    }
    
    boolean getBoolean(Properties props, String key) {
        String s = props.getProperty(key);
        return s!=null&&s.equals("true")?true:false;
    }
    
    /** Deletes the specified image, were 1<=n<=nImages. */
    public void deleteSlice(int n) {
        if (n<1 || n>nImages)
            throw new IllegalArgumentException("Argument out of range: "+n);
        if (nImages<1) return;
        for (int i=n; i<nImages; i++)
            info[i-1] = info[i];
        info[nImages-1] = null;
        nImages--;
    }
    
    /** Returns an ImageProcessor for the specified image,
     * were 1<=n<=nImages. Returns null if the stack is empty.
     */
    public ImageProcessor getProcessor(int n) {
        if (n<1 || n>nImages) {
            throw new IllegalArgumentException("Argument out of range: "+n);
        }
        if (IJ.debugMode) { IJ.log("FileInfoVirtualStack: "+n+", "+info[n-1].getOffset()); }
        //if (n>1) IJ.log("  "+(info[n-1].getOffset()-info[n-2].getOffset()));
        info[n-1].nImages = 1; // why is this needed?
        MyFileOpener mfo = new MyFileOpener(info[n-1]);
        ImagePlus imp = mfo.open(false,imc.getUnchecked(n+rotation));
        return imp.getProcessor();
    }
    
    /* set parallelShift */
    public void setParallelShift(double[][][] in) {
        parallelShift = in;
        imc.invalidateAll();
        grid.invalidateAll();
    }
    
    /* set Flag for alignment */
    public void setAlign(boolean flag) {
        flagAlign = flag;
        imc.invalidateAll();
        grid.invalidateAll();
    }
    
    /* set FlipX */
    public void setFlipX(boolean flag) {
        flagFlipX = flag;
        imc.invalidateAll();
        grid.invalidateAll();
    }
    
    /* set FlipY */
    public void setFlipY(boolean flag) {
        flagFlipY = flag;
        imc.invalidateAll();
        grid.invalidateAll();
    }
    
    /* set FlipX */
    public void setFlipZ(boolean flag) {
        flagFlipZ = flag;
        imc.invalidateAll();
        grid.invalidateAll();
    }
    
    /* set Logarighmic display */
    public void setLn(boolean flag) {
        flagLn = flag;
        imc.invalidateAll();
        grid.invalidateAll();
    }
    
    /* set interporation for Z */
    public void setInterpZ(boolean flag) {
        // preparing
        flagInterpZ = flag;
        int[] dims = myimp.getDimensions();
        int nznew = nzorig;
        if (flagInterpZ) { nznew = ((int)((nzorig-1)*zscale)) + 1; }
        
        // update nImages and info
        nImages = dims[2]*nznew*dims[4];
        FileInfo[] tmpinfo = new FileInfo[nImages];
        for (int i=0; i<nImages; i++) {
            if (i<info.length) { tmpinfo[i] = info[i]; }
            else               { tmpinfo[i] = info[i%dims[2]]; }
        }
        info = tmpinfo;
        myimp.setDimensions(dims[2],nznew,dims[4]);
        
        // update caliburation information
        Calibration cal = myimp.getCalibration();
        if (flagInterpZ) { cal.pixelDepth = cal.pixelHeight; }
        else             { cal.pixelDepth = cal.pixelHeight*zscale; }
        myimp.setCalibration(cal);
        
        // invalidate image cache
        imc.invalidateAll();
        grid.invalidateAll();
    }
    
    /* set radius for background subtraction */
    public void setRadiusSubBG(double r) {
        radiusSubBG = r;
        imc.invalidateAll();
        grid.invalidateAll();
    }
    
    /* set background subtraction */
    public void setSubBG(boolean flag) {
        flagSubBG = flag;
        imc.invalidateAll();
        grid.invalidateAll();
    }
    
    /* set Blur (gaussian) filter */
    public void setBlur(boolean flag) {
        flagBlur = flag;
        imc.invalidateAll();
        grid.invalidateAll();
    }
    
    
    /* set Median filter */
    public void setMedian(boolean flag) {
        flagMedian = flag;
        imc.invalidateAll();
        grid.invalidateAll();
    }
    
    /* toggle Flag for alignment */
    public void toggleAlign() { setAlign(!flagAlign); }
    
    /* toggle FlipX */
    public void toggleFlipX() { setFlipX(!flagFlipX); }
    
    /* toggle FlipY */
    public void toggleFlipY() { setFlipY(!flagFlipY); }
    
    /* toggle FlipZ */
    public void toggleFlipZ() { setFlipZ(!flagFlipZ); }
    
    /* toggle Logarighmic display */
    public void toggleLn() { setLn(!flagLn); }
    
    /* toggle interporation for Z */
    public void toggleInterpZ() { setInterpZ(!flagInterpZ); }
    
    /* toggle subtract background */
    public void toggleSubBG() { setSubBG(!flagSubBG); }
    
    /* toggle Blur (gaussian) filter */
    public void toggleBlur() { setBlur(!flagBlur); }
    
    /* toggle Median filter */
    public void toggleMedian() { setMedian(!flagMedian); }
    
    /* set rotation (by radius unit) */
    public void setRotation(double radian) {
        if (0<=radian && radian<2*Math.PI) { rotation = radian/(4*Math.PI); }
        else { rotation = (((radian/(2*Math.PI)%1)+1)%1)*0.5; }
    }
    
    /* get rotation by raduis unit */
    public double getRotation() {
        return rotation*4*Math.PI;
    }
    
    
     /*
      * Calc n-th image of r[rad] rotated hyperstack.
      * K is a key for cache and encoded by n + 0.5*r', where r' satisfies 0<=r'<=1.
      * Such r' can be obtained by r' = r/(2*pi) for 0<r<2*pi, or
      * mod(mod(r,2*pi)+2*pi,2*pi)/(4*pi) for any r of real value.
      * If r' equals to 0, the cache returns both Aligned and Fliped image.
      * If r' not equals to 0, the cache returns rotated images of that for r=0.
      */
    public Object calcImage(double key) {
        int n = (int)(key); // 1-based index for hyperstack position
        double radian = (key-n)*4*Math.PI; // get rotation radius by floor operation
        if (radian==0) { // no rotation
            //return getProcessorWithoutRotation(n).getPixels();
            return getSubImageInterpolatedZ(n);
        } else { // rotated
            ArrayList al = grid.getUnchecked(radian);
            switch (getBitDepth()){
                case 8:  return getImageRotatedByte(al,n);
                case 16: return getImageRotatedShort(al,n);
                case 32: return getImageRotatedFloat(al,n);
                default: return null;
            }
        }
    }
    
    
    public Object getSubImageInterpolatedZ(int n) {
        if (flagInterpZ) {
            int[] czt = myimp.convertIndexToPosition(n);
            double zNewInRaw = (czt[1]-1) / zscale + 1;
            int zNewInRawFloor = ((int)(zNewInRaw));
            double zMod = zNewInRaw - zNewInRawFloor;
            if (zMod==0 || zNewInRawFloor>=nzorig) {
                // required image is just the image or out of bounds (no data)
                int si1 = getStackIndexCZT(czt[0],zNewInRawFloor,czt[2]);
                return getSubImageEffected(si1);
            } else {
                boolean flagLnOld = flagLn;
                flagLn = false;
                int si1 = getStackIndexCZT(czt[0],zNewInRawFloor,  czt[2]);
                int si2 = getStackIndexCZT(czt[0],zNewInRawFloor+1,czt[2]);
                
                if (flagDebug) {
                    System.out.println("cc="+czt[0]+", cz="+czt[1]+", ct="+czt[2]
                    +", zNewInRaw="+zNewInRaw+", si1="+si1+", si2="+si2);
                }
                
                int w=getWidth(), h=getHeight();
                Object pixels1 = getSubImageEffected(si1);
                Object pixels2 = getSubImageEffected(si2);
                Object pixels3 = interp(pixels1,pixels2,zMod);
                ImageProcessor ip;
                switch (getBitDepth()) {
                    case 8:  ip = new  ByteProcessor(w,h, (byte[])pixels3); break;
                    case 16: ip = new ShortProcessor(w,h); ip.setPixels(pixels3); break;
                    case 32: ip = new FloatProcessor(w,h,(float[])pixels3); break;
                    default: ip = new FloatProcessor(w,h);
                }
                flagLn = flagLnOld;
                if (flagLn) { ip.log(); }
                return ip.getPixels();
            }
        } else { // no interpolation for z
            return getSubImageEffected(n);
        }
    }
    
    
    
    public Object getSubImageEffected(int n) {
        if (flagFlipZ) {
            int[] czt = convertIndexToPositionCZT(n);
            n = getStackIndexCZT(czt[0],nzorig-czt[1]+1,czt[2]);
        }
        if (n<1 || n>nImages)
            throw new IllegalArgumentException("Argument out of range: "+n);
        if (IJ.debugMode) IJ.log("FileInfoVirtualStack: "+n+", "+info[n-1].getOffset());
        //if (n>1) IJ.log("  "+(info[n-1].getOffset()-info[n-2].getOffset()));
        
        if (flagDebug) {
            System.out.println("n="+n);
        }
        
        int w=getWidth(), h=getHeight();
        Object pixels = getSubImage(n);
        ImageProcessor ip;
        
        // System.out.println("getSubImageEffected, bitdepth="+getBitDepth()); ////// for debug
        
        switch (getBitDepth()) {
            case 8:  ip = new  ByteProcessor(w,h, (byte[])pixels); break;
            case 16: ip = new ShortProcessor(w,h); ip.setPixels(pixels); break;
//             case 24: ip = new ColorProcessor(w,h,pixels);
            case 32: ip = new FloatProcessor(w,h,(float[])pixels); break;
            default: ip = new FloatProcessor(w,h);
        }
        if (flagSubBG & radiusSubBG>0) {
            if (getBitDepth()==32) { ip.snapshot(); } // avoid bug in BackgroundSubtractor
            bs.rollingBallBackground(ip,radiusSubBG,false,false,false,true,true);
        }
        if (flagAlign & parallelShift!=null) {
            // int[] czt = myimp.convertIndexToPosition(n);
            int[] czt = convertIndexToPositionCZT(n);
            double sx = parallelShift[0][czt[1]-1][czt[2]-1];
            double sy = parallelShift[1][czt[1]-1][czt[2]-1];
            ip.setInterpolationMethod(this.interpolationMethod);
            ip.translate(sx,sy);
        }
        if (flagMedian) { ip.filter(ip.MEDIAN_FILTER); }
        if (flagBlur) {ip.filter(ip.BLUR_MORE); }
        if (flagFlipX) { ip.flipHorizontal(); }
        if (flagFlipY) { ip.flipVertical(); }
        if (flagLn) { ip.log(); }
        return ip.getPixels();
    }
    
    
    public Object interp(Object p1, Object p2, double zMod) {
        switch (getBitDepth()) {
            case 8:  return interp( (byte[])p1, (byte[])p2,zMod);
            case 16: return interp((short[])p1,(short[])p2,zMod);
            case 32: return interp((float[])p1,(float[])p2,zMod);
            default: return null;
        }
    }
    
    public byte[] interp(byte[] p1, byte[] p2, double zMod) {
        byte[] p3 = new byte[p1.length];
        for (int i=0; i<p1.length; i++) {
            p3[i] = (byte)(((double)(p1[i]))*(1-zMod) + ((double)(p2[i]))*zMod);
        }
        return p3;
    }
    
    public short[] interp(short[] p1, short[] p2, double zMod) {
        short[] p3 = new short[p1.length];
        for (int i=0; i<p1.length; i++) {
            p3[i] = (short)(((double)(p1[i]))*(1-zMod) + ((double)(p2[i]))*zMod);
        }
        return p3;
    }
    
    public float[] interp(float[] p1, float[] p2, double zMod) {
        float[] p3 = new float[p1.length];
        for (int i=0; i<p1.length; i++) {
            p3[i] = (float)(((double)(p1[i]))*(1-zMod) + ((double)(p2[i]))*zMod);
        }
        return p3;
    }
    
    
    public int getStackIndexCZT(int cc, int cz, int ct) {
        // return ncorig*nzorig*(ct-1) + ncorig*(cz-1) + cc-1 + 1;
        int n = ( (ct-1)*nzorig + (cz-1))*ncorig +(cc-1) + 1; // 1-based index
        if (n<=0) { n = 1; }
        return n;
    }
    
    public int[] convertIndexToPositionCZT(int n) {
        int qt = (n-1)/(ncorig*nzorig); // quotient for t
        int n2 = (n-1)-qt*(ncorig*nzorig);
        int qz = n2/ncorig; // quotient for z
        int qc = n2-qz*ncorig;
        if (qc<0) { qc = 0; }
        if (qz<0) { qz = 0; }
        if (qt<0) { qt = 0; }
        int[] czt = {qc+1,qz+1,qt+1};
        return czt;
    }
    
    public ArrayList getGridRotated(double radian) {
        
        int[] dims = myimp.getDimensions();
        int height = dims[1];
        int depth  = dims[3];
        double offsety = ((double)(height))*0.5-0.5;
        double offsetz = ((double)( depth))*0.5-0.5;
        
        Calibration cal = myimp.getCalibration();
        double dy = cal.getY(1);
        double dz = cal.getZ(1);
        double dy_inv = 1/dy;
        double dz_inv = 1/dz;
        
        // obtain rotated grid; assuming x (width) direction as rotation axis;
        double sinr = Math.sin(-radian);
        double cosr = Math.cos(-radian);
        double[] ysin = new double[height];
        double[] ycos = new double[height];
        double[] zsin = new double[depth];
        double[] zcos = new double[depth];
        for (int cy=0; cy<height; cy++) {
            double tmpy = (cy - offsety)*dy;
            ysin[cy] = tmpy * sinr;
            ycos[cy] = tmpy * cosr;
        }
        for (int cz=0; cz<depth; cz++) {
            double tmpz = (cz - offsetz)*dz;
            zsin[cz] = tmpz * sinr;
            zcos[cz] = tmpz * cosr;
        }
        int[][] yg_floor = new int[depth][height];
        int[][] zg_floor = new int[depth][height];
        double[][] yg_mod   = new double[depth][height];
        double[][] zg_mod   = new double[depth][height];
        for (int cz=0; cz<depth; cz++) {
            for (int cy=0; cy<height; cy++) {
                double tmpy = (ycos[cy] - zsin[cz])*dy_inv + offsety;
                double tmpz = (ysin[cy] + zcos[cz])*dz_inv + offsetz;
                if (tmpy<0 || tmpy>=height-1 || tmpz<0 || tmpz>=depth-1) {
                    yg_floor[cz][cy] = -1;
                    zg_floor[cz][cy] = -1;
                    yg_mod[cz][cy] = 0;
                    zg_mod[cz][cy] = 0;
                } else {
                    yg_floor[cz][cy] = ((int)(tmpy)); // cast is faster than Math.floor
                    zg_floor[cz][cy] = ((int)(tmpz)); // cast is faster than Math.floor
                    yg_mod[cz][cy] = tmpy - yg_floor[cz][cy];
                    zg_mod[cz][cy] = tmpz - zg_floor[cz][cy];
                }
            }
        }
        
        ArrayList<Object> al = new ArrayList<Object>(4);
        al.add(yg_floor);
        al.add(zg_floor);
        al.add(yg_mod);
        al.add(zg_mod);
        return al;
    }
    
    public byte[] getImageRotatedByte(ArrayList al, int nSlice) {
        int[][] yg_floor = ((int[][])(al.get(0)));
        int[][] zg_floor = ((int[][])(al.get(1)));
        double[][] yg_mod = ((double[][])(al.get(2)));
        double[][] zg_mod = ((double[][])(al.get(3)));
        
        int[] dims = myimp.getDimensions();
        int width  = dims[0];
        int height = dims[1];
        
        int[] czt = myimp.convertIndexToPosition(nSlice);
        int cc = czt[0]; // 1-based index for channel
        int cz = czt[1]-1; // 0-based index for rotated grid
        int ct = czt[2]; // 1-based index for time
        
        // obtain rotated image
        byte[] im_rot = new byte[width*height];
        for (int cy=0; cy<height; cy++) {
            if (yg_floor[cz][cy]<0 || zg_floor[cz][cy]<0) { continue; }
            double idxSlice1 = myimp.getStackIndex(cc,zg_floor[cz][cy]+1,ct);
            double idxSlice2 = myimp.getStackIndex(cc,zg_floor[cz][cy]+2,ct);
            byte[] im1 = (byte[])(imc.getUnchecked(idxSlice1));
            byte[] im2 = (byte[])(imc.getUnchecked(idxSlice2));
            int offset = yg_floor[cz][cy]*width;
            for (int cx=0; cx<width; cx++) {
                double ip1 = yg_mod[cz][cy]*(im1[cx+width+offset]-im1[cx+offset]) + im1[cx+offset];
                double ip2 = yg_mod[cz][cy]*(im2[cx+width+offset]-im2[cx+offset]) + im2[cx+offset];
                double ip3 = zg_mod[cz][cy]*(ip2-ip1) + ip1;
                im_rot[cy*width+cx] = ((byte)(ip3));
            }
        }
        return im_rot;
    }
    
    public short[] getImageRotatedShort(ArrayList al, int nSlice) {
        int[][] yg_floor = ((int[][])(al.get(0)));
        int[][] zg_floor = ((int[][])(al.get(1)));
        double[][] yg_mod = ((double[][])(al.get(2)));
        double[][] zg_mod = ((double[][])(al.get(3)));
        
        int[] dims = myimp.getDimensions();
        int width  = dims[0];
        int height = dims[1];
        
        int[] czt = myimp.convertIndexToPosition(nSlice);
        int cc = czt[0]; // 1-based index for channel
        int cz = czt[1]-1; // 0-based index for rotated grid
        int ct = czt[2]; // 1-based index for time
        
        // obtain rotated image
        short[] im_rot = new short[width*height];
        for (int cy=0; cy<height; cy++) {
            if (yg_floor[cz][cy]<0 || zg_floor[cz][cy]<0) { continue; }
            double idxSlice1 = myimp.getStackIndex(cc,zg_floor[cz][cy]+1,ct); // z is 1-based index
            double idxSlice2 = myimp.getStackIndex(cc,zg_floor[cz][cy]+2,ct); // z is 1-based index
            short[] im1 = (short[])(imc.getUnchecked(idxSlice1));
            short[] im2 = (short[])(imc.getUnchecked(idxSlice2));
            int offset = yg_floor[cz][cy]*width;
            for (int cx=0; cx<width; cx++) {
                double ip1 = yg_mod[cz][cy]*(im1[cx+width+offset]-im1[cx+offset]) + im1[cx+offset];
                double ip2 = yg_mod[cz][cy]*(im2[cx+width+offset]-im2[cx+offset]) + im2[cx+offset];
                double ip3 = zg_mod[cz][cy]*(ip2-ip1) + ip1;
                im_rot[cy*width+cx] = ((short)(ip3));
            }
        }
        return im_rot;
    }
    
    public float[] getImageRotatedFloat(ArrayList al, int nSlice) {
        int[][] yg_floor = ((int[][])(al.get(0)));
        int[][] zg_floor = ((int[][])(al.get(1)));
        double[][] yg_mod = ((double[][])(al.get(2)));
        double[][] zg_mod = ((double[][])(al.get(3)));
        
        int[] dims = myimp.getDimensions();
        int width  = dims[0];
        int height = dims[1];
        
        int[] czt = myimp.convertIndexToPosition(nSlice);
        int cc = czt[0]; // 1-based index for channel
        int cz = czt[1]-1; // 0-based index for rotated grid
        int ct = czt[2]; // 1-based index for time
        
        // obtain rotated image
        float[] im_rot = new float[width*height];
        for (int cy=0; cy<height; cy++) {
            if (yg_floor[cz][cy]<0 || zg_floor[cz][cy]<0) { continue; }
            double idxSlice1 = myimp.getStackIndex(cc,zg_floor[cz][cy]+1,ct);
            double idxSlice2 = myimp.getStackIndex(cc,zg_floor[cz][cy]+2,ct);
            float[] im1 = (float[])(imc.getUnchecked(idxSlice1));
            float[] im2 = (float[])(imc.getUnchecked(idxSlice2));
            int offset = yg_floor[cz][cy]*width;
            for (int cx=0; cx<width; cx++) {
                double ip1 = yg_mod[cz][cy]*(im1[cx+width+offset]-im1[cx+offset]) + im1[cx+offset];
                double ip2 = yg_mod[cz][cy]*(im2[cx+width+offset]-im2[cx+offset]) + im2[cx+offset];
                double ip3 = zg_mod[cz][cy]*(ip2-ip1) + ip1;
                im_rot[cy*width+cx] = ((float)(ip3));
            }
        }
        return im_rot;
    }
    
    /** Returns the number of images in this stack. */
    public int getSize() {
        return nImages;
    }
    
    /** get info */
    public FileInfo[] getInfo() {
        return info;
    }
    
    /** Returns the label of the Nth image. */
    public String getSliceLabel(int n) {
        if (n<1 || n>nImages)
            throw new IllegalArgumentException("Argument out of range: "+n);
        if (info[0].sliceLabels==null || info[0].sliceLabels.length!=nImages)
            return null;
        else
            return info[0].sliceLabels[n-1];
    }
    
    public int getWidth() {
        return info[0].width;
    }
    
    public int getHeight() {
        return info[0].height;
    }
    
    public int[] getOriginalDimensions() {
        return new int[] {nxorig,nyorig,nzorig,ncorig,ntorig};
    }
    
    public int[] getOriginalDimensionsCZT() {
        return new int[] {nxorig,nyorig,ncorig,nzorig,ntorig};
    }
    
    /*
     * Return pixel values of subimage specified with slice and range.
     * This function works faster than getPixels by the following reasons.
     * 1. This function does not access unneccesary pixels.
     * 2. This function is specialized to get pixel values.
     *    Unnecessary processes (i.e. LUT) are omitted.
     * 3. This function utilizes java.nio.MappedByteBuffer (MemoryMapped I/O).
     */
    public Object getSubImage(int slice, int[] range) {
        // public Object getPixelsPartial_mbb(int slice, int[] range) {
        int bytePerPixel = info[0].getBytesPerPixel();
        int width = info[0].width;
        int height = info[0].height;
        int offset_width = width*bytePerPixel;
        int x_start = range[0];
        int y_start = range[1];
        int x_end = range[2];
        int y_end = range[3];
        int width_subim  = x_end - x_start + 1;
        int height_subim = y_end - y_start + 1;
        int width_subim_byte = width_subim*bytePerPixel;
        long tmpOffset = info[slice-1].longOffset;
        if (tmpOffset==0) { tmpOffset=info[slice-1].getOffset(); }
        
        long offset_to_image = tmpOffset-info[0].getOffset();
        int idx_map       = (int) (offset_to_image / mapping_size);
        int offset_in_map = (int) (offset_to_image % mapping_size);
        
        if (idx_map!=idx_mbb) { updateMappedByteBuffer(idx_map); }
        
        byte[] buf = new byte[width_subim_byte*height_subim];
        int y_loop = y_end - y_start;
        for (int cy=0; cy<=y_loop; cy++) {
            mbb.position(offset_in_map + offset_width*(cy+y_start) + x_start*bytePerPixel);
            if (mbb.remaining()>=width_subim_byte) {
                mbb.get(buf,width_subim_byte*cy,width_subim_byte);
            } else {
                int lenCurrent = mbb.remaining();
                int lenNext    = width_subim_byte - lenCurrent;
                mbb.get(buf,width_subim_byte*cy,lenCurrent);
                idx_map = idx_map+1;
                offset_in_map = 0;
                y_start = -(cy+1);
                updateMappedByteBuffer(idx_map);
                mbb.get(buf,width_subim_byte*cy+lenCurrent+1,lenNext);
            }
        }
        
        switch (bytePerPixel) {
            case 1:
                return buf;
            case 2:
                return convertByte2Short(buf);
            case 4:
                return convertByte2Float(buf);
        }
        return null;
    }
    
    public Object getSubImage(int slice) {
        // if (vsoir!=null) { return vsoir.getProcessor(slice).getPixels(); }
        if (od!=null) { return od.getPixels(slice); }
        
        int bytePerPixel = info[0].getBytesPerPixel();
        int width = info[0].width;
        int height = info[0].height;
        int imsize_byte = width*height*bytePerPixel;
        long tmpOffset = info[slice-1].longOffset;
        if (tmpOffset==0) { tmpOffset=info[slice-1].getOffset(); }
        
        //long offset_to_image = info[slice-1].longOffset-info[0].getOffset();
        long offset_to_image = tmpOffset-info[0].getOffset();
        int idx_map       = (int) (offset_to_image / mapping_size);
        int offset_in_map = (int) (offset_to_image % mapping_size);
        //MappedByteBuffer mbb = (MappedByteBuffer) mappings.get(idx_map);
        
        if (idx_map!=idx_mbb) { updateMappedByteBuffer(idx_map); }
        
        byte[] buf = new byte[imsize_byte];
        
        mbb.position(offset_in_map);
        if (mbb.remaining()>=imsize_byte) { mbb.get(buf); }
        else {
            int lenCurrent = mbb.remaining();
            int lenNext    = imsize_byte - lenCurrent;
            mbb.get(buf,0,lenCurrent);
            idx_map = idx_map+1;
            offset_in_map = 0;
            updateMappedByteBuffer(idx_map);
            //mbb.get(buf,lenCurrent+1,lenNext);
            mbb.get(buf,lenCurrent,lenNext);
        }
        
        switch (bytePerPixel) {
            case 1:
                return buf;
            case 2:
                return convertByte2Short(buf);
            case 4:
                return convertByte2Float(buf);
        }
        return null;
    }
    
    
    /* Convert byte arrays to byte arrays. Derived from ij.io.ImageReader.java */
    public Object convertByte2Byte(byte[] buffer) {
        FileInfo fi = info[0];
        int base = 0;
        int bufferSize = buffer.length;
        int bytesPerPixel = fi.getBytesPerPixel();
        int pixelsRead = bufferSize/bytesPerPixel;
        short[] pixels = new short[pixelsRead];
        if (fi.intelByteOrder) {
            if (fi.fileType==FileInfo.GRAY16_SIGNED)
                for (int i=base,j=0; i<(base+pixelsRead); i++,j+=2)
                    pixels[i] = (short)((((buffer[j+1]&0xff)<<8) | (buffer[j]&0xff))+32768);
            else
                for (int i=base,j=0; i<(base+pixelsRead); i++,j+=2)
                    pixels[i] = (short)(((buffer[j+1]&0xff)<<8) | (buffer[j]&0xff));
        } else {
            if (fi.fileType==FileInfo.GRAY16_SIGNED)
                for (int i=base,j=0; i<(base+pixelsRead); i++,j+=2)
                    pixels[i] = (short)((((buffer[j]&0xff)<<8) | (buffer[j+1]&0xff))+32768);
            else
                for (int i=base,j=0; i<(base+pixelsRead); i++,j+=2)
                    pixels[i] = (short)(((buffer[j]&0xff)<<8) | (buffer[j+1]&0xff));
        }
        return pixels;
    }
    
    
    /* Convert byte arrays to short arrays. Derived from ij.io.ImageReader.java */
    public Object convertByte2Short(byte[] buffer) {
        FileInfo fi = info[0];
        int base = 0;
        int bufferSize = buffer.length;
        int bytesPerPixel = fi.getBytesPerPixel();
        int pixelsRead = bufferSize/bytesPerPixel;
        short[] pixels = new short[pixelsRead];
        if (fi.intelByteOrder) {
            if (fi.fileType==FileInfo.GRAY16_SIGNED)
                for (int i=base,j=0; i<(base+pixelsRead); i++,j+=2)
                    pixels[i] = (short)((((buffer[j+1]&0xff)<<8) | (buffer[j]&0xff))+32768);
            else
                for (int i=base,j=0; i<(base+pixelsRead); i++,j+=2)
                    pixels[i] = (short)(((buffer[j+1]&0xff)<<8) | (buffer[j]&0xff));
        } else {
            if (fi.fileType==FileInfo.GRAY16_SIGNED)
                for (int i=base,j=0; i<(base+pixelsRead); i++,j+=2)
                    pixels[i] = (short)((((buffer[j]&0xff)<<8) | (buffer[j+1]&0xff))+32768);
            else
                for (int i=base,j=0; i<(base+pixelsRead); i++,j+=2)
                    pixels[i] = (short)(((buffer[j]&0xff)<<8) | (buffer[j+1]&0xff));
        }
        return pixels;
    }
    
    /* Convert byte arrays to float arrays. Derived from ij.io.ImageReader.java */
    public Object convertByte2Float(byte[] buffer) {
        FileInfo fi = info[0];
        int base = 0;
        int bufferSize = buffer.length;
        int bytesPerPixel = fi.getBytesPerPixel();
        int pixelsRead = bufferSize/bytesPerPixel;
        float[] pixels = new float[pixelsRead];
        int pmax = base + pixelsRead;
        int tmp;
        int j = 0;
        if (fi.intelByteOrder)
            for (int i=base; i<pmax; i++) {
            tmp = (((buffer[j+3]&0xff)<<24) | ((buffer[j+2]&0xff)<<16) | ((buffer[j+1]&0xff)<<8) | (buffer[j]&0xff));
            if (fi.fileType==FileInfo.GRAY32_FLOAT)
                pixels[i] = Float.intBitsToFloat(tmp);
            else if (fi.fileType==FileInfo.GRAY32_UNSIGNED)
                pixels[i] = (float)(tmp&0xffffffffL);
            else
                pixels[i] = tmp;
            j += 4;
            }
        else
            for (int i=base; i<pmax; i++) {
            tmp = (((buffer[j]&0xff)<<24) | ((buffer[j+1]&0xff)<<16) | ((buffer[j+2]&0xff)<<8) | (buffer[j+3]&0xff));
            if (fi.fileType==FileInfo.GRAY32_FLOAT)
                pixels[i] = Float.intBitsToFloat(tmp);
            else if (fi.fileType==FileInfo.GRAY32_UNSIGNED)
                pixels[i] = (float)(tmp&0xffffffffL);
            else
                pixels[i] = tmp;
            j += 4;
            }
        return pixels;
    }
    
    public void clearAndResizeCache(long maxnum){
        imc = CacheBuilder.newBuilder()
        .maximumSize(maxnum)
        .build( new CacheLoader<Double, Object>() {
            public Object load(Double key) { return calcImage(key); }
        });
    }
    
    public void setMappingSize(int msize) {
        mapping_size = msize*sliceBytesWithGap;
        updateMappedByteBuffer(0);
    }
    
    public void resetMappingSize() {
        mapping_size = (MAPPING_SIZE_MAX/sliceBytesWithGap)*sliceBytesWithGap;
        updateMappedByteBuffer(0);
    }
    
    public void updateMappedByteBuffer(int idx_new) {
        idx_mbb = idx_new;
        long offset = offsetFileInfo + idx_mbb*mapping_size;
        long size2 = Math.min(fileSize-offset, mapping_size);
        try {
            mbb = fch.map(FileChannel.MapMode.READ_ONLY,offset,size2);
        } catch (IOException e) {
            IJ.error("Virtual Stack","Error occurred in setup MappedByteBuffer");
        }
    }
    
}

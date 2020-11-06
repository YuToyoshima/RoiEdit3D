import ij.IJ;
import ij.VirtualStack;
import ij.ImagePlus;
import ij.CompositeImage;
import ij.io.FileInfo;
import ij.measure.Calibration;
import ij.process.ImageProcessor;
import ij.process.FloatProcessor;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.regex.Pattern;
import java.util.regex.Matcher;
import com.google.common.cache.*;


/** This plugin opens following files as a virtual stack.
 *  *.tif (multi-page TIFF)
 *  *.dcv (Hamamatsu Camera Image Data)
 *  *.sif (Andor Camera Image Data)
 *  *.pif (9pipe Image Data)
 *  *.mif (multi-image file)
 *  Files other than the last were handled MyFileInfoVirtualStack.java
 *  directly, and the last was handled by this class
 */
public class MyMultiFileInfoVirtualStack extends VirtualStack {
    
    ImagePlus myimp;
    int nImages;
    public int numChannels, width, height, slices, frames;
    public boolean[] flagFlipX4M, flagFlipY4M, flagRelative;
    public int[] numC, useC, padZ, skip;
    public double[] scaleX,scaleY,rotXY,transX,transY;
    public double[][] gx,gy;
    public String strMif, strPathMif;
    public String[][] strPathImage;
    public MyFileInfoVirtualStack[] arMfivs;
    public FileInfo fi;
    public static String strContent = "\\s*=\\s*([^;]+)\\s*";
    
    
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
    .build( new CacheLoader<Double,Object>() {
        public Object load(Double key) { return calcImage(key); }
    });
    
    public LoadingCache<Double,ArrayList> grid = CacheBuilder.newBuilder()
    .maximumSize(100)
    .build( new CacheLoader<Double,ArrayList>() {
        public ArrayList load(Double radian) { return getGridRotated(radian); }
    });
    
    
    public MyMultiFileInfoVirtualStack(String path, boolean show) {
        strPathMif = path;
        File file = new File(strPathMif);
        String strDir = file.getParent()+file.separator;
        String strName = file.getName();
        strMif = "";
        
        if  (strName.endsWith(".mif")) { // for mif file
            readMif(strPathMif);  // load mif file and get numChannels;            
            init();
            parseMif(); // parse mif parameters
            for (int i=0; i<numChannels; i++) {
                if (strPathImage[i].length == 1) { // single file
                    arMfivs[i] = new MyFileInfoVirtualStack(strPathImage[i][0], false);
                } else { // splitted (multi) file
//                     arMfivs[i] = new MySplitFileInfoVirtualStack(strPathImage[i], false);
                    arMfivs[i] = new MySplitFileInfoVirtualStack(strPathImage[i], false, slices, frames);
                }
            }
        } else { // for single-color (dcv, sif, pif) and multi-color (tiff) images
            MyFileInfoVirtualStack mfivs = 
                new MyFileInfoVirtualStack(strPathMif, false);
            numChannels = mfivs.ncorig; // get numChannels
            init();
            for (int i=0; i<numChannels; i++) { 
                arMfivs[i] = mfivs; 
                numC[i] = numChannels;
                useC[i] = i+1; 
            }
        }
        open();
        if (show) { myimp.show(); }
    }
    
    public void init() {
        strPathImage = new String[numChannels][];
        flagRelative = new boolean[numChannels];
        arMfivs = new MyFileInfoVirtualStack[numChannels];
        
        width = 0;
        height = 0;
        slices = 0;
        frames = 0;
        flagSubBG = false;
        flagFlipX4M = new boolean[numChannels];
        flagFlipY4M = new boolean[numChannels];
        scaleX = new double[numChannels];
        scaleY = new double[numChannels];
        rotXY  = new double[numChannels];
        transX = new double[numChannels];
        transY = new double[numChannels];
        numC   = new int[numChannels];
        useC   = new int[numChannels];
        padZ   = new int[numChannels];
        skip   = new int[numChannels];
        Arrays.fill(scaleX,1);
        Arrays.fill(scaleY,1);
        Arrays.fill(numC,1);
        Arrays.fill(useC,1);
    }
    
    public void open() {
        if (width==0)  { width  = arMfivs[0].getImage().getWidth(); }
        if (height==0) { height = arMfivs[0].getImage().getHeight(); }
        if (slices==0) { slices = arMfivs[0].getOriginalDimensions()[2]; }
        
        setSubBG(flagSubBG);
        setFlipX4M(flagFlipX4M);
        setFlipY4M(flagFlipY4M);
        
        // make fileinfo
        File fileMif = new File(strPathMif);
        fi = new FileInfo();
        fi.width = width;
        fi.height = height;
        fi.offset = 0;
        fi.fileType = fi.GRAY32_FLOAT;
        fi.fileName = fileMif.getName();
        fi.directory = fileMif.getParent()+fileMif.separator;
        fi.nImages = 1;
        
        // take min of nImages of image files as nImages of mif
        nImages = Integer.MAX_VALUE;
        for (int i=0; i<numChannels; i++) {
            if ( !(strPathMif.endsWith(".mif")) || strPathImage[i].length == 1) { // single file
                ImagePlus tmpimp = arMfivs[i].getImage();
                int[] dims = tmpimp.getDimensions();
                //int n = dims[3]*dims[4]; // Z*T
                //int n = dims[3]*dims[4]/numC[i]; // Z*T/numC
                int n = dims[2]*dims[3]*dims[4]/numC[i]; // C*Z*T/numC
                if (nImages>n) { nImages = n; }
            }
        }
        if (nImages == Integer.MAX_VALUE) { // no single file channel 
            nImages = slices*frames; 
        } 
        if (frames == 0) { frames = nImages/slices; }  
        nImages *= numChannels;
        
        calcGridXY();
        setBitDepth(32); // default 32-bit image
        imc.invalidateAll();
        grid.invalidateAll();
        myimp = new ImagePlus(fi.fileName, this);
        myimp.setFileInfo(fi);
        
//         int frames = nImages/(numChannels*slices);
        if (numChannels*slices*frames==nImages) {
            myimp.setDimensions(numChannels, slices, frames);
            myimp.setOpenAsHyperStack(true);
        }
        nxorig = width;
        nyorig = height;
        ncorig = numChannels;
        nzorig = slices;
        ntorig = frames;
        
        for (int i=0; i<numChannels; i++) {
            int nImages = arMfivs[i].getSize();
            int tmpNcorig = nImages/(nzorig*ntorig);
            arMfivs[i].ncorig = tmpNcorig;
            arMfivs[i].nzorig = nzorig;
            arMfivs[i].ntorig = ntorig;
            if(tmpNcorig*nzorig*ntorig==nImages) {
                ImagePlus imp = arMfivs[i].getImage();
                imp.setDimensions(tmpNcorig,nzorig,ntorig);
            }
        }
        
        if (numChannels>1) {
            myimp = new CompositeImage(myimp, CompositeImage.COMPOSITE);
        }
        
    }
    
    
    public void readMif(String path) {
        StringBuilder builder = new StringBuilder();
        try {
            BufferedReader reader = new BufferedReader(new FileReader(path));
            String str = reader.readLine();
            while (str != null){
                builder.append(str + System.getProperty("line.separator"));
                str = reader.readLine();
            }
        } catch(Exception e) { IJ.error("Virtual Stack",""+e); }
        strMif = builder.toString();
        
        String strContent = "\\s*=\\s*([^;]+)\\s*";
        
        // count number of image files
        Matcher mat = Pattern.compile("(?i)(?m)^path\\d*"+strContent).matcher(strMif);
        numChannels = 0;
        while (mat.find()) { numChannels++; }
    }
    
    
    public String getMifParam(String str, String orig) {
        Matcher mat = Pattern.compile("(?i)(?m)^"+str+strContent).matcher(strMif);
        if (mat.find()) { return mat.group(1); }
        else { return orig; }
    }
    
    public double getMifParam(String str, double orig) {
        Matcher mat = Pattern.compile("(?i)(?m)^"+str+strContent).matcher(strMif);
        if (mat.find()) { return Double.parseDouble(mat.group(1)); }
        else { return orig; }
    }
    
    public int getMifParam(String str, int orig) {
        Matcher mat = Pattern.compile("(?i)(?m)^"+str+strContent).matcher(strMif);
        if (mat.find()) { return Integer.parseInt(mat.group(1)); }
        else { return orig; }
    }
    
    public boolean getMifParam(String str, boolean orig) {
        Matcher mat = Pattern.compile("(?i)(?m)^"+str+strContent).matcher(strMif);
        if (mat.find()) { return Boolean.parseBoolean(mat.group(1)); }
        else { return orig; }
    }
    
    
    public void parseMif() {
        // get image path
        File fileMif = new File(strPathMif);
        for (int i=0; i<numChannels; i++) {
            // String substr = getMifParam("path"+(i+1),strPathImage[i]).trim();
            
            strPathImage[i] = getMifParam("path"+(i+1),"").split(",",0);
            for (int j=0; j<strPathImage[i].length; j++) {
                String substr = strPathImage[i][j].trim();
                File file = new File(substr);
                if (!file.isAbsolute()) { // convert relative path to absolute path
                    flagRelative[i] = true;
                    substr = fileMif.getParent()+fileMif.separator+substr;                    
                }
                strPathImage[i][j] = substr;
            }
        }
        
        // get width, height, slices, frames, flagSubBG
        width     = getMifParam("width", width);
        height    = getMifParam("height",height);
        slices    = getMifParam("slices",slices);
        frames    = getMifParam("frames",frames);
        flagSubBG = getMifParam("subbg", flagSubBG);
        
        // get parameters for rigid transformation and crop
        for (int i=0; i<numChannels; i++) {
            flagFlipX4M[i] = getMifParam("flipX"+(i+1), flagFlipX4M[i]);
            flagFlipY4M[i] = getMifParam("flipY"+(i+1), flagFlipY4M[i]);
            scaleX[i]      = getMifParam("scaleX"+(i+1),scaleX[i]);
            scaleY[i]      = getMifParam("scaleY"+(i+1),scaleY[i]);
            rotXY[i]       = getMifParam("rotXY"+(i+1), rotXY[i]);
            transX[i]      = getMifParam("transX"+(i+1),transX[i]);
            transY[i]      = getMifParam("transY"+(i+1),transY[i]);
            numC[i]        = getMifParam("numC"+(i+1),  numC[i]);
            useC[i]        = getMifParam("useC"+(i+1),  useC[i]);
            padZ[i]        = getMifParam("padZ"+(i+1),  padZ[i]);
            skip[i]        = getMifParam("skip"+(i+1),  skip[i]);
        }
    }
    
    
    public void writeMif(String pathSave) {
        try {
            File fileSave = new File(pathSave);
            BufferedWriter bw = new BufferedWriter(new FileWriter(fileSave));
            bw.write("width = "+width+";\r\n");
            bw.write("height = "+height+";\r\n");
            bw.write("slices = "+nzorig+";\r\n");
            bw.write("frames = "+ntorig+";\r\n");
            bw.write("subBG = "+flagSubBG+";\r\n");
            for (int i=0; i<numChannels; i++) {
                String strPath = arMfivs[i].getPath();
                if (flagRelative[i]) {
                    File fileMif = new File(strPathMif);
                    strPath = strPath.replace(fileMif.getParent()+fileMif.separator,"");
                }
                bw.write("path"+(i+1)+" = "+strPath+";\r\n");
                bw.write("flipX"+(i+1)+" = "+flagFlipX4M[i]+";\r\n");
                bw.write("flipY"+(i+1)+" = "+flagFlipY4M[i]+";\r\n");
                bw.write("scaleX"+(i+1)+" = "+scaleX[i]+";\r\n");
                bw.write("scaleY"+(i+1)+" = "+scaleY[i]+";\r\n");
                bw.write("rotXY"+(i+1)+" = "+rotXY[i]+";\r\n");
                bw.write("transX"+(i+1)+" = "+transX[i]+";\r\n");
                bw.write("transY"+(i+1)+" = "+transY[i]+";\r\n");
                bw.write("numC"+(i+1)+" = "+numC[i]+";\r\n");
                bw.write("useC"+(i+1)+" = "+useC[i]+";\r\n");
                bw.write("padZ"+(i+1)+" = "+padZ[i]+";\r\n");
                bw.write("skip"+(i+1)+" = "+skip[i]+";\r\n");
            }
            bw.close();
        } catch(Exception e) { IJ.error("Virtual Stack",""+e); }
    }
    
    
    public ImagePlus getImage() {
        return myimp;
    }
    
    public String getPath() { return strPathMif; }
    
    /** Returns an ImageProcessor for the specified image,
     * were 1<=n<=nImages. Returns null if the stack is empty.
     */
    public ImageProcessor getProcessor(int n) {
        if (n<1 || n>nImages) {
            throw new IllegalArgumentException("Argument out of range: "+n);
        }
        if (flagDebug) { System.out.println("getProcessor, n="+n); }
        MyFileOpener mfo = new MyFileOpener(fi);
        ImagePlus imp = mfo.open(false,imc.getUnchecked(n+rotation));
        return imp.getProcessor();
    }
    
    
    public void deleteSlice(int n) {
        IJ.error("Virtual Stack","Cannot delete slice !");
        return;
    }
    
    
    public void clearAndResizeCache(long maxnum){
        imc = CacheBuilder.newBuilder()
        .maximumSize(maxnum)
        .build( new CacheLoader<Double, Object>() {
            public Object load(Double key) { return calcImage(key); }
        });
        for (int i=0; i<numChannels; i++) { arMfivs[i].clearAndResizeCache(maxnum); }
    }
    
    public void setMappingSize(int msize) {
        for (int i=0; i<numChannels; i++) { arMfivs[i].setMappingSize(msize); }
    }
    
    public void resetMappingSize() {
        for (int i=0; i<numChannels; i++) { arMfivs[i].resetMappingSize(); }
    }
    
    
    
    public void setFlipX4M(boolean[] flag) {
        flagFlipX4M = flag;
        for (int i=0; i<numChannels; i++) { arMfivs[i].setFlipX(flag[i]); }
        imc.invalidateAll();
    }
    
    public void setFlipY4M(boolean[] flag) {
        flagFlipY4M = flag;
        for (int i=0; i<numChannels; i++) { arMfivs[i].setFlipY(flag[i]); }
        imc.invalidateAll();
    }
    
    public void setScaleX(double[] in) { scaleX = in; calcGridXY(); imc.invalidateAll(); }
    public void setScaleY(double[] in) { scaleY = in; calcGridXY(); imc.invalidateAll(); }
    public void setRotXY(double[] in)  {  rotXY = in; calcGridXY(); imc.invalidateAll(); }
    public void setTransX(double[] in) { transX = in; calcGridXY(); imc.invalidateAll(); }
    public void setTransY(double[] in) { transY = in; calcGridXY(); imc.invalidateAll(); }
    public void setNumC(int[] in)      {   numC = in; imc.invalidateAll(); }
    public void setUseC(int[] in)      {   useC = in; imc.invalidateAll(); }
    public void setPadZ(int[] in)      {   padZ = in; imc.invalidateAll(); }
    public void setSkip(int[] in)      {   skip = in; imc.invalidateAll(); }
    public void setWidth(int in)   { width = in;  fi.width = in; imc.invalidateAll(); }
    public void setHeight(int in) { height = in; fi.height = in; imc.invalidateAll(); }
    
    
    public boolean[] getFlipX4M() { return flagFlipX4M; }
    public boolean[] getFlipY4M() { return flagFlipY4M; }
    public double[] getScaleX()  { return scaleX; }
    public double[] getScaleY()  { return scaleY; }
    public double[] getRotXY()  { return rotXY; }
    public double[] getTransX() { return transX; }
    public double[] getTransY() { return transY; }
    public int[] getNumC()      { return numC; }
    public int[] getUseC()      { return useC; }
    public int[] getPadZ()      { return padZ; }
    public int[] getSkip()      { return skip; }
    public int getWidth()       { return width;  }
    public int getHeight()      { return height; }
    
    public int getNumSlices()   { return nzorig; }
    public int getNumChannels() { return numChannels; }
    public int getNumFrames()   { return ntorig; }
    
    public int getSize()        { return nImages; }
    public String getSliceLabel(int n) { return null; }
    
    public int[] getOriginalDimensions() {
        return new int[] {nxorig,nyorig,nzorig,ncorig,ntorig};
    }
    
    public int[] getOriginalDimensionsCZT() {
        return new int[] {nxorig,nyorig,ncorig,nzorig,ntorig};
    }
    
    
    ////// set parmeters common in MyFileInfoVirtualStack //////
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
    
    /* set FlipZ */
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
//         FileInfo[] tmpinfo = new FileInfo[nImages];
//         for (int i=0; i<nImages; i++) {
//             if (i<info.length) { tmpinfo[i] = info[i]; }
//             else               { tmpinfo[i] = info[i%dims[2]]; }
//         }
//         info = tmpinfo;
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
        for (int i=0; i<numChannels; i++) { arMfivs[i].setRadiusSubBG(r); }
        imc.invalidateAll();
        grid.invalidateAll();
    }
    
    /* set background subtraction */
    public void setSubBG(boolean flag) {
        flagSubBG = flag;
        for (int i=0; i<numChannels; i++) { arMfivs[i].setSubBG(flag); }
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
    
    
    public double[][][] getParallelShift() { return parallelShift; }
    public boolean getAlign() { return flagAlign; }
    public boolean getFlipX() { return flagFlipX; }
    public boolean getFlipY() { return flagFlipY; }
    public boolean getFlipZ() { return flagFlipZ; }
    public boolean getLn() { return flagLn; }
    public boolean getInterpZ() { return flagInterpZ; }
    public double getRadiusSubBG() { return radiusSubBG; }
    public boolean getSubBG() { return flagSubBG; }
    public boolean getBlur() { return flagBlur; }
    public boolean getMedian() { return flagMedian; }
    
    
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
    
    public float[] calcImage(double key) {
        int n = (int)(key); // 1-based index for hyperstack position
        double radian = (key-n)*4*Math.PI; // get rotation radius by floor operation
        if (radian==0) { // no rotation
            return getSubImageInterpolatedZ(n);
        } else { // rotated
            ArrayList al = grid.getUnchecked(radian);
            return getImageRotatedFloat(al,n);
        }
    }
    
    
    public float[] getSubImageInterpolatedZ(int n) {
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
                
                float[] pixels1 = getSubImageEffected(si1);
                float[] pixels2 = getSubImageEffected(si2);
                float[] pixels3 = interp(pixels1,pixels2,zMod);
                FloatProcessor ip = new FloatProcessor(width,height,pixels3);
                flagLn = flagLnOld;
                if (flagLn) { ip.log(); }
                return (float[])(ip.getPixels());
            }
        } else { // no interpolation for z
            return getSubImageEffected(n);
        }
    }
    
    
    
    public float[] getSubImageEffected(int n) {
        if (flagFlipZ) {
            int[] czt = convertIndexToPositionCZT(n);
            n = getStackIndexCZT(czt[0],nzorig-czt[1]+1,czt[2]);
        }
        if (n<1 || n>nImages)
            throw new IllegalArgumentException("Argument out of range: "+n);
        if (flagDebug) { System.out.println("n="+n); }
        
        FloatProcessor ip = getSubImageProcessor(n);
//         if (flagSubBG & radiusSubBG>0) {
//             if (getBitDepth()==32) { ip.snapshot(); } // avoid bug in BackgroundSubtractor
//             bs.rollingBallBackground(ip,radiusSubBG,false,false,false,true,true);
//         }
        if (flagAlign & parallelShift!=null) {
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
        return (float[])(ip.getPixels());
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
        
        ArrayList al = new ArrayList(4);
        al.add(yg_floor);
        al.add(zg_floor);
        al.add(yg_mod);
        al.add(zg_mod);
        return al;
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
    
    
    public FloatProcessor getSubImageProcessor(int key) {
        int channel = (key-1)%numChannels;
        int n       = (key-channel-1)/numChannels+1;
        MyFileInfoVirtualStack mfivs = arMfivs[channel];
        // int nq = (n-1)*mfivs.ncorig+useC[channel];
        int nq = (n-1)*numC[channel]+useC[channel];
        
        if (flagDebug) {
            System.out.println("getSubImageProcessor; key="+key+
            ", channel="+channel+", n="+n+", nq="+nq+", skip="+skip[channel]);
        }
        
        // if z is smaller than padZ, image is padded with 0.
        if (padZ[channel] > 0) {
            int[] czt = mfivs.convertIndexToPositionCZT(nq);
            if (czt[1]<=padZ[channel]) {
                return new FloatProcessor(width,height);
            }
        }
        
        Object obj = mfivs.imc.getUnchecked((double)nq+skip[channel]);
        float[] pixels;
        switch(mfivs.getBitDepth()) {
            case 8:
                byte[] arrByte = (byte[])obj;
                pixels = new float[arrByte.length];
                for (int i=0; i<arrByte.length; i++) {
                    pixels[i] = (float)arrByte[i];
                }
                break;
            case 16:
                FileInfo[] info = mfivs.getInfo();
                boolean flagUnsigned = (info[0].fileType == FileInfo.GRAY16_UNSIGNED);
                short[] arrShort = (short[])obj;
                pixels = new float[arrShort.length];
//                 for (int i=0; i<arrShort.length; i++) {
//                     if (flagUnsigned && arrShort[i]<0) {
//                         pixels[i] = (float)arrShort[i] + 65536; // convert to unsigned
//                     } else {
//                         pixels[i] = (float)arrShort[i];
//                     }
//                 }
                if (flagUnsigned) {
                    for (int i=0; i<arrShort.length; i++) {
                        pixels[i] = (float) (arrShort[i] & 0xFFFF); // convert to unsigned
                    }
                } else {
                    for (int i=0; i<arrShort.length; i++) {
                        pixels[i] = (float) arrShort[i];
                    }
                }
                break;
            default: pixels = (float[])obj;
        }
        int w = mfivs.getWidth();
        int h = mfivs.getHeight();
        float[] pixels2 = rigid(w,h,pixels,channel);
        FloatProcessor ip = new FloatProcessor(width,height,pixels2);
//         if (flagFlipX4M[channel]) { ip.flipHorizontal(); }
//         if (flagFlipY4M[channel]) { ip.flipVertical(); }
        return ip;
    }
    
    
    public float[] rigid(int nx, int ny, float[] pixels, int channel) {
        float[] retpixels = new float[width*height];
        
        for (int j=0; j<height; j++) {
            for (int i=0; i<width; i++) {
                double cx = gx[channel][j*width+i];
                double cy = gy[channel][j*width+i];
                int fx = (int)(Math.floor(cx));
                int fy = (int)(Math.floor(cy));
                double dx = cx - fx;
                double dy = cy - fy;
                
                if( fx<0 || fx>nx-1 || fy<0 || fy>ny-1 ) {
                    retpixels[j*width+i] = 0; // zero fill
                    continue;
                }
                
                double p00=0, p01=0, p10=0, p11=0;
                p00 = pixels[(fy+0)*nx+(fx+0)];
                if (fx+1<nx) { p01 = pixels[(fy+0)*nx+(fx+1)]; }
                if (fy+1<ny) { p10 = pixels[(fy+1)*nx+(fx+0)]; }
                if (fx+1<nx && fy+1<ny) { p11 = pixels[(fy+1)*nx+(fx+1)]; }
                
                double p0 = (1-dx)*p00 + dx*p01;
                double p1 = (1-dx)*p10 + dx*p11;
                retpixels[j*width+i] = (float)((1-dy)*p0 + dy*p1);
            }
        }
        return retpixels;
    }
    
    
    public void calcGridXY() {
        gx = new double[numChannels][width*height];
        gy = new double[numChannels][width*height];
        for (int i=0; i<numChannels; i++) {
            double sint = Math.sin(rotXY[i]);
            double cost = Math.cos(rotXY[i]);
            double isx = 1/scaleX[i];
            double isy = 1/scaleY[i];
            for (int cy=0; cy<height; cy++) {
                for (int cx=0; cx<width; cx++) {
                    double tmpx1 = cx - transX[i];
                    double tmpy1 = cy - transY[i];
                    double tmpx2 =  tmpx1*cost + tmpy1*sint;
                    double tmpy2 = -tmpx1*sint + tmpy1*cost;
                    gx[i][cy*width+cx] = tmpx2*isx;
                    gy[i][cy*width+cx] = tmpy2*isy;
                }
            }
        }
    }
    
    
}
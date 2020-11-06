import ij.ImagePlus;
import ij.VirtualStack;
import ij.io.FileInfo;
import java.io.IOException;
import java.nio.ByteBuffer;
import jp.co.olympus.viewer.exception.IDAException;
import jp.co.olympus.viewer.ui.ImageCloseHandler;


/*
 * .oir file handler for MyFileInfoVirtualStack
 * assuming followings:
 *   isUsingVirtualStack = true;
 *   isStitching = false;
 *   group, level, and area are all 0; (no group, no level, no area)
 */

public class OirDecoder {
    public String directory, name;
    public IDADecoder id;
    public IDAFileInfo fi;
    public CustomImagePlus imp;
    public ImagePlus locImp;
    public CustomVirtualStack stack;
    public String key = "1_1";
    public Object[] keylist;
    public IDA_HAREA hArea;
    
    public OirDecoder(String directory, String name) throws IOException {
        this.directory = directory;
        this.name = name;
        
        try {
            // from OlympusViewer _Viewer.java run()
            id = new IDADecoder(directory, name);
            id.open();
            IDAFileInfo[][][] tmpIdaFileInfo = id.getIDAFileInfo();
            IDAFileInfo[] idaFileInfo = tmpIdaFileInfo[0][0];
            fi = idaFileInfo[0];
            id.InitTiffDimMap(key);
            
            // from OlympusViewer _Viewer.java openImage()
            id.populateImageInformation(idaFileInfo);
            int group = 0;
            int level = 0;
            int i = 0; // loop counter for areaCount
            int areaMappedtoRowCol[][] = new int[1][1];
            CustomBoolean isValidArea = new CustomBoolean();
            stack = id.createVirtualStack(group,level,i,fi,areaMappedtoRowCol,isValidArea);
            
            // from OlympusViewer _Viewer.java showImage()
            int area = 0;
            ImageCloseHandler ich = new ImageCloseHandler();
            imp = new CustomImagePlus(name,stack,id,fi,group,level,area);
            ImagePlus.addImageListener(ich);
            imp.setOpenAsHyperStack(true);
            imp.setFileInfo(fi);
            id.setCalibration(imp,fi);
            imp.setLUTs(id.luts);
            locImp = new CustomComposite(imp);
            
            // for getPixels
            keylist = id.GLtiffDimMap.get(key).keySet().toArray();
            hArea = id.getArea(0,0,0);
        }
        catch (IDAException e) {
            e.printStackTrace();
            throw new IOException("IDAEception occrued in OirDecoder",e);
        }
    }
    
    public FileInfo getFileInfo() { return (FileInfo)fi; }
    public IDAFileInfo getIDAFileInfo() { return fi; }
    public IDADecoder getIDADecoder() { return id; }
    public ImagePlus getImage() { return locImp; }
    public CustomImagePlus getCustomImagePlus() { return imp; }
    public VirtualStack getVirtualStack() {return (VirtualStack)stack; }
    public CustomVirtualStack getCustomVirtualStack() {return stack; }
    public IDA_HIMAGE getNewIdaHimage() { return new IDA_HIMAGE(); }
    
    public Object getPixels(int n) {
        try {
            AxisChannelStuct acs = id.GLtiffDimMap.get(key).get(keylist[n-1].toString());
            IDA_HIMAGE hImage = new IDA_HIMAGE();
            IDA_AXIS_INFO[] axes = acs.axes;
            String channelId = acs.channelId;
            id.idalJniLayer.GetImage(id.hAccessor,hArea,channelId,axes,id.nAxes,hImage);
            ByteBuffer byteBuffer = id.getImageBody(hImage, fi);
            byteBuffer.position(0);
            Object pixels = id.getPixels(fi,byteBuffer);
            return pixels;
        } catch (Exception e) {
            e.printStackTrace();
            return null;
        }
    }
    
    public void close() { id.idalJniLayer.ReleaseArea(id.hAccessor,hArea); }
    
}
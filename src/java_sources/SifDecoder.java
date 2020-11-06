import ij.io.*;
import java.io.*;
import java.util.ArrayList;
import org.apache.commons.lang3.ArrayUtils;


// copyed from ij.io.TiffDecoder and modified for sif format
public class SifDecoder {
    
    // metadata types
    static final String MAGIC_NUMBER = "Andor Technology Multi-Channel File";        
    
    String directory;
    String name;
    RandomAccessStream in;
    boolean littleEndian = true;
    
    public SifDecoder(String directory, String name) {
        this.directory = directory;
        this.name = name;
    }
    
    public FileInfo[] getSifInfo() throws IOException {
        if (in==null)
            in = new RandomAccessStream(new RandomAccessFile(new File(directory, name), "r"));        
        String strFirst = readLine();
        if (!strFirst.equals(MAGIC_NUMBER)) {
            System.err.println("The file seems not to be a SIF file. Exit.");
            in.close();
            return null;
        }
        
        String strQuery = "Pixel number";
        while (true) {
            String strCurrent = readLine();
            if (strCurrent == null) {
                System.err.println("The file seems not to be a SIF file. Exit.");
                in.close();
                return null;
            }
            if (strCurrent.startsWith(strQuery)) { break; }
        }
        readLine(); // ignore one line
        String[] strImageArea = readLine().split(" "); // will return 10 elements
        String[] strFrameArea = readLine().split(" "); // will return 7 elements
        int pixelPerImage = Integer.parseInt(strImageArea[9]);
        long pixelInVideo = Long.parseLong(strImageArea[8]);
        int leftPixel   = Integer.parseInt(strFrameArea[1]);
        int topPixel    = Integer.parseInt(strFrameArea[2]);
        int rightPixel  = Integer.parseInt(strFrameArea[3]);
        int bottomPixel = Integer.parseInt(strFrameArea[4]);
        int[] imageArea = new int[] {
            Integer.parseInt(strImageArea[2]),
            Integer.parseInt(strImageArea[5]),
            Integer.parseInt(strImageArea[7]),
            Integer.parseInt(strImageArea[4]),
            Integer.parseInt(strImageArea[3]),
            Integer.parseInt(strImageArea[6]) };
        int[] frameArea = new int[] {leftPixel,bottomPixel,rightPixel,topPixel};
        int vBin = Integer.parseInt(strFrameArea[6]);
        int hBin = Integer.parseInt(strFrameArea[5]);
        int[] frameBins = new int[] {vBin,hBin};
        
        int width  = (rightPixel - leftPixel + 1)/hBin;
        int height = (topPixel - bottomPixel + 1)/vBin;
        int pixelPerFrame = width*height;
        int f0 = Integer.parseInt(strImageArea[7]);
        int f1 = Integer.parseInt(strImageArea[6]);
        int numFrames = f1-f0+1;
        in.skip(numFrames*11+2+numFrames*21); // skip time stamp
        
        
        FileInfo fi = new FileInfo();
        fi.intelByteOrder = littleEndian;
        fi.offset = in.getFilePointer();
        fi.width  = width;
        fi.height = height;
        fi.nImages = numFrames;
        fi.fileType = FileInfo.GRAY32_FLOAT; // fixed to float image
        fi.stripLengths = new int[] { fi.width*fi.height*fi.getBytesPerPixel() };
        fi.samplesPerPixel = 1;
        fi.rowsPerStrip = fi.height;
        fi.fileFormat = FileInfo.TIFF;
        fi.fileName = name;
        fi.directory = directory;
        fi.compression = FileInfo.COMPRESSION_NONE;
        fi.whiteIsZero = false;
        
        fi.lutSize = 0;
        fi.gapBetweenImages = 0;
        
        in.seek(0);
        fi.inputStream = in;
        
        return new FileInfo[] {fi};
        
    }
    
    public int getInt() throws IOException {
        int b1 = in.read();
        int b2 = in.read();
        int b3 = in.read();
        int b4 = in.read();
        if (littleEndian)
            return ((b4 << 24) + (b3 << 16) + (b2 << 8) + (b1 << 0));
        else
            return ((b1 << 24) + (b2 << 16) + (b3 << 8) + b4);
    }
    
    public long getUnsignedInt() throws IOException {
        return (long)getInt()&0xffffffffL;
    }
    
    public int getShort() throws IOException {
        int b1 = in.read();
        int b2 = in.read();
        if (littleEndian)
            return ((b2<<8) + b1);
        else
            return ((b1<<8) + b2);
    }
    
    public long readLong() throws IOException {
        if (littleEndian)
            return ((long)getInt()&0xffffffffL) + ((long)getInt()<<32);
        else
            return ((long)getInt()<<32) + ((long)getInt()&0xffffffffL);
    }
    
    public double readDouble() throws IOException {
        return Double.longBitsToDouble(readLong());
    }
    
    public String readLine() throws IOException {
        Byte[] bytes = readLineByte();
        if (bytes==null) { return null; } // error reading
        else if (bytes.length==0 || bytes[0]==null) { return ""; } // string of zero length
        else return new String(ArrayUtils.toPrimitive(bytes));
    }
    
    public Byte[] readLineByte() throws IOException {
        ArrayList<Byte> list = new ArrayList<Byte>();
        int b;
        while (true) {
            b = in.read();
            if (b==0x0A || b<0) { break; }
            else {list.add((byte)b);}
        }
        if (b<0 && list.size()==0) { return null; } // error reading
        else { return list.toArray(new Byte[1]); }
    }
    
}

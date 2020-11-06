import ij.io.*;
import java.io.*;
import java.util.ArrayList;
import org.apache.commons.lang3.ArrayUtils;


// copyed from ij.io.TiffDecoder and modified for dcv format
public class DcvDecoder {
    
    // metadata types
    static final String MAGIC_NUMBER = "DCIMG";
    
    String directory;
    String name;
    RandomAccessStream in;
    boolean littleEndian = true;
    
    public DcvDecoder(String directory, String name) {
        this.directory = directory;
        this.name = name;
    }
    
    public FileInfo[] getDcvInfo() throws IOException {
        if (in==null)
            in = new RandomAccessStream(new RandomAccessFile(new File(directory, name), "r"));
        String strFirst = readLine();
        if (!strFirst.startsWith(MAGIC_NUMBER)) {
            System.err.println("The file seems not to be a DCV file. Exit.");
            in.close();
            return null;
        }
        
        in.seek(11);
        int version = (int) (getUnsignedInt());
        int numberFrames, width, height, bytePerFrame, gapBetweenImages;
        
        if (version==0) {
            gapBetweenImages = 0;
            in.seek(36);
            numberFrames = (int) (getUnsignedInt());
            //int offsetAddress = (int) (getUnsignedInt());
            //in.seek(48);
            //int fileSize = (int) (getUnsignedInt());
            //in.seek(120);
            //int offsetFromHeader = (int) (getUnsignedInt());
            in.seek(164);
            width = (int) (getUnsignedInt());
            in.seek(172);
            height = (int) (getUnsignedInt());
            bytePerFrame = (int) (getUnsignedInt());
            //in.seek(192);
            //int offsetImageEnd = ((int) (getUnsignedInt())) + offsetAddress;
            in.seek(232);
        } else { // add the section By Yu Toyohshima, 2018/6/6
            gapBetweenImages = 32;
            in.seek(36);
            numberFrames = (int) (getUnsignedInt());
            in.seek(184);
            width = (int) (getUnsignedInt());
            in.seek(188);
            height = (int) (getUnsignedInt());
            in.seek(196);
            bytePerFrame = (int) (getUnsignedInt());
            in.seek(1024);
        }
        
                
        FileInfo fi = new FileInfo();
        fi.intelByteOrder = littleEndian;
        fi.offset = in.getFilePointer();
        //fi.width  = width;
        //fi.width = 2048/(1200/width);
        fi.width = bytePerFrame/(height*2); // captured width is 2048 or 1024 (2x2 binning)
        fi.height = height;
        fi.nImages = numberFrames;
        fi.fileType = FileInfo.GRAY16_UNSIGNED; // fixed to float image        
        fi.stripLengths = new int[] { fi.width*fi.height*fi.getBytesPerPixel() };
        fi.samplesPerPixel = 1;
        fi.rowsPerStrip = fi.height;
        fi.fileFormat = FileInfo.TIFF;
        fi.fileName = name;
        fi.directory = directory;
        fi.compression = FileInfo.COMPRESSION_NONE;
        fi.whiteIsZero = false;        
        
        fi.description = "givenwidth="+width+"\n";
        
        fi.lutSize = 0;
        fi.gapBetweenImages = gapBetweenImages;
        
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

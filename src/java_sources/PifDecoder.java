import ij.io.*;
import java.io.*;
import java.util.ArrayList;
import org.apache.commons.lang3.ArrayUtils;


// copyed from ij.io.TiffDecoder and modified for pif format
public class PifDecoder {
    
    // metadata types
    static final int MAGIC_NUMBER = 0x50494600;  // "PIF "
    
    String directory;
    String name;
    RandomAccessStream in;
    boolean littleEndian = true;
    
    public PifDecoder(String directory, String name) {
        this.directory = directory;
        this.name = name;
    }
    
    public FileInfo[] getPifInfo() throws IOException {
        if (in==null)
            in = new RandomAccessStream(new RandomAccessFile(new File(directory, name), "r"));
        if (in.readInt()!=MAGIC_NUMBER) {
            System.err.println("The file seems not to be a PIF file. Exit.");
            in.close();
            return null;
        }
        FileInfo fi = new FileInfo();
        fi.intelByteOrder = littleEndian;
        fi.offset = (int) (getUnsignedInt());
        fi.width  = (int) (getUnsignedInt());
        fi.height = (int) (getUnsignedInt());
        long nz = getUnsignedInt();
        long t0 = getUnsignedInt();
        long t1 = getUnsignedInt();
        long nt = t1-t0;
        fi.nImages = (int) (nt*nz);
        
        int pixelType = (int) getUnsignedInt();
        // long PixelSize = getUnsignedInt(); // ignored
        
        switch (pixelType) {
            case 0:
                fi.fileType = FileInfo.GRAY8;
                break;
            case 1:
                System.err.println("ImageJ cannot handle signed 8-bit correctly.");
                fi.fileType = FileInfo.GRAY8;
                break;
            case 2:
                fi.fileType = FileInfo.GRAY16_UNSIGNED;
                break;
            case 3:
                fi.fileType = FileInfo.GRAY16_SIGNED;
                break;
            case 4:
                fi.fileType = FileInfo.GRAY32_INT;
                break;
            case 5:
                fi.fileType = FileInfo.GRAY32_FLOAT;
                break;
            case 6:
                fi.fileType = FileInfo.GRAY64_FLOAT;
                break;
            default:
                System.out.println("PifDecoder cannot handle the image type.");
                fi.fileType = FileInfo.BITMAP;
                break;
        }
        fi.stripLengths = new int[] { fi.width*fi.height*fi.getBytesPerPixel() };
        fi.samplesPerPixel = 1;
        fi.rowsPerStrip = fi.height;
        fi.fileFormat = FileInfo.TIFF;
        fi.fileName = name;
        fi.directory = directory;
        fi.compression = FileInfo.COMPRESSION_NONE;
        fi.whiteIsZero = false;
        
        fi.description = "ImageJ=1.47v\nimages="+fi.nImages+"\nchannels=1\n"
                +"slices="+nz+"\nframes="+nt+"\nhyperstack=true\n";
        
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

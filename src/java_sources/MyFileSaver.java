//package ij.io;
import ij.io.*;
import java.io.*;
import ij.*;

public class MyFileSaver {
    public FileInfo fi;
    public DataOutputStream out;

    public MyFileSaver() {    
    }

    public boolean saveAsTiffStack_open(String path,int[] imageSize) {
        return saveAsTiffStack_open(path,imageSize,FileInfo.GRAY32_FLOAT);
    }

    public boolean saveAsTiffStack_open(String path,int[] imageSize,int imageType) {
        // fi.virtualStack = null;
        // fi.info = imp.getInfoProperty();
        // fi.description = getDescriptionString();
        // if (imp.isComposite()) saveDisplayRangesAndLuts(imp, fi);
        
        int width    = imageSize[0];
        int height   = imageSize[1];
        int channels = imageSize[2];
        int slices   = imageSize[3];
        int frames   = imageSize[4];
        
        fi = new FileInfo();
        fi.width   = width;
        fi.height  = height;
        fi.nImages = channels*slices*frames;
        fi.fileType = imageType;

        StringBuffer sb = new StringBuffer(100);
        sb.append("ImageJ="+ImageJ.VERSION+"\n");
        if (fi.nImages>1 && fi.fileType!=FileInfo.RGB48)
            sb.append("images="+fi.nImages+"\n");
        if (channels>1)
            sb.append("channels="+channels+"\n");
        if (slices>1)
            sb.append("slices="+slices+"\n");
        if (frames>1)
            sb.append("frames="+frames+"\n");
        sb.append("hyperstack=true\n");
        sb.append("mode=composite\n");
        sb.append((char)0);
        fi.description = new String(sb);
        
        try {
            MyTiffEncoder file = new MyTiffEncoder(fi);
            out = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(path)));
            file.write_singleIFD(out);
            return true;
        } catch (IOException e) {
            String msg = e.getMessage();
            if (msg.length()>100)
            msg = msg.substring(0, 100);
            IJ.error("MyFileSaver","An error occured writing the file.\n \n" + msg);
            return false;
        }        
    }

    public boolean saveAsTiffStack_append(Object obj,int nImages) {
        int nImages_orig = fi.nImages;
        fi.pixels = obj;
        fi.nImages = nImages;
        try {
            new ij.io.ImageWriter(fi).write(out);
            fi.nImages = nImages_orig;
            return true;
        } catch (IOException e) {
            fi.nImages = nImages_orig;
            String msg = e.getMessage();
            if (msg.length()>100)
            msg = msg.substring(0, 100);
            IJ.error("MyFileSaver","An error occured writing the file.\n \n" + msg);
            return false;
        }
    }

    public boolean saveAsTiffStack_close(){
        try {
            out.close();
            return true;
        } catch (IOException e) {
            String msg = e.getMessage();
            if (msg.length()>100)
            msg = msg.substring(0, 100);
            IJ.error("MyFileSaver","An error occured writing the file.\n \n" + msg);
            return false;
        }        
    }
    
}

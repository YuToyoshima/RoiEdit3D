import ij.ImagePlus;
import ij.IJ;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;
import java.io.File;

public class OpenMIF_ implements PlugIn {

	public void run(String arg) {
		OpenDialog od = new OpenDialog("Open MIF/SIF/DCV/PIF/TIF...", arg);
		String file = od.getFileName();
		if (file == null) return;
		String directory = od.getDirectory();
	    MyMultiFileInfoVirtualStack mmfivs = new MyMultiFileInfoVirtualStack(directory+File.separator+file,false);
		ImagePlus imp = mmfivs.getImage();
        if (imp != null ) {
			imp.show();
		} else {
			IJ.showMessage("Open MIF...", "Failed.");
		}
	}
}
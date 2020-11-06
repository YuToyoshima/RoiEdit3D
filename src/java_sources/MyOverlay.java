//package ij.gui;
import ij.gui.*;
import java.util.*;
import java.awt.*;

public class MyOverlay extends ij.gui.Overlay {

    public static final int XY_MODE=0, XZ_MODE=1, YZ_MODE=2, PROJ_MODE=3;

    public TreeMap<String,Roi> roimap = new TreeMap<String,Roi>();
    public ArrayList<Roi> roilist = new ArrayList<Roi>();
    public boolean flagChanged = true; // if true, re-construct roilist before use it.
    public MyRoiManager2 myRoiManager;
    public int mode;

    // constructors
    public MyOverlay(MyRoiManager2 myRoiManager, int mode) {
        super();
        this.myRoiManager = myRoiManager;
        this.mode = mode;
    }

    // methods
    public void put(String name, Roi roi) {
        roimap.put(name,roi);
        flagChanged = true;
    }

    public void remove(String name) {
        roimap.remove(name);
        flagChanged = true;
    }

    public TreeMap<String,Roi> getRoimap() {
        return roimap;
    }

    public void setRoimap(TreeMap<String,Roi> roimap) {
        this.roimap = roimap;
        flagChanged= true;
    }

    @Override
    public void clear() {
        super.clear();
        roimap.clear();
        flagChanged = true;
    }

    public void reconstruct() {
        //ij.IJ.log("MyOverlay.reconstruct_before(); mode="+mode+", size="+roimap.size());
        flagChanged = false;
        roilist.clear();
        roilist.addAll(Arrays.asList(super.toArray())); // Overlay.list
        roilist.addAll(roimap.values());
        //ij.IJ.log("MyOverlay.reconstruct_after(); mode="+mode+", size="+roimap.size());
    }

    @Override
    public Roi get(int index) {
        if (flagChanged) { reconstruct(); }
        return roilist.get(index);
    }

    @Override
    public int size() {
        //ij.IJ.log("MyOverlay.size(); mode="+mode+", size="+roimap.size());
        if (myRoiManager.updateRoimap(mode)) { reconstruct(); }
        return roilist.size();
    }

    @Override
    public Roi[] toArray() {
        if (flagChanged) { reconstruct(); }
        int tmpsize = roilist.size();
        Roi[] array = new Roi[tmpsize];
        System.arraycopy(roilist.toArray(array),0,array,0,tmpsize);
        return array;
    }

    @Override
    public String toString() {
        if (flagChanged) { reconstruct(); }
        StringBuilder sb = new StringBuilder(super.toString());
        return sb.append(roilist.toArray()).toString();
    }
}

//package ij.gui;
import ij.gui.*;
import com.google.common.collect.*;
import ij.*;
import java.util.*;
import java.awt.*;
import org.apache.commons.lang3.ArrayUtils;


/*
This class uses guava (ver 18.0).
Matlab(R2012b) uses old google-collect.jar and cannot use guava collectly.
If you use this class through Matlab, put guava-18.0.jar to the front of
java static class path in Matlab.
 
From R2012b, classes can be added to the front of the static Java class
path by including the following line in javaclasspath.txt:
<before>
 
Make javaclasspath.txt that contains following
-----------FROM HERE
<before>
C:\matlab\R2012b\3rd_party\guava-18.0.jar
-----------TO HERE
and put it where the return value of prefdir() (matlab function).
 
Ref:
http://stackoverflow.com/questions/16366059/best-way-to-override-matlabs-default-static-javaclasspath
 
 *
 * roitableがTreeBasedTableになっているが、ここをキャッシュつきにしたい。
 * MyTableModel3とうまく紐付けて、on demandに必要なROIだけを生成・破棄するとよい。
 *
 */

public class MyRoiManager {
    // static variables
    public static final int XY_MODE=0, XZ_MODE=1, YZ_MODE=2, PROJ_MODE=3;
    
    // class variables
    public TreeBasedTable<Integer,String,MyEllipsoid> roitable = TreeBasedTable.create();
    //public TreeMultimap<Integer,String> selectiontable = TreeMultimap.create();
    public TreeBasedTable<Integer,String,Integer> selectiontable = TreeBasedTable.create();
    public TreeMap<String,Roi> roimap_xy, roimap_xz, roimap_yz, roimap_proj;
    public Color strokeColor    = Color.YELLOW;
    public Color selectionColor = Color.CYAN;
    public MyOverlay myOverlay_xy, myOverlay_xz, myOverlay_yz, myOverlay_proj;
    public MyOrthogonal_Views mov;
    public MyEllipsoid me_lastSelection;
    public ImagePlus imp, imp_proj;
    public int modeCZT;
    public boolean[] active_last;
    public int[] nx_last = new int[4];
    public int[] ny_last = new int[4];
    public int[] nz_last = new int[4];
    public int[] nt_last = new int[4];
    public boolean[] flagForceUpdate = new boolean[4];
    
    
    public MyRoiManager(ImagePlus imp) {
        this.imp = imp;
        ImageCanvas ic = imp.getCanvas();
        myOverlay_xy = new MyOverlay(this, XY_MODE);
        roimap_xy = myOverlay_xy.getRoimap();
        if (ic!=null) { ic.setShowAllList(myOverlay_xy); }
        setModeCZT();
        updateAndDraw();
    }
    
    public MyRoiManager(MyOrthogonal_Views mov) {
        this.mov = mov;
        this.imp = mov.getImage();
        myOverlay_xy = new MyOverlay(this, XY_MODE);
        myOverlay_xz = new MyOverlay(this, XZ_MODE);
        myOverlay_yz = new MyOverlay(this, YZ_MODE);
        roimap_xy = myOverlay_xy.getRoimap();
        roimap_xz = myOverlay_xz.getRoimap();
        roimap_yz = myOverlay_yz.getRoimap();
        mov.setOverlays(myOverlay_xy,myOverlay_xz,myOverlay_yz);
        setModeCZT();
        updateAndDraw();
    }
    
    
    public void setProj(ImagePlus imp_proj) {
        this.imp_proj = imp_proj;
        ImageCanvas ic = imp_proj.getCanvas();
        myOverlay_proj = new MyOverlay(this, PROJ_MODE);
        roimap_proj = myOverlay_proj.getRoimap();
        if (ic!=null) { ic.setShowAllList(myOverlay_proj); }
        updateAndDraw();
    }
    
    
    public void add(MyEllipsoid[] mes) {
        for (int p=0; p<mes.length; p++) {
            add_inner(mes[p]);
        }
        updateAndDraw();
    }
    
    public void add(MyEllipsoid me) {
        add_inner(me);
        updateAndDraw();
    }
    
    public void add_inner(MyEllipsoid me) {
        put_inner(me.frame,me.name,me);
    }
    
    public void put(int[] frames, String[] names, MyEllipsoid[] mes) {
        for (int p=0; p<frames.length; p++) {
            put_inner(frames[p],names[p],mes[p]);
        }
        updateAndDraw();
    }
    
    public void put (int frame, String name, MyEllipsoid me) {
        put_inner(frame,name,me);
        updateAndDraw();
    }
    
    public void put_inner(int frame, String name, MyEllipsoid me) {
        me.setModeCZT(modeCZT);
        roitable.put(frame,name,me);
    }
    
    public MyEllipsoid get(int frame, String name) {
        return roitable.get(frame,name);
    }
    
    public void remove(int[] frames, String[] names) {
        for (int p=0; p<frames.length; p++) {
            remove_inner(frames[p],names[p]);
        }
        updateAndDraw();
    }
    
    public void remove(int frame, String name) {
        remove_inner(frame,name);
        updateAndDraw();
    }
    
    public void remove_inner(int frame, String name) {
        MyEllipsoid tmproi = roitable.get(frame,name);
        if (tmproi != null) {
            roitable.remove(frame,name);
            if (roimap_xy!=null)   {   roimap_xy.remove(name); }
            if (roimap_xz!=null)   {   roimap_xz.remove(name); }
            if (roimap_yz!=null)   {   roimap_yz.remove(name); }
            if (roimap_proj!=null) { roimap_proj.remove(name); }
        }
        selectiontable.remove(frame,name);
    }
    
    public void removeAll(int[] frame) {
        for (int i=0; i<frame.length; i++) {
            removeAll_inner(frame[i]);
        }
        updateAndDraw();
    }
    
    public void removeAll(int frame) {
        removeAll_inner(frame);
        updateAndDraw();
    }
    
    public void removeAll_inner(int frame) {
        roitable.row(frame).clear();
        //selectiontable.removeAll(frame);
        selectiontable.row(frame).clear();
    }
    
    public void clear() {
        roitable.clear();
        if (myOverlay_xy!=null)   {   myOverlay_xy.clear(); }
        if (myOverlay_xz!=null)   {   myOverlay_xz.clear(); }
        if (myOverlay_yz!=null)   {   myOverlay_yz.clear(); }
        if (myOverlay_proj!=null) { myOverlay_proj.clear(); }
        updateAndDraw();
    }
    
    
    // update MyEllipseRoi and return wheteher updated or not.
    public boolean updateRoimap(int mode) {
        if (imp==null) { return false; }
        int nx = 1;
        int ny = 1;
        int nz = 1;
        int nt = imp.getT();
        boolean[] active;
        boolean flagActiveUpdated = true;
        //ij.IJ.log("MyRoiManager.update().check; mode="+mode+",roitable.size()="+roitable.row(nt).size());
        
        if (imp.isComposite()) {
            active = ((CompositeImage)imp).getActiveChannels();
        } else {
            active = new boolean[7];
            active[0] = true;
        }
        if (active_last != null) {
            flagActiveUpdated= (   active[0] != active_last[0]
            || active[1] != active_last[1]
            || active[2] != active_last[2]
            || active[3] != active_last[3]
            || active[4] != active_last[4]
            || active[5] != active_last[5]
            || active[6] != active_last[6] );
        }
        active_last = active;
        
        // judge whether update is necessary
        if (mov!=null) { // z-stack
            int[] loc = mov.getCrossLoc(); // 0-based positions of [x,y,z]
            nx = loc[0]+1;
            ny = loc[1]+1;
            nz = loc[2]+1;
        }
        if (   nx == nx_last[mode]
        && ny == ny_last[mode]
        && nz == nz_last[mode]
        && nt == nt_last[mode]
        && !flagActiveUpdated
        && !flagForceUpdate[mode]) {
            return false;
        }
        nx_last[mode] = nx;
        ny_last[mode] = ny;
        nz_last[mode] = nz;
        nt_last[mode] = nt;
        
        // update
        //ij.IJ.log("MyRoiManager.update().update; mode="+mode);
        Collection<MyEllipsoid> tmproitable = roitable.row(nt).values();
        Iterator<MyEllipsoid> iter = tmproitable.iterator();
        MyEllipsoid tmpMyEllipsoid;
        MyEllipseRoi tmpMyEllipseRoi;
        TreeMap<String,Roi> tmpRoimap;
        int pos;
        
        switch (mode) {
            case PROJ_MODE:
                tmpRoimap = roimap_proj;
                pos = 1;
                break;
            case XY_MODE:
                tmpRoimap = roimap_xy;
                pos = nz;
                break;
            case XZ_MODE:
                tmpRoimap = roimap_xz;
                pos = ny;
                break;
            case YZ_MODE:
                tmpRoimap = roimap_yz;
                pos = nx;
                break;
            default:
                return false;
        }
        
        tmpRoimap.clear();
        while(iter.hasNext()) {
            tmpMyEllipsoid = iter.next();
            int cc = tmpMyEllipsoid.getCPosition();
            if ( cc==0 || active[cc-1] ) {
                tmpMyEllipseRoi = tmpMyEllipsoid.getMyEllipseAt(pos,mode);
                if (tmpMyEllipseRoi != null) {
                    tmpRoimap.put(tmpMyEllipseRoi.getNameOrig(),tmpMyEllipseRoi);
                    //ij.IJ.log("MyRoiManager.update().put; mode="+mode+", size="+tmpRoimap.size());
                }
            }
        }
        
        //ij.IJ.log("MyRoiManager.update().finalize; mode="+mode);
        flagForceUpdate[mode] = false;
        return true;
    }
    
    
    
    public void setNameDisp(int frame, String name, String nameDisp) {
        roitable.get(frame,name).setNameDisp(nameDisp);
    }
    
    public void setNameDisp(int[] frame, String[] name, String[] nameDisp) {
        for (int p=0; p<frame.length; p++) {
            setNameDisp(frame[p],name[p],nameDisp[p]);
        }
        updateAndDraw();
    }
    
    
    public void setFlagUseNameDisp(int frame, String name, boolean tf) {
        roitable.get(frame,name).setFlagUseNameDisp(tf);
    }
    
    public void setFlagUseNameDisp(int[] frame, String[] name, boolean[] tf) {
        for (int p=0; p<frame.length; p++) {
            setFlagUseNameDisp(frame[p],name[p],tf[p]);
        }
        updateAndDraw();
    }
    
    public void setFlagUseNameDisp(boolean tf) {
        Collection<MyEllipsoid> tmproitable = roitable.values();
        Iterator<MyEllipsoid> iter = tmproitable.iterator();
        MyEllipsoid tmpMyEllipsoid;
        while(iter.hasNext()) {
            tmpMyEllipsoid = iter.next();
            tmpMyEllipsoid.setFlagUseNameDisp(tf);
        }
        updateAndDraw();
    }
    
    public void setFlagUseFixedName(int frame, String name, boolean tf) {
        roitable.get(frame,name).setFlagUseFixedName(tf);
    }
    
    public void setFlagUseFixedName(int[] frame, String[] name, boolean[] tf) {
        for (int p=0; p<frame.length; p++) {
            setFlagUseFixedName(frame[p],name[p],tf[p]);
        }
        updateAndDraw();
    }
    
    public void setFlagUseFixedName(boolean tf) {
        Collection<MyEllipsoid> tmproitable = roitable.values();
        Iterator<MyEllipsoid> iter = tmproitable.iterator();
        MyEllipsoid tmpMyEllipsoid;
        while(iter.hasNext()) {
            tmpMyEllipsoid = iter.next();
            tmpMyEllipsoid.setFlagUseFixedName(tf);
        }
        updateAndDraw();
    }
    
    
    
    public void setColor(Color strokeColor, Color selectionColor) {
        this.strokeColor    = strokeColor;
        this.selectionColor = selectionColor;
        
        if (roitable.size()==0) { return; }
        Iterator<MyEllipsoid> liter = roitable.values().iterator();
        while (liter.hasNext()) {
            liter.next().setStrokeColor(strokeColor);
        }
        
        if (selectiontable.size()!=0) {
            //Iterator<Map.Entry<Integer,String>> siter = selectiontable.entries().iterator();
            Iterator<Table.Cell<Integer,String,Integer>> siter = selectiontable.cellSet().iterator();
            while (siter.hasNext()) {
                //Map.Entry<Integer,String> e = siter.next();
                Table.Cell<Integer,String,Integer> e = siter.next();
                //MyEllipsoid me = roitable.get(e.getKey(),e.getValue());
                MyEllipsoid me = roitable.get(e.getRowKey(),e.getColumnKey());
                if (me!=null) { me.setStrokeColor(selectionColor); }
            }
        }
        
        updateAndDraw();
    }
    
    
    public void addSelection(int[] frames, String[] names,int[] rowindices) {
        for (int p=0; p<frames.length; p++) {
            addSelection_inner(frames[p],names[p],rowindices[p]);
        }
        moveToLastSelection();
        updateAndDraw();
    }
    
    public void addSelection(int[] frames, String[] names) {
        for (int p=0; p<frames.length; p++) {
            addSelection_inner(frames[p],names[p],-1);
        }
        moveToLastSelection();
        updateAndDraw();
    }
    
    public void addSelection(int frame, String name, int rowindex) {
        addSelection_inner(frame,name,rowindex);
        moveToLastSelection();
        updateAndDraw();
    }
    
    public void addSelection(int frame, String name) {
        addSelection_inner(frame,name,-1);
        moveToLastSelection();
        updateAndDraw();
    }
    
    public void addSelection_inner(int frame, String name, int rowindex) {
        // if (selectiontable.containsEntry(frame,name)) { return; }
        if (selectiontable.contains(frame,name)) { return; }
        selectiontable.put(frame,name,rowindex);
        MyEllipsoid me = roitable.get(frame,name);
        if (me!=null) {
            me.setStrokeColor(selectionColor);
            me_lastSelection = me;
        }
    }
    
    public void removeSelection(int[] frames, String[] names) {
        for (int p=0; p<frames.length; p++) {
            this.removeSelection_inner(frames[p],names[p]);
        }
        updateAndDraw();
    }
    
    public void removeSelection(int frame, String name) {
        removeSelection_inner(frame,name);
        updateAndDraw();
    }
    
    public void removeSelection_inner(int frame, String name) {
        selectiontable.remove(frame,name);
        MyEllipsoid me = roitable.get(frame,name);
        if (me!=null) { me.setStrokeColor(strokeColor); }
    }
    
    public void clearSelection() {
        if (selectiontable.size()==0) { return; }
        //Iterator<Map.Entry<Integer,String>> siter = selectiontable.entries().iterator();
        Iterator<Table.Cell<Integer,String,Integer>> siter = selectiontable.cellSet().iterator();
        while (siter.hasNext()) {
            //Map.Entry<Integer,String> e = siter.next();
            Table.Cell<Integer,String,Integer> e = siter.next();
            //MyEllipsoid me = roitable.get(e.getKey(),e.getValue());
            MyEllipsoid me = roitable.get(e.getRowKey(),e.getColumnKey());
            if (me!=null) { me.setStrokeColor(strokeColor); }
        }
        selectiontable.clear();
        updateAndDraw();
        
    }
    
    public boolean isSelect(int frame, String name) {
        //return selectiontable.containsEntry(frame,name);
        return selectiontable.contains(frame,name);
    }
    
    public double[] getSelectedRowIndices() {
        int[] indices =  ArrayUtils.toPrimitive(
        selectiontable.values().toArray(new Integer[0]));
        return convertToDouble(indices);
    }
    
    public double[] getSelectedRowIndices(int frame){
        int[] indices = ArrayUtils.toPrimitive(
        selectiontable.row(frame).values().toArray(new Integer[0]));
        return convertToDouble(indices);
    }
    
    public double[] convertToDouble(int[] in) {
        double[] ret = new double[in.length];
        for (int p=0; p<in.length; p++) { ret[p] = (double)(in[p]); }
        return ret;
    }
    
    public String[] getSelectedRowNames() {
        return selectiontable.columnKeySet().toArray(new String[0]);
    }
    
    public String[] SelectedRowNames(int frame){
        return selectiontable.row(frame).keySet().toArray(new String[0]);
    }
    
    public MyOverlay[] getOverlays() {
        MyOverlay[] ret = new MyOverlay[4];
        ret[0] = myOverlay_xy;
        ret[1] = myOverlay_xz;
        ret[2] = myOverlay_yz;
        ret[3] = myOverlay_proj;
        return ret;
    }
    
    public void setModeCZT() {
        //    z    t    c    modeCZT
        //    o    o    o     4
        //    o    o    x     3
        //    o    x    o     2
        //    o    x    x    -1
        //    x    o    o     1
        //    x    o    x    -2
        //    x    x    o    -3
        //    x    x    x    -4
        if (imp.getNSlices()>1) {
            if (imp.getNFrames()>1) {
                if (imp.getNChannels()>1) {     modeCZT =  4; }
                else {                modeCZT =  3; }
            } else {
                if (imp.getNChannels()>1) {     modeCZT =  2; }
                else {                modeCZT = -1; }
            }
        } else {
            if (imp.getNFrames()>1) {
                if (imp.getNChannels()>1) {     modeCZT =  1; }
                else {                modeCZT = -2; }
            } else {
                if (imp.getNChannels()>1) {     modeCZT = -3; }
                else {                modeCZT = -4; }
            }
        }
    }
    
    
    
    public void showAll(boolean flagShowAll) {
        mov = MyOrthogonal_Views.getInstance();
        if (mov!=null) {
            if (flagShowAll) {
                mov.setOverlays(myOverlay_xy,myOverlay_xz,myOverlay_yz);
            } else {
                mov.setOverlays(new MyOverlay(this,XY_MODE),
                new MyOverlay(this,XZ_MODE),
                new MyOverlay(this,YZ_MODE));
            }
        } else {
            ImageCanvas ic = imp.getCanvas();
            if (ic!=null) {
                if (flagShowAll) { ic.setShowAllList(myOverlay_xy); }
                else             { ic.setShowAllList(new MyOverlay(this,XY_MODE)); }
            }
        }
        
        if (imp_proj!=null) {
            ImageCanvas ic_proj = imp_proj.getCanvas();
            if (ic_proj!=null) {
                if(flagShowAll) { ic_proj.setShowAllList(myOverlay_proj); }
                else            { ic_proj.setShowAllList(new MyOverlay(this,PROJ_MODE)); }
            }
        }
        updateAndDraw();
    }
    
    public void showName(boolean flagShowName) {
        if(myOverlay_xy!=null)   {   myOverlay_xy.drawBackgrounds(flagShowName); }
        if(myOverlay_xz!=null)   {   myOverlay_xz.drawBackgrounds(flagShowName); }
        if(myOverlay_yz!=null)   {   myOverlay_yz.drawBackgrounds(flagShowName); }
        if(myOverlay_proj!=null) { myOverlay_proj.drawBackgrounds(flagShowName); }
        if(myOverlay_xy!=null)   {   myOverlay_xy.drawNames(flagShowName); }
        if(myOverlay_xz!=null)   {   myOverlay_xz.drawNames(flagShowName); }
        if(myOverlay_yz!=null)   {   myOverlay_yz.drawNames(flagShowName); }
        if(myOverlay_proj!=null) { myOverlay_proj.drawNames(flagShowName); }
        updateAndDraw();
    }
    
    public void updateAndDraw() { // force update
        flagForceUpdate[0] = true;
        flagForceUpdate[1] = true;
        flagForceUpdate[2] = true;
        flagForceUpdate[3] = true;
        if(mov != null) { mov.imageUpdated(imp); }
        else { imp.updateAndDraw(); }
        if (imp_proj!=null) { imp_proj.updateAndDraw(); }
    }
    
    public void moveToLastSelection(){
        if (me_lastSelection==null) { return; }
        if (modeCZT>0) {
            //imp.setC(me_lastSelection.channel);
            imp.setT(me_lastSelection.frame);
            if (mov != null) {
                mov.setCrossLoc((int)(me_lastSelection.xpos),
                (int)(me_lastSelection.ypos),
                (int)(me_lastSelection.zpos));
            } else {
                imp.setZ(me_lastSelection.slice);
            }
        } else {
            imp.setPosition(me_lastSelection.position);
        }
        if (imp_proj != null) {
            //imp_proj.setC(me_lastSelection.channel);
            imp_proj.setT(me_lastSelection.frame);
        }
        me_lastSelection = null;
    }
    
}

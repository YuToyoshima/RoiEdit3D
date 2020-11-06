//package ij.gui;
import ij.gui.*;
import java.awt.*;
import java.awt.image.*;
import ij.*;
import org.apache.commons.lang3.math.NumberUtils;


/** Elliptical region of interest. */
/** The method "equals" was modified */
public class MyEllipseRoi extends EllipseRoi {
    public int channel_orig, slice_orig, frame_orig;
    public String nameDisp;
    public boolean flagUseNameDisp = false; // if true, getName() returns nameDisp
    public boolean flagUseFixedName = false; // if true, getName() returns nameDisp only if isNumberString(name)
    
    public MyEllipseRoi(double x1, double y1, double x2, double y2, double aspectRatio) {
        super(x1,y1,x2,y2,aspectRatio);
    }
    
    @Override
    /** Checks whether two rectangles are equal. */
    public boolean equals(Object obj) {
        if (obj instanceof MyEllipseRoi) {
            MyEllipseRoi roi2 = (MyEllipseRoi)obj;
            if (this.getCPosition() != roi2.getCPosition()) { return false; }
            if (this.getZPosition() != roi2.getZPosition()) { return false; }
            if (this.getTPosition() != roi2.getTPosition()) { return false; }
            if (this.getPosition()  != roi2.getPosition())  { return false; }
            if (!this.getName().equals(roi2.getName()))     { return false; }
            return true;
        } else {
            return false;
        }
    }
    
    @Override
    public void setPosition(int channel, int slice, int frame) {
        super.setPosition(channel,slice,frame);
        channel_orig = channel;
        slice_orig   = slice;
        frame_orig   = frame;
    }
    
    public int getCPositionOrig() {
        return channel_orig;
    }
    
    public int getZPositionOrig() {
        return slice_orig;
    }
    
    public int getTPositionOrig() {
        return frame_orig;
    }
    
    @Override
    /* to move the roi should be prohbited */
    public int isHandle(int sx, int sy) {
        return -1;
    }
    
    @Override
    /* to select the roi should be prohibited */
    public boolean contains(int x, int y) {
        return false;
    }
    
    @Override
    /* to set name and nameDisp at once */
    public void setName(String s) {
        setNameOrig(s);
        setNameDisp(s);
    }
    
    @Override
    /* return name or nameDisp depending on flagUseNameDisp and flagUseFixeName */
    public String getName() {
        if (flagUseNameDisp && nameDisp!=null && 
            !(flagUseFixedName && !NumberUtils.isNumber(super.getName()))) {
            return getNameDisp();
        } else {
            return getNameOrig();
        }
    }
    
    public void setNameDisp(String s) {
        nameDisp = s;
    }
    
    public String getNameDisp() {
        return nameDisp;
    }
    
    public void setNameOrig(String s) {
        super.setName(s);
    }
    
    public String getNameOrig() {
        return super.getName();
    }
    
    public void setFlagUseNameDisp(boolean tf) {
        flagUseNameDisp = tf;
    }

    public boolean getFlagUseNameDisp() {
        return flagUseNameDisp;
    }
    
    public void setFlagUseFixedName(boolean tf) {
        flagUseFixedName = tf;
    }

    public boolean getFlagUseFixedName() {
        return flagUseFixedName;
    }
}

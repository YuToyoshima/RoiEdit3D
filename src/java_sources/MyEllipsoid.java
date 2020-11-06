//package ij.gui;
import ij.gui.*;
import java.util.*;
import java.awt.Color;
import ij.*;
import org.apache.commons.lang3.math.NumberUtils;

public class MyEllipsoid implements Comparable {
    // static variables
    public static final int XY_MODE=0, XZ_MODE=1, YZ_MODE=2, PROJ_MODE=3;

    // class variables
    public MyEllipseRoi roi_xy, roi_xz, roi_yz, roi_proj;
    public Color color = Color.YELLOW;
    public String name, nameDisp;
    public boolean flagUseNameDisp = false; // if true, getName() returns nameDisp
    public boolean flagUseFixedName = false; // if true, getName() returns nameDisp only if isNumberString(name)
    
    public double[] params_xy, params_xz, params_yz;
    public double xpos, ypos, zpos, 
                  x_last, y_last, z_last,
                  xc_disp, yc_disp, zc_disp, 
                  xc_slice, yc_slice, zc_slice,
                  nsigma2;
    public int channel, frame, slice, position=0;
    public int modeCZT = 4; // specify image dimension
        //    z    t    c    modeCZT
        //    o    o    o     4
        //    o    o    x     3
        //    o    x    o     2
        //    o    x    x    -1
        //    x    o    o     1
        //    x    o    x    -2
        //    x    x    o    -3
        //    x    x    x    -4


    // constructor
    public MyEllipsoid(double[] mu, double[] sigma, double nsigma,
            double z_ratio, int channel, int frame, String name) {
        this.channel = channel;
        this.frame = frame;
        init(mu,sigma,nsigma,z_ratio);
        setName(name);
    }

    public MyEllipsoid(double[] mu, double[] sigma, double nsigma, 
            double z_ratio, String name) {
        this(mu,sigma,nsigma,z_ratio, 0, 1, name);        
    }

    public MyEllipsoid(double[] mu, double[] sigma, double nsigma,
            double z_ratio, int channel, int frame, String name, Color color) {
        this(mu,sigma,nsigma,z_ratio,channel,frame,name);
        this.setStrokeColor(color);
    }
    
    
    public MyEllipsoid(double[] mu, double[] sigma, double nsigma,
            double z_ratio, int channel, int frame, String name,
            String nameDisp) {
        this(mu,sigma,nsigma,z_ratio,channel,frame,name);
        setNameDisp(nameDisp);
    }

    public MyEllipsoid(double[] mu, double[] sigma, double nsigma, 
            double z_ratio, String name, String nameDisp) {
        this(mu,sigma,nsigma,z_ratio,name);
        setNameDisp(nameDisp);
    }

    public MyEllipsoid(double[] mu, double[] sigma, double nsigma,
            double z_ratio, int channel, int frame, String name, 
            String nameDisp, Color color) {
        this(mu,sigma,nsigma,z_ratio,channel,frame,name,color);
        setNameDisp(nameDisp);
    }

    
    // static methods (multi-create)
    public static MyEllipsoid[] create(double[][] mu, double[][] sigma, double nsigma,
            double z_ratio, int channel, int[] frame, String[] name, Color color) {
        int siz = frame.length;
        MyEllipsoid[] ret = new MyEllipsoid[siz];
        for (int p=0; p<siz; p++) {
            ret[p] = new MyEllipsoid(mu[p],sigma[p],nsigma,
                        z_ratio,channel,frame[p],name[p],color);
        }
        return ret;
    }

    public static MyEllipsoid[] create(double[][] mu, double[][] sigma, double nsigma,
            double z_ratio, int channel, int[] frame, String[] name) {
        int siz = frame.length;
        MyEllipsoid[] ret = new MyEllipsoid[siz];
        for (int p=0; p<siz; p++) {
            ret[p] = new MyEllipsoid(mu[p],sigma[p],nsigma,z_ratio,
                                        channel,frame[p],name[p]);
        }
        return ret;
    }

    public static MyEllipsoid create(double[] mu, double[] sigma, double nsigma,
             double z_ratio, int channel, int frame, String name, Color color) {        
        return new MyEllipsoid(mu,sigma,nsigma,z_ratio,channel,frame,name,color);
    }

    public static MyEllipsoid create(double[] mu, double[] sigma, double nsigma,
            double z_ratio, int channel, int frame, String name) {        
        return new MyEllipsoid(mu,sigma,nsigma,z_ratio,channel,frame,name);
    }

    
    public static MyEllipsoid[] create(double[][] mu, double[][] sigma, double nsigma,
            double z_ratio, int channel, int[] frame, String[] name, 
            String[] nameDisp, Color color) {
        int siz = frame.length;
        MyEllipsoid[] ret = new MyEllipsoid[siz];
        for (int p=0; p<siz; p++) {
            ret[p] = new MyEllipsoid(mu[p],sigma[p],nsigma,
                        z_ratio,channel,frame[p],name[p],nameDisp[p],color);
        }
        return ret;
    }

    public static MyEllipsoid[] create(double[][] mu, double[][] sigma, double nsigma,
            double z_ratio, int channel, int[] frame, String[] name, String[] nameDisp) {
        int siz = frame.length;
        MyEllipsoid[] ret = new MyEllipsoid[siz];
        for (int p=0; p<siz; p++) {
            ret[p] = new MyEllipsoid(mu[p],sigma[p],nsigma,z_ratio,
                                        channel,frame[p],name[p],nameDisp[p]);
        }
        return ret;
    }

    public static MyEllipsoid create(double[] mu, double[] sigma, double nsigma,
             double z_ratio, int channel, int frame, String name, String nameDisp, 
             Color color) {        
        return new MyEllipsoid(mu,sigma,nsigma,z_ratio,channel,frame,name,nameDisp,color);
    }

    public static MyEllipsoid create(double[] mu, double[] sigma, double nsigma,
            double z_ratio, int channel, int frame, String name, String nameDisp) {        
        return new MyEllipsoid(mu,sigma,nsigma,z_ratio,channel,frame,name,nameDisp);
    }
    

    // methods
    public void init(double[] mu, double[] sigma, double nsigma, double z_ratio) {

        nsigma2 = 2*nsigma*nsigma;

        xpos = mu[0]-1; // convert from 1-base to 0-base
        ypos = mu[1]-1; // convert from 1-base to 0-base
        zpos = mu[2]-1; // convert from 1-base to 0-base
        slice = (int)(Math.round(zpos))+1;

        xc_disp = xpos + 0.5; // for display
        yc_disp = ypos + 0.5; // for display
        zc_disp = zpos + 0.5; // for display

        xc_slice = xpos + 1; // for slice calculation
        yc_slice = ypos + 1; // for slice calculation
        zc_slice = zpos + 1; // for slice calculation

        double S_11 = sigma[0];
        double S_12 = sigma[1];
        double S_13 = sigma[2];
        double S_22 = sigma[3];
        double S_23 = sigma[4];
        double S_33 = sigma[5];

        double detS =    S_11*S_22*S_33 
                     + 2*S_12*S_13*S_23 
                     -   S_11*S_23*S_23 
                     -   S_12*S_12*S_33 
                     -   S_13*S_22*S_13;
        double detS_inv = 1/detS;

        double Sinv_11 = detS_inv*(S_22*S_33 - S_23*S_23);
        double Sinv_22 = detS_inv*(S_11*S_33 - S_13*S_13);
        double Sinv_33 = detS_inv*(S_11*S_22 - S_12*S_12);
        double Sinv_12 = detS_inv*(S_13*S_23 - S_12*S_33);
        double Sinv_13 = detS_inv*(S_12*S_23 - S_13*S_22);
        double Sinv_23 = detS_inv*(S_12*S_13 - S_11*S_23);

        double xrmax = Math.ceil(Math.sqrt(nsigma2*S_11)) + 2; // magrin = 2
        double yrmax = Math.ceil(Math.sqrt(nsigma2*S_22)) + 2; // magrin = 2
        double zrmax = Math.ceil(Math.sqrt(nsigma2*S_33)) + 2; // magrin = 2

        // for XY image
        params_xy = getParams(Sinv_11,Sinv_22,Sinv_12,Sinv_13,Sinv_23,Sinv_33);
        params_xy[11] = zrmax;
        
        // for XZ and YZ image, z-ratio should be applied.
        double z_ratio_inv = 1/z_ratio;
        Sinv_13 *= z_ratio_inv;
        Sinv_23 *= z_ratio_inv;
        Sinv_33 *= z_ratio_inv*z_ratio_inv;
        zc_disp *= z_ratio;

        // for XZ image
        params_xz = getParams(Sinv_11,Sinv_33,Sinv_13,Sinv_12,Sinv_23,Sinv_22);
        params_xz[11] = yrmax;
        
        // for YZ image
        params_yz = getParams(Sinv_22,Sinv_33,Sinv_23,Sinv_12,Sinv_13,Sinv_11);
        params_yz[11] = xrmax;
        
        x_last = 0;
        y_last = 0;
        z_last = 0;
    }
    
    
    public double[] getParams(double Sinv_11, double Sinv_22, double Sinv_12,
                               double Sinv_13, double Sinv_23, double Sinv_33) {
        // diagonalization by orthogonal matrix ( convert to x'-y' coordinate system)
        
        double Lsqrt_z = Math.sqrt( 4*Sinv_12*Sinv_12 
                                   + (Sinv_11 - Sinv_22)*(Sinv_11 - Sinv_22));

        double L1_z = (Sinv_11 + Sinv_22 - Lsqrt_z)*0.5;
        double L2_z = (Sinv_11 + Sinv_22 + Lsqrt_z)*0.5;

        double R_11_z, R_12_z, R_21_z, R_22_z;
        if (Math.abs(Sinv_12)<Double.MIN_NORMAL) { // already diagonalized
            R_11_z = 0;
            R_12_z = 1;
            R_21_z = 1;
            R_22_z = 0;
        } else {
            double N1inv_z = 1/Math.sqrt((Sinv_22-L1_z)*(Sinv_22-L1_z) + Sinv_12*Sinv_12);
            double N2inv_z = 1/Math.sqrt((Sinv_22-L2_z)*(Sinv_22-L2_z) + Sinv_12*Sinv_12);

            R_11_z = N1inv_z*(L1_z-Sinv_22);
            R_12_z = N2inv_z*(L2_z-Sinv_22);
            R_21_z = N1inv_z*(Sinv_12);
            R_22_z = N2inv_z*(Sinv_12);
        }

        // normal form of ellipses in x'-y' cooridnate system (completion of the square)
        double L1inv_z = 1/L1_z;
        double L2inv_z = 1/L2_z;

        double center_coefficient_1_z = (Sinv_13*R_11_z+Sinv_23*R_21_z)*L1inv_z;
        double center_coefficient_2_z = (Sinv_13*R_12_z+Sinv_23*R_22_z)*L2inv_z;

        double Rhs_tmp_xz = (Sinv_22*Sinv_13 - Sinv_12*Sinv_23)*Sinv_13;
        double Rhs_tmp_yz = (Sinv_11*Sinv_23 - Sinv_12*Sinv_13)*Sinv_23;

        double Rhs_num_z = Rhs_tmp_xz + Rhs_tmp_yz;
        double Rhs_den_z = Sinv_11*Sinv_22 - Sinv_12*Sinv_12;
        double Rhs_tmp2_z = Rhs_num_z/Rhs_den_z-Sinv_33;

        // obtain coordinates of edge points on x'-y' coordinate systems, and
        // convert them to x-y coodinate system

        double Lninv_z, R_1n_z, R_2n_z;
        if (L1_z<L2_z) { Lninv_z = L1inv_z; R_1n_z = R_11_z; R_2n_z = R_21_z; }
        else {           Lninv_z = L2inv_z; R_1n_z = R_12_z; R_2n_z = R_22_z; }

        double aspectRatio = Math.sqrt(Math.min(L1_z,L2_z)/Math.max(L1_z,L2_z));

        double[] params = new double[12];
        params[0]  = center_coefficient_1_z;
        params[1]  = center_coefficient_2_z;
        params[2]  = Rhs_tmp2_z;
        params[3]  = Lninv_z;
        params[4]  = R_1n_z;
        params[5]  = R_2n_z;
        params[6]  = R_11_z;
        params[7]  = R_12_z;
        params[8]  = R_21_z;
        params[9]  = R_22_z;
        params[10] = aspectRatio;
        //params[11] = zrmax

        return params;
    }

    // get x1,y1,x2,y2 of EllipseRoi
    public double[] getEdge(double zd, double[] params) {
        double center_coefficient_1_z = params[0];
        double center_coefficient_2_z = params[1];
        double Rhs_tmp2_z             = params[2];
        double Lninv_z                = params[3];
        double R_1n_z                 = params[4];
        double R_2n_z                 = params[5];
        double R_11_z                 = params[6];
        double R_12_z                 = params[7];
        double R_21_z                 = params[8];
        double R_22_z                 = params[9];

        double Rhs_z = zd*zd*Rhs_tmp2_z + nsigma2;
        double distance_tmp = Math.sqrt(Rhs_z*Lninv_z);
        double distance_1 = distance_tmp*R_1n_z;
        double distance_2 = distance_tmp*R_2n_z;
        double center2_1 = -zd*center_coefficient_1_z;
        double center2_2 = -zd*center_coefficient_2_z;
        double center_1 = center2_1*R_11_z + center2_2*R_12_z;
        double center_2 = center2_1*R_21_z + center2_2*R_22_z;

        double[] edge = new double[4];
        edge[0] = (center_1 + distance_1); // x_1 for z=z0, edge_11
        edge[1] = (center_2 + distance_2); // y_1 for z=z0, edge_12
        edge[2] = (center_1 - distance_1); // x_2 for z=z0, edge_21
        edge[3] = (center_2 - distance_2); // y_2 for z=z0, edge_22
        return edge;
    }

    public MyEllipseRoi getMyEllipseAt(int pos, int mode) {
        return getMyEllipseAt((double)pos, mode);
    }
    
    public MyEllipseRoi getMyEllipseAt(double pos, int mode) {
        if (!hasEllipseAt(pos,mode)) { return null; }
        double zd;
        double[] edge;
        switch (mode) {
            case PROJ_MODE:
                if (roi_proj!=null) { return roi_proj; } // already calculated
                zd = 0;
                edge = getEdge(zd,params_xy);
                roi_proj = new MyEllipseRoi(edge[0]+xc_disp,edge[1]+yc_disp,
                                            edge[2]+xc_disp,edge[3]+yc_disp,
                                            params_xy[10]);
                updatePosition(pos,mode);
                roi_proj.setName(name);
                roi_proj.setNameDisp(nameDisp);
                roi_proj.setFlagUseNameDisp(flagUseNameDisp);
                roi_proj.setFlagUseFixedName(flagUseFixedName);
                roi_proj.setStrokeColor(color);
                return roi_proj;
            case XY_MODE:
                if (z_last==pos && roi_xy!=null) { return roi_xy; } // already calculated
                zd = pos - zc_slice;
                edge = getEdge(zd,params_xy);
                roi_xy = new MyEllipseRoi(edge[0]+xc_disp,edge[1]+yc_disp,
                                            edge[2]+xc_disp,edge[3]+yc_disp,
                                            params_xy[10]);
                updatePosition(pos,mode);
                roi_xy.setName(name);
                roi_xy.setNameDisp(nameDisp);
                roi_xy.setFlagUseNameDisp(flagUseNameDisp);
                roi_xy.setFlagUseFixedName(flagUseFixedName);
                roi_xy.setStrokeColor(color);
                z_last = pos;
                return roi_xy;
            case XZ_MODE:
                if (y_last==pos && roi_xz!=null) { return roi_xz; } // already calculated
                zd = pos - yc_slice;
                edge = getEdge(zd,params_xz);
                roi_xz = new MyEllipseRoi(edge[0]+xc_disp,edge[1]+zc_disp,
                                            edge[2]+xc_disp,edge[3]+zc_disp,
                                            params_xz[10]);
                updatePosition(pos,mode);
                roi_xz.setName(name);
                roi_xz.setNameDisp(nameDisp);
                roi_xz.setFlagUseNameDisp(flagUseNameDisp);     
                roi_xz.setFlagUseFixedName(flagUseFixedName);     
                roi_xz.setStrokeColor(color);
                y_last = pos;
                return roi_xz;
            case YZ_MODE:
                if (x_last==pos && roi_yz!=null) { return roi_yz; } // already calculated
                zd = pos - xc_slice;
                edge = getEdge(zd,params_yz);
                roi_yz = new MyEllipseRoi(edge[1]+zc_disp,edge[0]+yc_disp,
                                            edge[3]+zc_disp,edge[2]+yc_disp,
                                            params_yz[10]);
                updatePosition(pos,mode);
                roi_yz.setName(name);
                roi_yz.setNameDisp(nameDisp);
                roi_yz.setFlagUseNameDisp(flagUseNameDisp);
                roi_yz.setFlagUseFixedName(flagUseFixedName);
                roi_yz.setStrokeColor(color);
                x_last = pos;
                return roi_yz;
            default:
                return null;
        }
    }

    public boolean hasEllipseAt(int pos, int mode) {
        return hasEllipseAt((double)pos, mode);
    }
    
    public boolean hasEllipseAt(double pos, int mode){
        if (pos < 1) { return false; } // to avoid non-positive position
        if ( modeCZT==1 || modeCZT==-2 || modeCZT==-3 || modeCZT==-4 ) { // no Z-slices.
            switch (mode) {
                case PROJ_MODE:
                    break;
                case XY_MODE:
                    return pos==1;
                case XZ_MODE: 
                    return false;
                case YZ_MODE:
                    return false;
            }
        }
        double zd, Rhs_z, Rhs_tmp2_z;
        switch (mode) {
            case PROJ_MODE:
                zd = 0;
                Rhs_tmp2_z = params_xy[2];
                break;
            case XY_MODE:
                zd = pos - zc_slice;
                Rhs_tmp2_z = params_xy[2];
                break;
            case XZ_MODE:
                zd = pos - yc_slice;
                Rhs_tmp2_z = params_xz[2];
                break;
            case YZ_MODE:
                zd = pos - xc_slice;
                Rhs_tmp2_z = params_yz[2];
                break;
            default:
                return false;
        }
        Rhs_z = zd*zd*Rhs_tmp2_z + nsigma2;
        return Rhs_z >= 0; // to avoid complex numbers
    }

    public void setName(String name) {
        this.name = name;
        this.nameDisp = name;
        if (roi_proj!=null) { roi_proj.setName(name); }
        if (roi_xy!=null)   {   roi_xy.setName(name); }
        if (roi_xz!=null)   {   roi_xz.setName(name); }
        if (roi_yz!=null)   {   roi_yz.setName(name); }        
    }

    public String getName() {
        if (flagUseNameDisp && nameDisp!=null && 
            !(flagUseFixedName && !NumberUtils.isNumber(name))) {
            return getNameDisp();
        } else {
            return getNameOrig();
        }
    }

    
    public void setNameDisp(String name) {        
        this.nameDisp = name;
        if (roi_proj!=null) { roi_proj.setNameDisp(name); }
        if (roi_xy!=null)   {   roi_xy.setNameDisp(name); }
        if (roi_xz!=null)   {   roi_xz.setNameDisp(name); }
        if (roi_yz!=null)   {   roi_yz.setNameDisp(name); }
    }
    
    public String getNameDisp() {
        return nameDisp;
    }
    
    public void setNameOrig(String name) {
        this.name = name;
        if (roi_proj!=null) { roi_proj.setNameOrig(name); }
        if (roi_xy!=null)   {   roi_xy.setNameOrig(name); }
        if (roi_xz!=null)   {   roi_xz.setNameOrig(name); }
        if (roi_yz!=null)   {   roi_yz.setNameOrig(name); }
    }
    
    public String getNameOrig() {
        return name;
    }

    public void setFlagUseNameDisp(boolean tf) {
        flagUseNameDisp = tf;
        if (roi_proj!=null) { roi_proj.setFlagUseNameDisp(tf); }
        if (roi_xy!=null)   {   roi_xy.setFlagUseNameDisp(tf); }
        if (roi_xz!=null)   {   roi_xz.setFlagUseNameDisp(tf); }
        if (roi_yz!=null)   {   roi_yz.setFlagUseNameDisp(tf); }
    }
    
    public boolean getFlagUseNameDisp() {
        return flagUseNameDisp;
    }
    
    public void setFlagUseFixedName(boolean tf) {
        flagUseFixedName = tf;
        if (roi_proj!=null) { roi_proj.setFlagUseFixedName(tf); }
        if (roi_xy!=null)   {   roi_xy.setFlagUseFixedName(tf); }
        if (roi_xz!=null)   {   roi_xz.setFlagUseFixedName(tf); }
        if (roi_yz!=null)   {   roi_yz.setFlagUseFixedName(tf); }
    }
    
    public boolean getFlagUseFixedName() {
        return flagUseFixedName;
    }
    

    public void setStrokeColor(Color c) {
        this.color = c;
        if (roi_proj!=null) { roi_proj.setStrokeColor(c); }
        if (roi_xy!=null)   {   roi_xy.setStrokeColor(c); }
        if (roi_xz!=null)   {   roi_xz.setStrokeColor(c); }
        if (roi_yz!=null)   {   roi_yz.setStrokeColor(c); }
    }

    public int getCPosition() {
        return channel;
    }

    public int getTPosition() {
        return frame;
    }

    public void setCPosition(int channel) {
        this.channel = channel;
        //this.updatePosition();
    }
    
    public void setTPosition(int frame) {
        this.frame = frame;
        //this.updatePosition();
    }

    public void setCTPosition(int channel, int frame) {
        this.channel = channel;
        this.frame = frame;
        //this.updatePosition();
    }   
    
    public void setModeCZT(int modeCZT) {
        this.modeCZT = modeCZT;

	switch (modeCZT) {
            case 4: // hyperstack (ztc)
            case 3: // hyperstack (zt)
            case 2: // hyperstack (zc)
            case 1: // hyperstack (tc)
            case -1: // z-only
                position = slice;
                break;
            case -2: // t-only
                position = frame;
                break;
            case -3: // c-only
            case -4: // single
                position = channel;
                break;
        }       
        //this.updatePosition();
    }
       
    
    public void updatePosition(double pos, int mode) {
        switch (mode) {
            case PROJ_MODE:
                switch (modeCZT) {
                    case 4: // hyperstack (ztc)
                        roi_proj.setPosition(channel,1,frame);
                        break;
                    case 3: // hyperstack (zt)
                        roi_proj.setPosition(frame);
                        break;
                    case 2: // hyperstack (zc)
                        roi_proj.setPosition(channel);
                        break;
                    case 1: // hyperstack (tc)
                        roi_proj.setPosition(channel,1,frame);
                        break;
                    case -1: // z-only
                        roi_proj.setPosition(0);
                        break;
                    case -2: // t-only
                        roi_proj.setPosition(frame);
                        break;
                    case -3: // c-only
                        roi_proj.setPosition(channel);
                        break;
                    case -4: // single
                        roi_proj.setPosition(channel);
                        break;
                }
                break;
            case XY_MODE:            	
                switch (modeCZT) {
                    case 4: // hyperstack (ztc)
                    case 3: // hyperstack (zt)
                    case 2: // hyperstack (zc)
                    case 1: // hyperstack (tc)
                        roi_xy.setPosition(channel,(int)pos,frame);
                        break;
                    case -1: // z-only
                        roi_xy.setPosition((int)pos);
                        break;
                    case -2: // t-only
                        roi_xy.setPosition(frame);
                        break;
                    case -3: // c-only
                    case -4: // single
                        roi_xy.setPosition(channel);
                        break;
                }
                break;
            case XZ_MODE:
                //roi_xz.setPosition(channel,(int)pos,frame);
                roi_xz.setPosition(0);
                break;
            case YZ_MODE:
                //roi_yz.setPosition(channel,(int)pos,frame);
                roi_yz.setPosition(0);
                break;
        }
    }

    @Override
    public int compareTo(Object obj){
        MyEllipsoid me2 = (MyEllipsoid)obj;
        if (name == null) { return -1; }
        String name_me2 = me2.getNameOrig();
        if (name_me2 == null) { return 1; }
        return name.compareTo(name_me2);
    }
    
}

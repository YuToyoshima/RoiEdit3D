import javax.swing.table.DefaultTableModel;
import javax.swing.event.TableModelListener;
import javax.swing.event.TableModelEvent;
import java.util.regex.Pattern;
import java.util.regex.Matcher;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.Stack;
import java.util.TreeMap;
import java.util.Vector;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import com.google.common.cache.*;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.MultimapBuilder;
import org.apache.commons.lang3.math.NumberUtils;
import org.apache.commons.lang3.ArrayUtils;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;


/* MyTableModel4 is based on MyTableModel3 and introducing mapF2R. */
public class MyTableModel4 extends DefaultTableModel{
    
    private static final long serialVersionUID = 4L;
    public Object oldobj;
    public Object newobj;
    public TreeMap<String,String[]> map = new TreeMap<String,String[]>();
    public transient Stack<History> undoStack = new Stack<History>();
    public transient Stack<History> redoStack = new Stack<History>();
    
    public double[][][] parallelShift = null;
    public boolean flagAlign = false;
    public boolean flagFlipX = false;
    public boolean flagFlipY = false;
    public boolean flagFlipZ = false;
    public boolean flagInterpZ = false;
    public double rotation = 0; // 0<rotation<2*pi
    public double zscale = 1;
    public int nxorig = 1;
    public int nyorig = 1;
    public int nzorig = 1;
    public int ncorig = 1;
    public int ntorig = 1;
    
    public transient LoadingCache<Integer,double[]> muCached;
    public transient LoadingCache<Integer,double[]> sigmaCached;
    public transient LoadingCache<Double,double[]> rotmat;
    
    public transient ArrayListMultimap<Integer,Integer> mapF2R;
    public transient boolean flagDirtyMapF2R = true;
    
    /* constructor */
    public MyTableModel4(){
        super();
        cacheConstructor();
    }
    
    public MyTableModel4(Object[][] data, Object[] columnNames){
        super(data, columnNames);
        cacheConstructor();
    }
    
    public MyTableModel4(Vector data, Vector columnNames){
        super(data, columnNames);
        cacheConstructor();
    }
    
    public void cacheConstructor() {
        muCached = CacheBuilder.newBuilder()
        .maximumSize(1000000)
        .build( new CacheLoader<Integer, double[]>() {
            public double[] load(Integer key) { return calcMuEffected(key); }
        });
        
        sigmaCached = CacheBuilder.newBuilder()
        .maximumSize(1000000)
        .build( new CacheLoader<Integer, double[]>() {
            public double[] load(Integer key) { return calcSigmaEffected(key); }
        });
        
        rotmat = CacheBuilder.newBuilder()
        .maximumSize(1000)
        .build( new CacheLoader<Double,double[]>() {
            public double[] load(Double key) { return calcRotMat(key); }
        });
    }
    
    
    /* functions should be override */
    @Override
    public Class getColumnClass(int col){
        return getValueAt(0, col).getClass();
    }
    
    @Override
    public boolean isCellEditable (int row, int col) {
        return true;
    }
    
    @Override
    // capture old object as oldobj
    public void setValueAt(Object obj, int row, int col) {
        addUndo("update",row);
        Vector v = (Vector) (dataVector.get(row));
        oldobj = v.get(col);
        newobj = obj;
        checkInputIsFrame(col);
        if (obj instanceof Double) {
            double[] tmpd;
            int[] idxColMu = getColumnIndexWarn("mu");
            int[] idxColSigma = getColumnIndexWarn("sigma");
            int tmpIdxMu = ArrayUtils.indexOf(idxColMu,col);
            int tmpIdxSigma = ArrayUtils.indexOf(idxColSigma,col);
            if (tmpIdxMu>=0) {
                tmpd = muCached.getUnchecked(row);
                tmpd[tmpIdxMu] = (Double)obj;
                // muCached.put(row,tmpd); // unnecessary because pointer (tmpd) is overwritten
                updateMuRaw(row,tmpd);
            } else if (tmpIdxSigma>=0) {
                tmpd = sigmaCached.getUnchecked(row);
                tmpd[tmpIdxSigma] = (Double)obj; // unnecessary because pointer (tmpd) is overwritten
                // sigmaCached.put(row,tmpd);
                updateSigmaRaw(row,tmpd);
            } else { super.setValueAt(obj,row,col); }
        } else { super.setValueAt(obj,row,col); }
    }
    
    
    /* setValuesAt for boolean */
    public void setValuesAt(boolean[][] obj, boolean[] flag, String str, boolean[] flags) {
        int[] row = convertBooleanToRows(flag);
        setValuesAt(obj,row,str,flags);
    }
    
    public void setValuesAt(boolean[][] obj, boolean[] flag, int[] col, boolean[] flags) {
        int[] row = convertBooleanToRows(flag);
        setValuesAt(obj,row,col,flags);
    }
    
    public void setValuesAt(boolean[][] obj, int[] row, String str, boolean[] flags) {
        int[] col = getColumnIndexWarn(str);
        setValuesAt(obj,row,col,flags);
    }
    
    public void setValuesAt(boolean[][] obj, int[] row, int[] col, boolean[] flags) {
        // flags = [fire,undo,redo]
        if (flags.length<=1 || flags[1]==true) { addUndo("update",row); }
        for (int r=0; r<row.length; r++) {
            Vector v = (Vector)(dataVector.get(row[r]));
            for (int c=0; c<col.length; c++) { v.set(col[c],obj[r][c]); }
        }
        checkInputIsFrame(col);
        if (flags.length>0 && flags[0]==true) { fireTableRowsUpdatedFromRows(row); }
    }
    
    
    /* setValuesAt for double */
    public void setValuesAt(double[][] obj, boolean[] flag, String str, boolean[] flags) {
        int[] row = convertBooleanToRows(flag);
        setValuesAt(obj,row,str,flags);
    }
    
    public void setValuesAt(double[][] obj, boolean[] flag, int[] col, boolean[] flags) {
        int[] row = convertBooleanToRows(flag);
        setValuesAt(obj,row,col,flags);
    }
    
    public void setValuesAt(double[][] obj, int[] row, String str, boolean[] flags) {
        int[] col = getColumnIndexWarn(str);
        setValuesAt(obj,row,col,flags);
    }
    
    public void setValuesAt(double[][] obj, int[] row, int[] col, boolean[] flags) {
        // flags = [fire,undo,redo]
        if (flags.length<=1 || flags[1]==true) { addUndo("update",row); }
        checkInputIsFrame(col);
        boolean flagUpdateMu = false;
        boolean flagUpdateSigma = false;        
        int[] idxColMu = getColumnIndexWarn("mu");
        int[] idxColSigma = getColumnIndexWarn("sigma");
        int[] tmpIdxMu = new int[col.length];
        int[] tmpIdxSigma = new int[col.length];
        double[][] muUpdate = null;
        double[][] sigmaUpdate = null;
        for (int c=0; c<col.length; c++) {
            tmpIdxMu[c] = ArrayUtils.indexOf(idxColMu,col[c]);
            tmpIdxSigma[c] = ArrayUtils.indexOf(idxColSigma,col[c]);
            if (tmpIdxMu[c]>=0) { flagUpdateMu = true; } 
            if (tmpIdxSigma[c]>=0) { flagUpdateSigma = true; }
        }
        if (flagUpdateMu) { muUpdate = new double[row.length][idxColMu.length]; }
        if (flagUpdateSigma) { sigmaUpdate = new double[row.length][idxColSigma.length]; }
        for (int r=0; r<row.length; r++) {
            Vector v = (Vector) (dataVector.get(row[r]));
            if (flagUpdateMu) { muUpdate[r] = muCached.getUnchecked(row[r]); }
            if (flagUpdateSigma) { sigmaUpdate[r] = sigmaCached.getUnchecked(row[r]); }
            for (int c=0; c<col.length; c++) {
                if (tmpIdxMu[c]>=0) { 
                    muUpdate[r][tmpIdxMu[c]] = obj[r][c]; 
                } else if (tmpIdxSigma[c]>=0) {                    
                    sigmaUpdate[r][tmpIdxSigma[c]] = obj[r][c];                    
                } else {
                    v.set(col[c],obj[r][c]);
                }
            }
        }
        if (flagUpdateMu) { 
            updateMuRaw(row,muUpdate); 
            muCached.invalidateAll();
        }
        if (flagUpdateSigma) { 
            updateSigmaRaw(row,sigmaUpdate); 
            sigmaCached.invalidateAll();
        }        
        if (flags.length>0 && flags[0]==true) { fireTableRowsUpdatedFromRows(row); }
    }
    
    
    /* setValuesAt for String */
    public void setValuesAt(String[][] obj, boolean[] flag, String str, boolean[] flags) {
        int[] row = convertBooleanToRows(flag);
        setValuesAt(obj,row,str,flags);
    }
    
    public void setValuesAt(String[][] obj, boolean[] flag, int[] col, boolean[] flags) {
        int[] row = convertBooleanToRows(flag);
        setValuesAt(obj,row,col,flags);
    }
    
    public void setValuesAt(String[][] obj, int[] row, String str, boolean[] flags) {
        int[] col = getColumnIndexWarn(str);
        setValuesAt(obj,row,col,flags);
    }
    
    public void setValuesAt(String[][] obj, int[] row, int[] col, boolean[] flags) {
        // flags = [fire,undo,redo]
        if (flags.length<=1 || flags[1]==true) { addUndo("update",row); }
        for (int r=0; r<row.length; r++) {
            Vector v = (Vector)(dataVector.get(row[r]));
            for (int c=0; c<col.length; c++) { v.set(col[c],obj[r][c]); }
        }
        checkInputIsFrame(col);
        if (flags.length>0 && flags[0]==true) { fireTableRowsUpdatedFromRows(row); }
    }
    
    
    
    /* setValuesInFrame series */
    public void setValuesInFrame(boolean[][] obj, int frame, String str, boolean[] flags) {
        int[] col = getColumnIndexWarn(str);
        setValuesInFrame(obj,frame,col,flags);
    }
    
    public void setValuesInFrame(boolean[][] obj, int frame, int[] col, boolean[] flags) {
        int[] row = convertFrameToRows(frame);
        setValuesAt(obj,row,col,flags);
    }
    
    public void setValuesInFrame(double[][] obj, int frame, String str, boolean[] flags) {
        int[] col = getColumnIndexWarn(str);
        setValuesInFrame(obj,frame,col,flags);
    }
    
    public void setValuesInFrame(double[][] obj, int frame, int[] col, boolean[] flags) {
        int[] row = convertFrameToRows(frame);
        setValuesAt(obj,row,col,flags);
    }
    
    public void setValuesInFrame(String[][] obj, int frame, String str, boolean[] flags) {
        int[] col = getColumnIndexWarn(str);
        setValuesInFrame(obj,frame,col,flags);
    }
    
    public void setValuesInFrame(String[][] obj, int frame, int[] col, boolean[] flags) {
        int[] row = convertFrameToRows(frame);
        setValuesAt(obj,row,col,flags);
    }
    
    
    
    /* setValuesAt Silent series*/
    public void setValuesAtSilent(boolean[][] obj, boolean[] flag, String str) {
        setValuesAt(obj,flag,str,new boolean[3]);
    }
    
    public void setValuesAtSilent(boolean[][] obj, boolean[] flag, int[] col) {
        setValuesAt(obj,flag,col,new boolean[3]);
    }
    
    public void setValuesAtSilent(boolean[][] obj, int[] row, String str) {
        setValuesAt(obj,row,str,new boolean[3]);
    }
    
    public void setValuesAtSilent(boolean[][] obj, int[] row, int[] col) {
        setValuesAt(obj,row,col,new boolean[3]);
    }
    
    public void setValuesAtSilent(double[][] obj, boolean[] flag, String str) {
        setValuesAt(obj,flag,str,new boolean[3]);
    }
    
    public void setValuesAtSilent(double[][] obj, boolean[] flag, int[] col) {
        setValuesAt(obj,flag,col,new boolean[3]);
    }
    
    public void setValuesAtSilent(double[][] obj, int[] row, String str) {
        setValuesAt(obj,row,str,new boolean[3]);
    }
    
    public void setValuesAtSilent(double[][] obj, int[] row, int[] col) {
        setValuesAt(obj,row,col,new boolean[3]);
    }
    
    public void setValuesAtSilent(String[][] obj, boolean[] flag, String str) {
        setValuesAt(obj,flag,str,new boolean[3]);
    }
    
    public void setValuesAtSilent(String[][] obj, boolean[] flag, int[] col) {
        setValuesAt(obj,flag,col,new boolean[3]);
    }
    
    public void setValuesAtSilent(String[][] obj, int[] row, String str) {
        setValuesAt(obj,row,str,new boolean[3]);
    }
    
    public void setValuesAtSilent(String[][] obj, int[] row, int[] col) {
        setValuesAt(obj,row,col,new boolean[3]);
    }
    
    
    /* setValuesInFrame as Silent series */
    public void setValuesInFrameSilent(boolean[][] obj, int frame, String str) {
        setValuesInFrame(obj,frame,str,new boolean[3]);
    }
    
    public void setValuesInFrameSilent(boolean[][] obj, int frame, int[] col) {
        setValuesInFrame(obj,frame,col,new boolean[3]);
    }
    
    public void setValuesInFrameSilent(double[][] obj, int frame, String str) {
        setValuesInFrame(obj,frame,str,new boolean[3]);
    }
    
    public void setValuesInFrameSilent(double[][] obj, int frame, int[] col) {
        setValuesInFrame(obj,frame,col,new boolean[3]);
    }
    
    public void setValuesInFrameSilent(String[][] obj, int frame, String str) {
        setValuesInFrame(obj,frame,str,new boolean[3]);
    }
    
    public void setValuesInFrameSilent(String[][] obj, int frame, int[] col) {
        setValuesInFrame(obj,frame,col,new boolean[3]);
    }
    
    
    
    @Override
    public Object getValueAt(int row, int col) {
        int[] idxColMu = getColumnIndexWarn("mu");
        int[] idxColSigma = getColumnIndexWarn("sigma");
        int tmpIdxMu = ArrayUtils.indexOf(idxColMu,col);
        int tmpIdxSigma = ArrayUtils.indexOf(idxColSigma,col);
        
        if (tmpIdxMu>=0) { return muCached.getUnchecked(row)[tmpIdxMu]; }
        else if (tmpIdxSigma>=0) { return sigmaCached.getUnchecked(row)[tmpIdxSigma]; }
        else { return ((Vector) dataVector.get(row)).get(col); }
    }
    
    /* getValuesAt */
    public Object getValuesAt(boolean[] flag, String str) {
        int[] row = convertBooleanToRows(flag);
        return getValuesAt(row,str);
    }
    
    public Object getValuesAt(boolean[] flag, int[] col) {
        int[] row = convertBooleanToRows(flag);
        return getValuesAt(row,col);
    }
    
    public Object getValuesAt(int[] row, String str) {
        int[] col = getColumnIndexWarn(str);
        return getValuesAt(row,col);
    }
    
    public Object getValuesAt(int[] row, int[] col) {
        if (row.length == 0) { return null; }
        Object obj = ((Vector) dataVector.get(0)).get(col[0]);
        if      (obj instanceof Boolean) {
            boolean[][] a = new boolean[row.length][col.length];
            for (int r=0; r<row.length; r++) {
                Vector v = (Vector)(dataVector.get(row[r]));
                for (int c=0; c<col.length; c++) { a[r][c] = (Boolean)(v.get(col[c])); }
            }
            return a;
        }
        else if (obj instanceof Double)  {
            double[][] a = new double[row.length][col.length];
            double[] tmpd;
            int[] idxColMu = getColumnIndexWarn("mu");
            int[] idxColSigma = getColumnIndexWarn("sigma");
            int[] tmpIdxMu = new int[col.length];
            int[] tmpIdxSigma = new int[col.length];
            for (int c=0; c<col.length; c++) {
                tmpIdxMu[c] = ArrayUtils.indexOf(idxColMu,col[c]);
                tmpIdxSigma[c] = ArrayUtils.indexOf(idxColSigma,col[c]);
            }
            for (int r=0; r<row.length; r++) {
                Vector v = (Vector)(dataVector.get(row[r]));
                for (int c=0; c<col.length; c++) {
                    if (tmpIdxMu[c]>=0) {
                        tmpd = muCached.getUnchecked(row[r]);
                        a[r][c] = tmpd[tmpIdxMu[c]];
                    }
                    else if (tmpIdxSigma[c]>=0) {
                        tmpd = sigmaCached.getUnchecked(row[r]);
                        a[r][c] = tmpd[tmpIdxSigma[c]];
                    } else {
                        a[r][c] = (Double)(v.get(col[c]));
                    }
                }
            }
            return a;
        }
        else if (obj instanceof String)  {
            String[][] a = new String[row.length][col.length];
            for (int r=0; r<row.length; r++) {
                Vector v = (Vector)(dataVector.get(row[r]));
                for (int c=0; c<col.length; c++) { a[r][c] = (String)(v.get(col[c])); }
            }
            return a;
        }
        return null;
    }
    
    public Object getValuesInFrame(int frame, String str) {
        int[] col = getColumnIndexWarn(str);
        return getValuesInFrame(frame,col);
    }
    
    public Object getValuesInFrame(int frame, int[] col) {
        return getValuesAt(convertFrameToRows(frame),col);
    }
    
    /* return row index of a roi specified by frame and name */
    public int getRowIndexOf(int frame, String name) {
        int cols[] = getColumnIndexWarn("Name");
        int[] rows = convertFrameToRows(frame);
        for (int cr=0; cr<rows.length; cr++) {
            Vector v = (Vector)(dataVector.get(rows[cr]));
            if (name.equals((String)(v.get(cols[0])))) { return rows[cr]; }
        }
        return -1;
    }
    
    /* return row index of a roi specified by frame and name */
    public int[] getRowIndexOf(int[] frames, String[] names) {
        int[] ret = new int[frames.length];
        for (int cq=0; cq<frames.length; cq++) {
            ret[cq] = getRowIndexOf(frames[cq],names[cq]);
        }
        return ret;
    }
    
    
    /* insert, remove, update */
    public void insertRowsAt(Object[][] obj, int row, boolean[] flags) {
        int[] rows = new int[obj.length];
        for (int i=0; i<rows.length; i++) { rows[i] = row+i; }
        insertRowsAt(obj, rows, flags);
    }
    
    public void insertRowsAt(Object[][] obj, int[] rows, boolean[] flags) {
        Vector<Vector> vv = new Vector<Vector>(obj.length);
        for (int i=0; i<obj.length; i++) { vv.add(new Vector(Arrays.asList(obj[i]))); }
        insertRowsAt(vv, rows, flags);
    }
    
    public void insertRowsAt(Vector<Vector> vv, int[] rows, boolean[] flags) {
        for (int i=0; i<rows.length; i++) {dataVector.add(rows[i],vv.get(i)); }
        // flags = [fire,undo,redo]
        if (flags.length>2 && flags[2]==true) { addRedo("insert",rows); }
        else if (flags.length<=1 || flags[1]==true) { addUndo("insert",rows); }
        invalidateCaches();
        flagDirtyMapF2R = true;
        
        // update mu & sigma
        updateMuRaw(rows);
        updateSigmaRaw(rows);
        muCached.invalidateAll();
        sigmaCached.invalidateAll();
        
        if (flags.length>0 && flags[0]==true) { fireTableRowsInsertedFromRows(rows); }
    }
    
    public void removeRowsAt(boolean[] flag, boolean[] flags) {
        int[] rows = convertBooleanToRows(flag);
        removeRowsAt(rows,flags);
    }
    
    public void removeRowsAt(int[] rows, boolean[] flags) {
        Arrays.sort(rows);
        // flags = [fire,undo,redo]
        if (flags.length>2 && flags[2]==true) { addRedo("remove",rows); }
        else if (flags.length<=1 || flags[1]==true) { addUndo("remove",rows); }
        for (int i=rows.length; i>0; i--) { dataVector.remove(rows[i-1]); }
        invalidateCaches();
        flagDirtyMapF2R = true;
        if (flags.length>0 && flags[0]==true) { fireTableRowsRemovedFromRows(rows); }
    }
    
    public void updateRowsAt(Object[][] obj, boolean[] flag, boolean[] flags) {
        int[] rows = convertBooleanToRows(flag);
        updateRowsAt(obj,rows,flags);
    }
    
    public void updateRowsAt(Object[][] obj, int[] rows, boolean[] flags) {
        Vector<Vector> vv = new Vector<Vector>(obj.length);
        for (int i=0; i<obj.length; i++) { vv.add(new Vector(Arrays.asList(obj[i]))); }
        updateRowsAt(vv, rows, flags);
    }
    
    public void updateRowsAt(Vector<Vector> vv, int[] rows, boolean[] flags) {
        // flags = [fire,undo,redo]
        if (flags.length>2 && flags[2]==true) { addRedo("update",rows); }
        else if (flags.length<=1 || flags[1]==true) { addUndo("update",rows); }
        for (int i=0; i<rows.length; i++) { dataVector.set(rows[i],vv.get(i)); }
        invalidateCaches();
        flagDirtyMapF2R = true;
        
        // update mu & sigma
        updateMuRaw(rows);
        updateSigmaRaw(rows);
        muCached.invalidateAll();
        sigmaCached.invalidateAll();
        
        if (flags.length>0 && flags[0]==true) { fireTableRowsUpdatedFromRows(rows); }
    }
    
    
    /* regular expression, strcmp */
    public boolean[][] strcmp(String str, boolean[] flag, String strcol) {
        int[] row = convertBooleanToRows(flag);
        return strcmp(str,row,strcol);
    }
    
    public boolean[][] strcmp(String str, boolean[] flag, int[] col) {
        int[] row = convertBooleanToRows(flag);
        return strcmp(str,row,col);
    }
    
    public boolean[][] strcmp(String str, int[] row, String strcol) {
        int[] col = getColumnIndexWarn(strcol);
        return strcmp(str,row,col);
    }
    
    public boolean[][] strcmp(String str, int[] row, int[] col) {
        if (row.length == 0) { return null; }
        boolean[][] a = new boolean[row.length][col.length];
        for (int r=0; r<row.length; r++) {
            Vector v = (Vector)(dataVector.get(row[r]));
            for (int c=0; c<col.length; c++) {
                a[r][c] = str.equals( (String)(v.get(col[c])) );
            }
        }
        return a;
    }
    
    
    /* regular expression, isNumberString */
    public boolean[][] isNumberString(boolean[] flag, String strcol) {
        int[] row = convertBooleanToRows(flag);
        return isNumberString(row,strcol);
    }
    
    public boolean[][] isNumberString(boolean[] flag, int[] col) {
        int[] row = convertBooleanToRows(flag);
        return isNumberString(row,col);
    }
    
    public boolean[][] isNumberString(int[] row, String strcol) {
        int[] col = getColumnIndexWarn(strcol);
        return isNumberString(row,col);
    }
    
    public boolean[][] isNumberString(int[] row, int[] col) {
        if (row.length == 0) { return null; }
        boolean[][] a = new boolean[row.length][col.length];
        for (int r=0; r<row.length; r++) {
            Vector v = (Vector)(dataVector.get(row[r]));
            for (int c=0; c<col.length; c++) {
                a[r][c] =  NumberUtils.isNumber( (String)(v.get(col[c])) );
            }
        }
        return a;
    }
    
    
    /* regular expression, regexpfind */
    public boolean[][] regexpfind(String regex, boolean[] flag, String strcol) {
        int[] row = convertBooleanToRows(flag);
        return regexpfind(regex,row,strcol);
    }
    
    public boolean[][] regexpfind(String regex, boolean[] flag, int[] col) {
        int[] row = convertBooleanToRows(flag);
        return regexpfind(regex,row,col);
    }
    
    public boolean[][] regexpfind(String regex, int[] row, String strcol) {
        int[] col = getColumnIndexWarn(strcol);
        return regexpfind(regex,row,col);
    }
    
    public boolean[][] regexpfind(String regex, int[] row, int[] col) {
        if (row.length == 0) { return null; }
        boolean[][] a = new boolean[row.length][col.length];
        Pattern p = Pattern.compile(regex);
        for (int r=0; r<row.length; r++) {
            Vector v = (Vector)(dataVector.get(row[r]));
            for (int c=0; c<col.length; c++) {
                a[r][c] =  p.matcher( (String)(v.get(col[c])) ).find();
            }
        }
        return a;
    }
    
    
    /* regular expression, regexptoken */
    public String[][] regexptoken(String regex, boolean[] flag, String strcol) {
        int[] row = convertBooleanToRows(flag);
        return regexptoken(regex,row,strcol);
    }
    
    public String[][] regexptoken(String regex, boolean[] flag, int[] col) {
        int[] row = convertBooleanToRows(flag);
        return regexptoken(regex,row,col);
    }
    
    public String[][] regexptoken(String regex, int[] row, String strcol) {
        int[] col = getColumnIndexWarn(strcol);
        return regexptoken(regex,row,col);
    }
    
    public String[][] regexptoken(String regex, int[] row, int[] col) {
        if (row.length == 0) { return null; }
        String[][] a = new String[row.length][col.length];
        Pattern p = Pattern.compile(regex);
        for (int r=0; r<row.length; r++) {
            Vector v = (Vector)(dataVector.get(row[r]));
            for (int c=0; c<col.length; c++) {
                Matcher m = p.matcher( (String)(v.get(col[c])) );
                if (m.find()) { a[r][c] = m.group(1); }
            }
        }
        return a;
    }
    
    
    /* save and load */
    public void save(String path) {
        try { saveKryo(path); }
        catch (Exception e) {
            e.printStackTrace();
            try { saveStd(path); }
            catch (Exception e2) { e2.printStackTrace(); }
        }
    }
    
    public static MyTableModel4 load(String path) {
        return load(path, false);
    }
    
    public static MyTableModel4 load(String path, boolean flagDebug) {
        MyTableModel4 mtm = null;
        try { mtm = loadKryo(path); }
        catch (Exception e) {
            if (flagDebug) {
                System.out.println("Error in loadKryo;");
                e.printStackTrace();
            }
            System.out.println("Switch to loadKryoRef;");
            try { mtm = loadKryoRef(path); }
            catch (Exception e2) {
                if (flagDebug) {
                    System.out.println("Error in loadKryoRef;");
                    e2.printStackTrace();
                }
                System.out.println("Switch to loadStd;");
                try { mtm = loadStd(path); }
                catch (Exception e3) { e3.printStackTrace(); }
            }     
        }
        return mtm;
    }
    
    public void saveStd(String path)
    throws java.io.FileNotFoundException, java.io.IOException {
        // remove listeners before save
        TableModelListener[] tmls = getTableModelListeners();
        for (int c=0; c<tmls.length; c++){
            removeTableModelListener(tmls[c]);
        }
        
        ObjectOutputStream oos =
        new ObjectOutputStream(
        new BufferedOutputStream(
        new FileOutputStream(path)));
        oos.writeObject(this);
        oos.close();
        
        // recover listeners
        for (int c=0; c<tmls.length; c++){
            addTableModelListener(tmls[c]);
        }
    }
    
    public static MyTableModel4 loadStd(String path) throws Exception {
        ObjectInputStream ois =
        new ObjectInputStream(
        new BufferedInputStream(
        new FileInputStream(path)));
        MyTableModel4 mtm;
        Object ro = ois.readObject();
        ois.close();
        if (ro instanceof MyTableModel4) {
            mtm = (MyTableModel4)ro;
        } else { // MyTableModel3
            MyTableModel3 mtm3 = (MyTableModel3)ro;
            mtm = new MyTableModel4(mtm3.getDataVector(),mtm3.getColumnIdentifiers());
            // copy field values fromMyTableModel3 object
            mtm.oldobj = mtm3.oldobj;
            mtm.newobj = mtm3.newobj;
            mtm.map = mtm3.map;
            mtm.parallelShift = mtm3.parallelShift;
            mtm.flagAlign = mtm3.flagAlign;
            mtm.flagFlipX = mtm3.flagFlipX;
            mtm.flagFlipY = mtm3.flagFlipY;
            mtm.flagFlipZ = mtm3.flagFlipZ;
            mtm.flagInterpZ = mtm3.flagInterpZ;
            mtm.rotation = mtm3.rotation;
            mtm.zscale = mtm3.zscale;
            mtm.nxorig = mtm3.nxorig;
            mtm.nyorig = mtm3.nyorig;
            mtm.nzorig = mtm3.nzorig;
            mtm.ncorig = mtm3.ncorig;
            mtm.ntorig = mtm3.ntorig;
        }
        mtm.cacheConstructor();
        mtm.flagDirtyMapF2R = true;
        return mtm;
    }
    
    public void saveKryo(String path) throws Exception {
        Kryo kryo = new Kryo(null);
        // remove listeners before save
        TableModelListener[] tmls = getTableModelListeners();
        for (int c=0; c<tmls.length; c++){
            removeTableModelListener(tmls[c]);
        }
        
        Output output = new Output(new FileOutputStream(path));
        kryo.writeObject(output,this);
        output.close();
        
        // recover listeners
        for (int c=0; c<tmls.length; c++){
            addTableModelListener(tmls[c]);
        }
    }
    
    public void saveKryoRef(String path) throws Exception {
        Kryo kryo = new Kryo();
        // remove listeners before save
        TableModelListener[] tmls = getTableModelListeners();
        for (int c=0; c<tmls.length; c++){
            removeTableModelListener(tmls[c]);
        }
        
        Output output = new Output(new FileOutputStream(path));
        kryo.writeObject(output,this);
        output.close();
        
        // recover listeners
        for (int c=0; c<tmls.length; c++){
            addTableModelListener(tmls[c]);
        }
    }
    
    public static MyTableModel4 loadKryo(String path) throws Exception {
        Kryo kryo = new Kryo(null);
        Input input = new Input(new FileInputStream(path));
        MyTableModel4 mtm = kryo.readObject(input,MyTableModel4.class);
        input.close();
        mtm.cacheConstructor();
        mtm.flagDirtyMapF2R = true;
        return mtm;
    }
    
    public static MyTableModel4 loadKryoRef(String path) throws Exception {
        Kryo kryo = new Kryo();
        Input input = new Input(new FileInputStream(path));
        MyTableModel4 mtm = kryo.readObject(input,MyTableModel4.class);
        input.close();
        mtm.cacheConstructor();
        mtm.flagDirtyMapF2R = true;
        return mtm;
    }
    
    
    
    /* undo and redo */
    public void addUndo(String strdo, int row) {
        int[] rows = {row};
        addUndo(strdo,rows);
    }
    
    public void addUndo(String strdo,int[] rows) {
        int[] idxColMu = getColumnIndexWarn("mu");
        int[] idxColSigma = getColumnIndexWarn("sigma");
        Vector<Vector> vv = new Vector<Vector>();
        for (int r=0; r<rows.length; r++) {
            //vv.add((Vector) dataVector.get(rows[r])); }
            //vv.add((Vector)(((Vector) dataVector.get(rows[r])).clone())); }
            Vector v = (Vector)(((Vector) dataVector.get(rows[r])).clone());
            double[] tmpmu = muCached.getUnchecked(rows[r]);
            double[] tmpsigma = sigmaCached.getUnchecked(rows[r]);
            for (int c=0; c<idxColMu.length; c++) { v.set(idxColMu[c],tmpmu[c]); }
            for (int c=0; c<idxColSigma.length; c++) { v.set(idxColSigma[c],tmpsigma[c]); }
            vv.add(v);
        }
        undoStack.push(new History(strdo,rows,vv));
        
        // redoStack should not be cleared when this function was called from redo
        StackTraceElement[] elms = Thread.currentThread().getStackTrace();
        for (int i=0; i<elms.length; i++) {
            if (elms[i].getMethodName().equals("redo")) { return; }
        }
        redoStack.clear();
    }
    
    public void addRedo(String strdo, int row) {
        int[] rows = {row};
        addRedo(strdo,rows);
    }
    
    public void addRedo(String strdo,int[] rows) {
        int[] idxColMu = getColumnIndexWarn("mu");
        int[] idxColSigma = getColumnIndexWarn("sigma");
        Vector<Vector> vv = new Vector<Vector>();
        for (int r=0; r<rows.length; r++) {
            //vv.add((Vector) dataVector.get(rows[r])); }
            //vv.add((Vector)(((Vector) dataVector.get(rows[r])).clone())); }
            Vector v = (Vector)(((Vector) dataVector.get(rows[r])).clone());
            double[] tmpmu = muCached.getUnchecked(rows[r]);
            double[] tmpsigma = sigmaCached.getUnchecked(rows[r]);
            for (int c=0; c<idxColMu.length; c++) { v.set(idxColMu[c],tmpmu[c]); }
            for (int c=0; c<idxColSigma.length; c++) { v.set(idxColSigma[c],tmpsigma[c]); }
            vv.add(v);
        }
        redoStack.push(new History(strdo,rows,vv));
    }
    
    public void undo(boolean flag_fire) {
        if (undoStack.empty()) { return; }
        boolean[] flags = new boolean[] {flag_fire,false,true};
        History h = undoStack.pop();
        if      (h.strdo.equals("update")) { updateRowsAt(h.data,h.rows,flags); }
        else if (h.strdo.equals("insert")) { removeRowsAt(h.rows,flags); }
        else if (h.strdo.equals("remove")) { insertRowsAt(h.data,h.rows,flags); }
    }
    
    public void redo(boolean flag_fire) {
        if (redoStack.empty()) { return; }
        boolean[] flags = new boolean[] {flag_fire,true,false};
        History h = redoStack.pop();
        if      (h.strdo.equals("update")) { updateRowsAt(h.data,h.rows,flags); }
        else if (h.strdo.equals("insert")) { removeRowsAt(h.rows,flags); }
        else if (h.strdo.equals("remove")) { insertRowsAt(h.data,h.rows,flags); }
    }
    
    
    /* convert char to string */
    public class AutoConverterCharToString implements TableModelListener {
        public void tableChanged(TableModelEvent e) {
            if (e.getType()==e.DELETE) { return; }
            int col = e.getColumn();
            int c0 = 0;
            int c1 = getColumnCount();
            ArrayList<Integer> colarray = new ArrayList<Integer>();
            if (col!=e.ALL_COLUMNS) { // column was specified
                c0 = col;
                c1 = col+1;
            }
            for (int i=c0; i<c1; i++) {
                Object obj = getValueAt(0,i);
                if (obj instanceof String || obj instanceof Character) {
                    colarray.add(i);
                }
            }
            if (colarray.size() == 0) { return; }
            int r0 = e.getFirstRow();
            int r1 = e.getLastRow();
            int maxrow = getRowCount()-1;
            if (r0 == e.HEADER_ROW) { // all data could be changed
                r0 = 0;
                r1 = maxrow;
            }
            if (r0 < 0) { r0 = 0; }
            if (r1 > maxrow) { r1 = maxrow; }
            for (int j=r0; j<=r1; j++) {
                Vector v = ((Vector) dataVector.get(j));
                for (int coltarget: colarray) {
                    Object obj = v.get(coltarget);
                    if (obj instanceof Character) {
                        v.set(coltarget,String.valueOf(obj));
                    }
                }
            }
        } // end of tableChanged function
    } // end of AutoConverterCharToString class
    
    public void addAutoConverterCharToString() { // add char_to_string converter
        TableModelListener[] tmls = getTableModelListeners();
        for (int i=0; i<tmls.length; i++) {
            if (tmls[i] instanceof AutoConverterCharToString) {return;} // do nothing
        }
        addTableModelListener( new AutoConverterCharToString() );
    }
    
    public void convertCharToString() {
        AutoConverterCharToString acct = new AutoConverterCharToString();
        TableModelEvent e = new TableModelEvent(this);
        acct.tableChanged(e);
    }
    
    public void convertCharToString(int firstRow, int lastRow) {
        AutoConverterCharToString acct = new AutoConverterCharToString();
        TableModelEvent e = new TableModelEvent(this,firstRow,lastRow);
        acct.tableChanged(e);
    }
    
    public void convertCharToString(int firstRow, int lastRow, int col) {
        AutoConverterCharToString acct = new AutoConverterCharToString();
        TableModelEvent e = new TableModelEvent(this,firstRow,lastRow,col);
        acct.tableChanged(e);
    }
    
    
    /* auxiliary, convert column names to index */
    public Vector getColumnIdentifiers() {
        return columnIdentifiers;
    }
    
    public int[] getColumnIndex(String str) {
        String[] arr;
        if (map.containsKey(str)) { arr = map.get(str); }
        else { arr = new String[1]; arr[0] = str; }
        int[] ret = new int[arr.length];
        Arrays.fill(ret,-1);
        for (int q=0; q<ret.length; q++) {
            for (int p=0; p<columnIdentifiers.size(); p++) {
                if (arr[q].equals(columnIdentifiers.get(p))) {
                    ret[q] = p;
                    break;
                }
            }
        }
        return ret;
    }
    
    public int[] getColumnIndexWarn(String str) {
        int[] ret = getColumnIndex(str);
        for (int p=0; p<ret.length; p++) {
            if (ret[p]==-1) {
                throw new IllegalArgumentException("input string is not found in column name");
            }
        }
        return ret;
    }
    
    public void setMap(String str, String[] arr) {
        map.put(str,arr);
    }
    
    public int[] convertBooleanToRows(boolean[] flag) {
        if (flag.length!=dataVector.size()) {
            System.err.println("Error: boolean length mismatch!");
            return null;
        }
        ArrayList<Integer> row = new ArrayList<Integer>(flag.length);
        for (int i=0; i<flag.length; i++) {
            if (flag[i]) row.add(i);
        }
        if (row.size()==0) return new int[0];
        return ArrayUtils.toPrimitive(row.toArray(new Integer[row.size()]));
    }
    
    public int[] convertFrameToRows(int frame) {
        if (flagDirtyMapF2R) { initMapF2R(); }
        int[] rows = ArrayUtils.toPrimitive(mapF2R.get(frame).toArray(new Integer[0]));
        return rows;
        
//         int[] tmpIdxColFrame = getColumnIndexWarn("frame");
//         int idxColFrame = tmpIdxColFrame[0];
//         ArrayList<Integer> rowarr = new ArrayList<Integer>();
//         for (int cr=0; cr<dataVector.size(); cr++) {
//             double tmpf = (Double)(((Vector)(dataVector.get(cr))).get(idxColFrame));
//             if (((int)tmpf)==frame) { rowarr.add(cr); }
//         }
//         int[] rows = ArrayUtils.toPrimitive(rowarr.toArray(new Integer[rowarr.size()]));
//         return rows;
    }
    
    public void initMapF2R() {
        if (!flagDirtyMapF2R) { return; }
        mapF2R =  ArrayListMultimap.create();
        int[] tmpIdxColFrame = getColumnIndexWarn("frame");
        int idxColFrame = tmpIdxColFrame[0];
        for (int cr=0; cr<dataVector.size(); cr++) {
            double tmpf = (Double)(((Vector)(dataVector.get(cr))).get(idxColFrame));
            mapF2R.put((int)tmpf,cr);
        }
        flagDirtyMapF2R = false;
    }
    
    public void checkInputIsFrame(int col) {
        checkInputIsFrame(new int[] {col});
    }
    
    public void checkInputIsFrame(int[] col){
        int[] idxColFrame = getColumnIndexWarn("frame");
        if (ArrayUtils.indexOf(col,idxColFrame[0])>=0) {
            flagDirtyMapF2R = true;
        }
    }
    
    
    /* fire table change */
    public void fireTableRowsUpdatedFromRows(int[] rows) {
        int[] minmaxrow = minmax(rows);
        fireTableRowsUpdated(minmaxrow[0],minmaxrow[1]);
    }
    
    public void fireTableRowsRemovedFromRows(int[] rows) {
        int[] minmaxrow = minmax(rows);
        fireTableRowsDeleted(minmaxrow[0],minmaxrow[1]);
    }
    
    public void fireTableRowsInsertedFromRows(int[] rows) {
        int[] minmaxrow = minmax(rows);
        fireTableRowsInserted(minmaxrow[0],minmaxrow[1]);
    }
    
    public int[] minmax(int[] rows) {
        int[] minmaxrow = {Integer.MAX_VALUE,Integer.MIN_VALUE};
        for (int i : rows) {
            if (minmaxrow[0]>i) minmaxrow[0] = i;
            if (minmaxrow[1]<i) minmaxrow[1] = i;
        }
        return minmaxrow;
    }
    
    
    
    /******** effector section ********/
    
    /* Calculation of rotation matrix and offset.
     *
     * cr = cos(rad);
     * sr = sin(rad);
     * zs = zscale; // or zs = 1 if flagInterpZ equals true (i.e. already corrected)
     * mu = [x,y,z]';
     * sigma = [s11,s12,s13;
     *          s12,s22,s23;
     *          s13,s23,s33];
     * offset = [ox;oy;oz];
     * scaling_mat = [1, 0,  0;
     *                0, 1,  0;
     *                0, 0, zs];
     * rotxmat = [1, 0,    0;
     *            0, cr, -zs;
     *            0, zs,  cr];
     * mu_rot = (inv(scaling_mat)*(rotxmat*(scaling_mat*(mu-offset))))+offset;
     *        = rotxpos*(mu-offset)+offset;
     * rotxpos = [1,     0,      0;
     *            0,    cr, -zs*sr;
     *            0, sr/zs,     cr];
     * sigma_rot = inv(scaling_mat)*(rotxmat*(scaling_mat*sigma*scaling_mat)*rotxmat')*inv(scaling_mat)
                 = rotxpos*sigma*rotxpos';
     *           = rotxsigma*sigma;
     *
     * mu_rot[0] = mu[0];
     * mu_rot[1] = cr   *mu[1] + (-sr*zs)*mu[2] + oy - cr*oy + oz*sr*zs;
     * mu_rot[2] = sr/zs*mu[1] +   cr    *mu[2] + oz - cr*oz - oy*sr/zs;
     *
     * mu[0] =   mu_rot[0];
     * mu[1] =   cr    *mu_rot[1] + (sr*zs)*mu_rot[2] + oy - cr*oy - oz*sr*zs;
     * mu[2] = (-sr/zs)*mu_rot[1] +     cr *mu_rot[2] + oz - cr*oz + oy*sr/zs;
     *
     * s11_rot = s11;
     * s12_rot =      cr*s12 + (-sr*zs)*s13;
     * s13_rot = (sr/zs)*s12 +       cr*s13;
     * s22_rot =        cr*cr *s22 + (-2*cr*sr*zs)*s23 + sr*sr*zs*zs*s33;
     * s23_rot =    (cr*sr/zs)*s22 + (cr*cr-sr*sr)*s23 + (-cr*sr*zs)*s33;
     * s33_rot = (sr*sr/zs/zs)*s22 +  (2*cr*sr/zs)*s23 +       cr*cr*s33;
     *
     * s11 = s11_rot;
     * s12 = cr*s12_rot + (sr*zs)*s13_rot;
     * s13 = (-sr/zs)*s12_rot + cr*s13_rot;
     * s22 =        cr*cr *s22_rot +  (2*cr*sr*zs)*s23_rot + (sr*sr*zs*zs)*s33_rot
     * s23 =   (-cr*sr/zs)*s22_rot + (cr*cr-sr*sr)*s23_rot +    (cr*sr*zs)*s33_rot
     * s33 = (sr*sr/zs/zs)*s22_rot + (-2*cr*sr/zs)*s23_rot +       (cr*cr)*s33_rot
     *
     * ret = [oy - cr*oy + oz*sr*zs,
     *        oz - cr*oz - oy*sr/zs, // offsets
     *        cr,
     *       -sr*zs,
     *        sr/zs,
     *        cr, // rotxpos
     *        cr*cr,
     *       -2*cr*sr*zs,
     *        sr*sr*zs*zs,
     *        cr*sr/zs,
     *        cr*cr-sr*sr,
     *       -cr*sr*zs,
     *        sr*sr/zs/zs,
     *        2*cr*sr/zs,
     *        cr*cr, // rotxsigma
     *        oy - cr*oy - oz*sr*zs,
     *        oz - cr*oz + oy*sr/zs, // offsets for inverse
     *        cr,
     *        sr*zs,
     *       -sr/zs,
     *        cr, // rotxpos for inverse
     *        cr*cr,
     *        2*cr*sr*zs,
     *        sr*sr*zs*zs,
     *       -cr*sr/zs,
     *        cr*cr-sr*sr,
     *        cr*sr*zs,
     *        sr*sr/zs/zs,
     *       -2*cr*sr/zs,
     *        cr*cr]; // rotxsigma for inverse
     *
     * mu_rot[0] = mu[0]
     * mu_rot[1] = ret[0] + ret[2]*mu[1] + ret[3]*mu[2];
     * mu_rot[2] = ret[1] + ret[4]*mu[1] + ret[5]*mu[2];
     *
     * s11_rot = s11;
     * s12_rot = ret[ 2]*s12 + ret[ 3]*s13;
     * s13_rot = ret[ 4]*s12 + ret[ 5]*s13;
     * s22_rot = ret[ 6]*s22 + ret[ 7]*s23 + ret[ 8]*s33;
     * s23_rot = ret[ 9]*s22 + ret[10]*s23 + ret[11]*s33;
     * s33_rot = ret[12]*s22 + ret[13]*s23 + ret[14]*s33;
     *
     * mu[0] = mu_rot[0]
     * mu[1] = ret[15] + ret[17]*mu_rot[1] + ret[18]*mu_rot[2];
     * mu[2] = ret[16] + ret[19]*mu_rot[1] + ret[20]*mu_rot[2];
     *
     * s11 = s11_rot;
     * s12 = ret[17]*s12_rot + ret[18]*s13_rot;
     * s13 = ret[19]*s12_rot + ret[20]*s13_rot;
     * s22 = ret[21]*s22_rot + ret[22]*s23_rot + ret[23]*s33_rot;
     * s23 = ret[24]*s22_rot + ret[25]*s23_rot + ret[26]*s33_rot;
     * s33 = ret[27]*s22_rot + ret[28]*s23_rot + ret[29]*s33_rot;
     *
     */
    public double[] calcRotMat(double radian) {
        double cr = Math.cos(radian);
        double sr = Math.sin(radian);
        double zs = zscale;
        double nz = nzorig;
        if (flagInterpZ) { zs = 1; nz = (nzorig-1)*zscale+1; }
        double oy = nyorig*0.5+0.5;
        double oz = nz    *0.5+0.5;
        
        double[] ret = {
            oy - cr*oy + oz*sr*zs,
            oz - cr*oz - oy*sr/zs, // offsets
            cr,
            -sr*zs,
            sr/zs,
            cr, // rotxpos
            cr*cr,
            -2*cr*sr*zs,
            sr*sr*zs*zs,
            cr*sr/zs,
            cr*cr-sr*sr,
            -cr*sr*zs,
            sr*sr/zs/zs,
            2*cr*sr/zs,
            cr*cr, // rotxsigma
            oy - cr*oy - oz*sr*zs,
            oz - cr*oz + oy*sr/zs, // offsets for inverse
            cr,
            sr*zs,
            -sr/zs,
            cr, // rotxpos for inverse
            cr*cr,
            2*cr*sr*zs,
            sr*sr*zs*zs,
            -cr*sr/zs,
            cr*cr-sr*sr,
            cr*sr*zs,
            sr*sr/zs/zs,
            -2*cr*sr/zs,
            cr*cr}; // rotxsigma for inverse
        return ret;
    }
    
    /* apply effect to mu */
    public double[] calcMuEffected(int key) {
        int[] idxColMu = getColumnIndexWarn("mu");
        int idxColFrame = getColumnIndexWarn("frame")[0];
        double[] mu = new double[3];
        Vector v = (Vector) (dataVector.get(key));
        for (int p=0; p<3; p++) { mu[p] = (Double)(v.get(idxColMu[p])); } // load
        if (flagAlign) { // calc shift amount for alignment
            int fq = ((Number)(v.get(idxColFrame))).intValue();
            int zq_floor  = Math.min(Math.max((int)(mu[2]),1),nzorig);
            int zq_floor2 = Math.min(Math.max(zq_floor+1,  1),nzorig);
            double zq_dist = mu[2] - zq_floor;
            mu[0] += (1-zq_dist)*parallelShift[0][zq_floor -1][fq-1]
            +  zq_dist *parallelShift[0][zq_floor2-1][fq-1];
            mu[1] += (1-zq_dist)*parallelShift[1][zq_floor -1][fq-1]
            +  zq_dist *parallelShift[1][zq_floor2-1][fq-1];
        }
        if (flagFlipX) { mu[0] = nxorig - mu[0] + 1; } // flipX
        if (flagFlipY) { mu[1] = nyorig - mu[1] + 1; } // flipY
        if (flagFlipZ) { mu[2] = nzorig - mu[2] + 1; } // flipZ
        if (flagInterpZ) { mu[2] = (mu[2]-1)*zscale + 1; } // interpZ
        if (rotation!=0) { // rotation with x axis
            double[] ret = rotmat.getUnchecked(rotation);
            double tmpmu1 = mu[1];
            double tmpmu2 = mu[2];
            mu[1] = ret[0] + ret[2]*tmpmu1 + ret[3]*tmpmu2;
            mu[2] = ret[1] + ret[4]*tmpmu1 + ret[5]*tmpmu2;
        }
        return mu;
    }
    
    /* remove effect from affected mu */
    public double[] calcMuRaw(int row, double[] mu) {                
        int idxColFrame = getColumnIndexWarn("frame")[0];
        if (rotation!=0) {
            double[] ret = rotmat.getUnchecked(rotation);
            double tmpmu1 = mu[1];
            double tmpmu2 = mu[2];
            mu[1] = ret[15] + ret[17]*tmpmu1 + ret[18]*tmpmu2;
            mu[2] = ret[16] + ret[19]*tmpmu1 + ret[20]*tmpmu2;
        }
        if (flagInterpZ) { mu[2] = (mu[2]-1)/zscale + 1; } // interpZ
        if (flagFlipZ) { mu[2] = nzorig - mu[2] + 1; } // flipZ
        if (flagFlipY) { mu[1] = nyorig - mu[1] + 1; } // flipY
        if (flagFlipX) { mu[0] = nxorig - mu[0] + 1; } // flipX
        if (flagAlign) { // calc shift amount for alignment
            Vector v = (Vector) (dataVector.get(row));
            int fq = ((Number)(v.get(idxColFrame))).intValue();
            /* replaced because of integrity between calcMuEffected
            int zq_floor = (int)(mu[2]);
            int zq_floor2 = zq_floor+1;
            if (zq_floor <1) { zq_floor  = 1; }
            if (zq_floor2<1) { zq_floor2 = 1; }
            if (zq_floor >nzorig) { zq_floor  = nzorig; }
            if (zq_floor2>nzorig) { zq_floor2 = nzorig; }
             */
            int zq_floor  = Math.min(Math.max((int)(mu[2]),1),nzorig);
            int zq_floor2 = Math.min(Math.max(zq_floor+1,  1),nzorig);
            double zq_dist = mu[2] - zq_floor;
            mu[0] -= (1-zq_dist)*parallelShift[0][zq_floor -1][fq-1]
            +  zq_dist *parallelShift[0][zq_floor2-1][fq-1];
            mu[1] -= (1-zq_dist)*parallelShift[1][zq_floor -1][fq-1]
            +  zq_dist *parallelShift[1][zq_floor2-1][fq-1];
        }
        return mu;
    }
    
    /* apply effect to sigma */
    public double[] calcSigmaEffected(int key) {
        int[] idxColSigma = getColumnIndexWarn("sigma");
        double[] sigma = new double[6];
        Vector v = (Vector) (dataVector.get(key));
        for (int p=0; p<6; p++) { sigma[p] = (Double)(v.get(idxColSigma[p])); } // load
        if (flagFlipX) { sigma[1]*=-1; sigma[2]*=-1; } // flipX changes S_12 & S_13
        if (flagFlipY) { sigma[1]*=-1; sigma[4]*=-1; } // flipY changes S_12 & S_23
        if (flagFlipZ) { sigma[2]*=-1; sigma[4]*=-1; } // flipZ changes S_13 % S_23
        if (flagInterpZ) {  // interpZ
            sigma[2] *= zscale; // S_13
            sigma[4] *= zscale; // S_23
            sigma[5] *= zscale*zscale; // S_33
        }
        if (rotation!=0) { // rotation with x axis
            double[] ret = rotmat.getUnchecked(rotation);
            double S_12 = sigma[1];
            double S_13 = sigma[2];
            double S_22 = sigma[3];
            double S_23 = sigma[4];
            double S_33 = sigma[5];
            sigma[1] = ret[ 2]*S_12 + ret[ 3]*S_13;
            sigma[2] = ret[ 4]*S_12 + ret[ 5]*S_13;
            sigma[3] = ret[ 6]*S_22 + ret[ 7]*S_23 + ret[ 8]*S_33;
            sigma[4] = ret[ 9]*S_22 + ret[10]*S_23 + ret[11]*S_33;
            sigma[5] = ret[12]*S_22 + ret[13]*S_23 + ret[14]*S_33;
        }
        return sigma;
    }
    
    /* remove effect from affected sigma */
    public double[] calcSigmaRaw(double[] sigma) {
        if (rotation!=0) {
            double[] ret = rotmat.getUnchecked(rotation);
            double S_12_rot = sigma[1];
            double S_13_rot = sigma[2];
            double S_22_rot = sigma[3];
            double S_23_rot = sigma[4];
            double S_33_rot = sigma[5];
            sigma[1] = ret[17]*S_12_rot + ret[18]*S_13_rot;
            sigma[2] = ret[19]*S_12_rot + ret[20]*S_13_rot;
            sigma[3] = ret[21]*S_22_rot + ret[22]*S_23_rot + ret[23]*S_33_rot;
            sigma[4] = ret[24]*S_22_rot + ret[25]*S_23_rot + ret[26]*S_33_rot;
            sigma[5] = ret[27]*S_22_rot + ret[28]*S_23_rot + ret[29]*S_33_rot;
        }
        if (flagInterpZ) {  // interpZ
            double zscale_inv = 1/zscale;
            sigma[2] *= zscale_inv; // S_13
            sigma[4] *= zscale_inv; // S_23
            sigma[5] *= zscale_inv*zscale_inv; // S_33
        }
        if (flagFlipX) { sigma[1]*=-1; sigma[2]*=-1; } // flipX changes S_12 & S_13
        if (flagFlipY) { sigma[1]*=-1; sigma[4]*=-1; } // flipY changes S_12 & S_23
        if (flagFlipZ) { sigma[2]*=-1; sigma[4]*=-1; } // flipZ changes S_13 % S_23
        return sigma;
    }
    
    public void setMuRaw(double[] mu, int row) {
        int[] idxColMu = getColumnIndexWarn("mu");
        Vector v = (Vector) (dataVector.get(row));
        for (int p=0; p<idxColMu.length; p++) { v.set(idxColMu[p],mu[p]); }
    }
    
    public void setSigmaRaw(double[] sigma, int row) {
        int[] idxColSigma = getColumnIndexWarn("sigma");
        Vector v = (Vector)(dataVector.get(row));
        for (int p=0; p<idxColSigma.length; p++) { v.set(idxColSigma[p],sigma[p]); }
    }
    
    public void updateMuRaw(int[] rows, double[][] muEffected) {
        for (int p=0; p<rows.length; p++) { updateMuRaw(rows[p],muEffected[p]); }
    }
    
    public void updateSigmaRaw(int[] rows, double[][] sigmaEffected) {
        for (int p=0; p<rows.length; p++) { updateSigmaRaw(rows[p],sigmaEffected[p]); }
    }
    
    public void updateMuRaw(int row, double[] muEffected) {
        setMuRaw(calcMuRaw(row,muEffected),row);        
    }
    
    public void updateSigmaRaw(int row, double[] sigmaEffected) {
        setSigmaRaw(calcSigmaRaw(sigmaEffected),row);
    }
    
    /* updateMuRaw and updateSigmaRaw without muEffected or sigmaEffected; 
     dataVector is considered as muEffected or sigmaEffected */
    public void updateMuRaw(int[] rows) {
        for (int p=0; p<rows.length; p++) { updateMuRaw(rows[p]); }
    }
    
    public void updateSigmaRaw(int[] rows) {
        for (int p=0; p<rows.length; p++) { updateSigmaRaw(rows[p]); }
    }
    
    public void updateMuRaw(int row) {
        int[] idxColMu = getColumnIndexWarn("mu");
        double[] muEffected = new double[idxColMu.length];
        Vector v = (Vector) (dataVector.get(row));
        for (int p=0; p<idxColMu.length; p++) { 
            muEffected[p] = (Double)(v.get(idxColMu[p])); 
        }
        setMuRaw(calcMuRaw(row,muEffected),row);        
    }
    
    public void updateSigmaRaw(int row) {
        int[] idxColSigma = getColumnIndexWarn("sigma");
        double[] sigmaEffected = new double[idxColSigma.length];
        Vector v = (Vector) (dataVector.get(row));
        for (int p=0; p<idxColSigma.length; p++) { 
            sigmaEffected[p] = (Double)(v.get(idxColSigma[p])); 
        }
        setSigmaRaw(calcSigmaRaw(sigmaEffected),row);
    }
    
    
    /* set parallelShift */
    public void setParallelShift(double[][][] in) {
        parallelShift = in;
        invalidateCaches();
    }
    
    /* set Flag for alignment */
    public void setAlign(boolean flag) {
        flagAlign = flag;
        invalidateCaches();
    }
    
    /* set FlipX */
    public void setFlipX(boolean flag) {
        flagFlipX = flag;
        invalidateCaches();
    }
    
    /* set FlipY */
    public void setFlipY(boolean flag) {
        flagFlipY = flag;
        invalidateCaches();
    }
    
    /* set FlipX */
    public void setFlipZ(boolean flag) {
        flagFlipZ = flag;
        invalidateCaches();
    }
    
    /* set interporation for Z */
    public void setInterpZ(boolean flag) {
        flagInterpZ = flag;
        invalidateCaches();
    }
    
    /* toggle Flag for alignment */
    public void toggleAlign() { setAlign(!flagAlign); }
    
    /* toggle FlipX */
    public void toggleFlipX() { setFlipX(!flagFlipX); }
    
    /* toggle FlipY */
    public void toggleFlipY() { setFlipY(!flagFlipY); }
    
    /* toggle FlipZ */
    public void toggleFlipZ() { setFlipZ(!flagFlipZ); }
    
    /* toggle interporation for Z */
    public void toggleInterpZ() { setInterpZ(!flagInterpZ); }
    
    /* set rotation (by radius unit) */
    public void setRotation(double radian) {
        if (0<=radian && radian<2*Math.PI) { rotation = radian; }
        else { rotation = (((radian/(2*Math.PI)%1)+1)%1)*2*Math.PI; }
        invalidateCaches();
    }
    
    /* get rotation by raduis unit */
    public double getRotation() {
        return rotation;
    }
    
    public void invalidateCaches() {
        muCached.invalidateAll();
        sigmaCached.invalidateAll();
        rotmat.invalidateAll();
    }
    
}



class History {
    public String strdo;
    public int[] rows;
    public Vector<Vector> data;
    public History(String strdo, int[] rows, Vector<Vector> data) {
        this.strdo = strdo;
        this.rows = rows;
        this.data = data;
    }
}

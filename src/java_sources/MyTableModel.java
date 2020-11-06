import javax.swing.table.DefaultTableModel;
import javax.swing.event.TableModelListener;
import javax.swing.event.TableModelEvent;
import java.util.Vector;
import java.lang.reflect.Array;
import java.util.regex.Pattern;
import java.util.regex.Matcher;
import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;
import java.util.Collections;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import org.apache.commons.lang3.math.NumberUtils;
import org.apache.commons.lang3.ArrayUtils;

public class MyTableModel extends DefaultTableModel{
    
    //private static final long serialVersionUID = -5815851557796090351L;
    private static final long serialVersionUID =  6994326089187124206L;
    public Object oldobj;
    public Object newobj;
    
    public MyTableModel(){
        super();
    }
    
    public MyTableModel(Object[][] data, Object[] columnNames){
        super(data, columnNames);
    }
    
    public MyTableModel(Vector data, Vector columnNames){
        super(data, columnNames);
    }
    
    @Override
    public Class getColumnClass(int col){
        return getValueAt(0, col).getClass();
    }
    
    @Override
    public boolean isCellEditable (int row, int col) {
        return true;
    }
    
    public Vector getColumnIdentifiers() {
        return columnIdentifiers;
    }
    
    public int getColumnIndex(String str) {
        for (int p=0; p<columnIdentifiers.size(); p++) {
            if (str.equals(columnIdentifiers.get(p))) { return p; }
        }
        return -1;
    }
    
    @Override
    // capture old object as oldobj
    public void setValueAt(Object obj, int row, int col) {
        Vector v = (Vector) (dataVector.get(row));
        oldobj = v.get(col);
        newobj = obj;
        super.setValueAt(obj,row,col);
    }
    
    public void fireTableRowsUpdatedFromRows(int[] row) {
        int minrow = Integer.MAX_VALUE;
        int maxrow = Integer.MIN_VALUE;
        for (int i : row) {
            if (minrow>i) minrow = i;
            if (maxrow<i) maxrow = i;
        }
        fireTableRowsUpdated(minrow,maxrow);
    }
        
    
    // setValuesAt, generic
    public<T> void setValuesAt(T[] obj, int[] row, int[] column) {
        setValuesAtSilent(obj,row,column);
        fireTableRowsUpdatedFromRows(row);
    }
    
    public<T> void setValuesAt(T[] obj, int[] row, int column) {
        setValuesAtSilent(obj,row,column);
        fireTableRowsUpdatedFromRows(row);
    }
    
    public<T> void setValuesAtSilent(T[] obj, int[] row, int[] column) {        
        for (int i=0; i<row.length; i++) {
            ((Vector) dataVector.get(row[i])).set(column[i],obj[i]);
        }
    }
    
    public<T> void setValuesAtSilent(T[] obj, int[] row, int column) {
        for (int i=0; i<row.length; i++) {
            ((Vector) dataVector.get(row[i])).set(column,obj[i]);
        }
    }
    
    // setValuesAt, for double
    public void setValuesAt(double[] obj, int[] row, int[] column) {
        setValuesAtSilent(obj,row,column);
        fireTableRowsUpdatedFromRows(row);
    }
    
    public void setValuesAt(double[] obj, int[] row, int column) {
        setValuesAtSilent(obj,row,column);
        fireTableRowsUpdatedFromRows(row);
    }
    
    public void setValuesAtSilent(double[] obj, int[] row, int[] column) {        
        for (int i=0; i<row.length; i++) {
            ((Vector) dataVector.get(row[i])).set(column[i],obj[i]);
        }
    }
    
    public void setValuesAtSilent(double[] obj, int[] row, int column) {
        for (int i=0; i<row.length; i++) {
            ((Vector) dataVector.get(row[i])).set(column,obj[i]);
        }
    }
    
    // setValuesAt, for boolean
    public void setValuesAt(boolean[] obj, int[] row, int[] column) {
        setValuesAtSilent(obj,row,column);
        fireTableRowsUpdatedFromRows(row);
    }
    
    public void setValuesAt(boolean[] obj, int[] row, int column) {
        setValuesAtSilent(obj,row,column);
        fireTableRowsUpdatedFromRows(row);
    }
    
    public void setValuesAtSilent(boolean[] obj, int[] row, int[] column) {        
        for (int i=0; i<row.length; i++) {
            ((Vector) dataVector.get(row[i])).set(column[i],obj[i]);
        }
    }
    
    public void setValuesAtSilent(boolean[] obj, int[] row, int column) {
        for (int i=0; i<row.length; i++) {
            ((Vector) dataVector.get(row[i])).set(column,obj[i]);
        }
    }
    
    
    
    // setValuesAt, generic, boolean
    public<T> void setValuesAt(T[] obj, boolean[] flag, int column) {
        int[] row = convertBooleanToRows(flag);
        setValuesAtSilent(obj,row,column);
        fireTableRowsUpdatedFromRows(row);
    }
    
    public<T> void setValuesAtSilent(T[] obj, boolean[] flag, int column) {
        int[] row = convertBooleanToRows(flag);
        setValuesAtSilent(obj,row,column);
    }
    
    // setValuesAt, for double, boolean
    public void setValuesAt(double[] obj, boolean[] flag, int column) {
        int[] row = convertBooleanToRows(flag);
        setValuesAtSilent(obj,row,column);
        fireTableRowsUpdatedFromRows(row);
    }
    
    public void setValuesAtSilent(double[] obj, boolean[] flag, int column) {
        int[] row = convertBooleanToRows(flag);
        setValuesAtSilent(obj,row,column);
    }
    
    // setValuesAt, for boolean, boolean
    public void setValuesAt(boolean[] obj, boolean[] flag, int column) {
        int[] row = convertBooleanToRows(flag);
        setValuesAtSilent(obj,row,column);
        fireTableRowsUpdatedFromRows(row);
    }
    
    public void setValuesAtSilent(boolean[] obj, boolean[] flag, int column) {
        int[] row = convertBooleanToRows(flag);
        setValuesAtSilent(obj,row,column);
    }
    
    
    
    
    public void insertRowsSilent(Object[] obj, int row) {
        Vector v = new Vector(Arrays.asList(obj));
        dataVector.add(row,v);
    }
    
    public void insertRowsSilent(Object[][] obj, int row) {
        Vector<Vector> vv = new Vector<Vector>(obj.length);
        for (int i=0; i<obj.length; i++) {
            vv.add(new Vector(Arrays.asList(obj[i])));
        }
        dataVector.addAll(row,vv);
    }
    
    public void removeRowsSilent(int row) {
        dataVector.remove(row);
    }
    
    public void removeRowsSilent(int[] rows) {
        Arrays.sort(rows);
        for (int i=rows.length; i>0; i--) {
            removeRowsSilent(rows[i-1]);
        }
    }
    
    public void updateRowsSilent(Object[] obj, int row) {
        Vector v = new Vector(Arrays.asList(obj));
        dataVector.set(row,v);
    }
    
    public void updateRowsSilent(Object[][] obj, int[] rows) {
        for (int i=0; i<rows.length; i++) {
            updateRowsSilent(obj[i],rows[i]);
        }
    }
    
    
    public boolean[] strcmp(String str, int column) {
        boolean[] flag = new boolean[dataVector.size()];
        for (int i=0; i<dataVector.size(); i++) {
            flag[i] = strcmp(str,i,column);
        }
        return flag;
    }
    
    public boolean[] strcmp(String str, int[] row, int column) {
        if (row.length == 0) { return null; }
        boolean[] flag = new boolean[row.length];
        for (int i=0; i<row.length; i++) {
            flag[i] = strcmp(str,row[i],column);
        }
        return flag;
    }
    
    public boolean[] strcmp(String str, boolean[] flag, int column) {
        int[] row = convertBooleanToRows(flag);
        return strcmp(str,row,column);
    }
    
    public boolean strcmp(String str, int row, int column) {          
        return str.equals(((Vector) dataVector.get(row)).get(column));
    }
    
    
    public boolean[] isNumberString(int column) {
        boolean[] flag = new boolean[dataVector.size()];
        for (int i=0; i<dataVector.size(); i++) {
            flag[i] = isNumberString(i,column);
        }
        return flag;
    }
    
    public boolean[] isNumberString(int[] row, int column) {
        if (row.length == 0) { return null; }
        boolean[] flag = new boolean[row.length];
        for (int i=0; i<row.length; i++) {
            flag[i] = isNumberString(row[i],column);
        }
        return flag;
    }
    
    public boolean[] isNumberString(boolean[] flag, int column) {
        int[] row = convertBooleanToRows(flag);
        return isNumberString(row,column);
    }
    
    public boolean isNumberString(int row, int column) {
        return NumberUtils.isNumber(
        (String)(((Vector) dataVector.get(row)).get(column)) );
    }

    
    public boolean[] regexpfind(String regex, int column) {
        Pattern p = Pattern.compile(regex);
        boolean[] flag = new boolean[dataVector.size()];
        for (int i=0; i<dataVector.size(); i++) {
            flag[i] = regexpfind(p,i,column);
        }
        return flag;
    }
    
    public boolean[] regexpfind(String regex, int[] row, int column) {
        if (row.length == 0) { return null; }
        Pattern p = Pattern.compile(regex);
        boolean[] flag = new boolean[row.length];
        for (int i=0; i<row.length; i++) {
            flag[i] = regexpfind(p,row[i],column);
        }
        return flag;
    }
    
    public boolean[] regexpfind(String regex, boolean[] flag, int column) {        
        int[] row = convertBooleanToRows(flag);
        return regexpfind(regex,row,column);
    }
    
    public boolean regexpfind(String regex, int row, int column) {
        Pattern p = Pattern.compile(regex);
        return regexpfind(p,row,column);        
    }
    
    public boolean regexpfind(Pattern p, int row, int column) {        
        return p.matcher((String) (((Vector) dataVector.get(row)).get(column))).find();
    }
    
    
    public String[] regexptoken(String regex, int column) {
        Pattern p = Pattern.compile(regex);
        String[] retstr = new String[dataVector.size()];
        for (int i=0; i<dataVector.size(); i++) {
            retstr[i] = regexptoken(p,i,column);
        }
//         if (retstr.length==1) { return retstr[0]; }
        return retstr;        
    }
    
    public String[] regexptoken(String regex, int[] row, int column) {
        if (row.length == 0) { return null; }
        Pattern p = Pattern.compile(regex);
        String[] retstr = new String[row.length];
        for (int i=0; i<row.length; i++) {
            retstr[i] = regexptoken(p,row[i],column);
        }
//         if (retstr.length==1) { return retstr[0]; }
        return retstr;
    }
    
    public String[] regexptoken(String regex, boolean[] flag, int column) {
        int[] row = convertBooleanToRows(flag);
        return regexptoken(regex,row,column);
    }
    
    public String regexptoken(String regex, int row, int column) {
        Pattern p = Pattern.compile(regex);
        return regexptoken(p,row,column);
    }
    
    public String regexptoken(Pattern p, int row, int column) {
        Matcher m = p.matcher((String) (((Vector) dataVector.get(row)).get(column)));
        if (m.find()) { return m.group(1); }
        return null;
    }
    
    

    public Object getColumnVector(int[] column) {
        Object obj = getColumnVector(column[0]);
        Object retarr = Array.newInstance(obj.getClass(),column.length);
        Array.set(retarr,0,obj);
        for (int i=1; i<column.length; i++) {
            Array.set(retarr,i,getColumnVector(column[i]));
        }
        return retarr;
    }
    
    public Object getColumnVector(int column) {
        Class C = getColumnClass(column);
        Object a = Array.newInstance(C,dataVector.size());
        for (int i=0; i<dataVector.size(); i++) {
            Array.set(a,i,((Vector) dataVector.get(i)).get(column));
        }
        return convertToPrimitives(a);
    }
    
    public Object getValuesAt(int[] row, int[] column) {
        if (row.length == 0) { return null; }
        Object obj = getValuesAt(row,column[0]);
        Object retarr = Array.newInstance(obj.getClass(),column.length);
        Array.set(retarr,0,obj);
        for (int i=1; i<column.length; i++) {
            Array.set(retarr,i,getValuesAt(row,column[i]));
        }
        return retarr;
    }
    
    public Object getValuesAt(int[] row, int column) {
        if (row.length == 0) { return null; }
        Class C = getColumnClass(column);
        Object a = Array.newInstance(C,row.length);
        for (int i=0; i<row.length; i++) {
            Array.set(a,i,((Vector) dataVector.get(row[i])).get(column));
        }
        return convertToPrimitives(a);
    }
    
//     public Object getValuesAt(int row, int column) {
//         return getValueAt(row,column);
//     }
    
    public Object getValuesAt(boolean[] flag, int[] column) {
        int[] row = convertBooleanToRows(flag);
        return getValuesAt(row,column);
    }
    
    public Object getValuesAt(boolean[] flag, int column) {
        int[] row = convertBooleanToRows(flag);
        return getValuesAt(row,column);
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
    
    public Object convertToPrimitives(Object a) {
        Object obj = Array.get(a,0);
//         if (Array.getLength(a)==1) { return obj; }
        if      (obj instanceof Boolean) { return ArrayUtils.toPrimitive((Boolean[])a); }
        else if (obj instanceof Byte)    { return ArrayUtils.toPrimitive((Byte[])   a); }
        else if (obj instanceof Short)   { return ArrayUtils.toPrimitive((Short[])  a); }
        else if (obj instanceof Integer) { return ArrayUtils.toPrimitive((Integer[])a); }
        else if (obj instanceof Long)    { return ArrayUtils.toPrimitive((Long[])   a); }
        else if (obj instanceof Float)   { return ArrayUtils.toPrimitive((Float[])  a); }
        else if (obj instanceof Double)  { return ArrayUtils.toPrimitive((Double[]) a); }
        else                             { return a; }
    }
    
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
    
    
    public void saveDataVector(String path) {
        try {
            ObjectOutputStream oos =
            new ObjectOutputStream(
            new BufferedOutputStream(
            new FileOutputStream(path)));
            oos.writeObject(getDataVector());
            oos.close();
        } catch (Throwable t) {
            t.printStackTrace();
        }
    }
    
    public void saveDataVectorAsync(final String path) {
        Thread thread = new Thread() { // asynchronized save data
            public void run() {
                saveDataVector(path);
            }
        };
        thread.start();
    }
    
    public void save(String path) {
        try {
            ObjectOutputStream oos =
            new ObjectOutputStream(
            new BufferedOutputStream(
            new FileOutputStream(path)));
            oos.writeObject(this);
            oos.close();
        } catch (Throwable t) {
            t.printStackTrace();
        }
    }
    
    public void saveAsync(final String path) {
        Thread thread = new Thread() { // asynchronized save data
            public void run() {
                save(path);
            }
        };
        thread.start();
    }
    
    public static MyTableModel load(String path) {
        try {
            ObjectInputStream ois =
            new ObjectInputStream(
            new BufferedInputStream(
            new FileInputStream(path)));
            MyTableModel mtm = (MyTableModel)ois.readObject();
            ois.close();
            return mtm;
        } catch (Throwable t) {
            t.printStackTrace();
            System.out.println("Error occured in MyTableModel.load()");
            return null;
        }
    }
    
}
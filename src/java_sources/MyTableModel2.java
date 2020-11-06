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
import org.apache.commons.lang3.math.NumberUtils;
import org.apache.commons.lang3.ArrayUtils;

public class MyTableModel2 extends DefaultTableModel{
    
    private static final long serialVersionUID = -7812759925336819450L;
    public Object oldobj;
    public Object newobj;
    public TreeMap<String,String[]> map = new TreeMap<String,String[]>();
    public Stack<History> undoStack = new Stack<History>();
    public Stack<History> redoStack = new Stack<History>();
    
    /* constructor */
    public MyTableModel2(){
        super();
    }
    
    public MyTableModel2(Object[][] data, Object[] columnNames){
        super(data, columnNames);
    }
    
    public MyTableModel2(Vector data, Vector columnNames){
        super(data, columnNames);
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
        super.setValueAt(obj,row,col);
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
        for (int r=0; r<row.length; r++) { for (int c=0; c<col.length; c++) {
            ((Vector) dataVector.get(row[r])).set(col[c],obj[r][c]); } }
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
        for (int r=0; r<row.length; r++) { for (int c=0; c<col.length; c++) {
            ((Vector) dataVector.get(row[r])).set(col[c],obj[r][c]); } }
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
        for (int r=0; r<row.length; r++) { for (int c=0; c<col.length; c++) {
            ((Vector) dataVector.get(row[r])).set(col[c],obj[r][c]); } }        
        if (flags.length>0 && flags[0]==true) { fireTableRowsUpdatedFromRows(row); }
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
            for (int r=0; r<row.length; r++) { for (int c=0; c<col.length; c++) {
                a[r][c] = (Boolean)(((Vector)(dataVector.get(row[r]))).get(col[c])); } }
            return a;
        }
        else if (obj instanceof Double)  {
            double[][] a = new double[row.length][col.length];
            for (int r=0; r<row.length; r++) { for (int c=0; c<col.length; c++) {
                a[r][c] = (Double)(((Vector)(dataVector.get(row[r]))).get(col[c])); } }
            return a;
        }
        else if (obj instanceof String)  {
            String[][] a = new String[row.length][col.length];
            for (int r=0; r<row.length; r++) { for (int c=0; c<col.length; c++) {
                a[r][c] = (String)(((Vector)(dataVector.get(row[r]))).get(col[c])); } }
            return a;
        }
        return null;
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
            for (int c=0; c<col.length; c++) {
                Vector v = (Vector)(dataVector.get(row[r]));
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
            for (int c=0; c<col.length; c++) {
                Vector v = (Vector)(dataVector.get(row[r]));
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
            for (int c=0; c<col.length; c++) {
                Vector v = (Vector)(dataVector.get(row[r]));
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
            for (int c=0; c<col.length; c++) {
                Vector v = (Vector)(dataVector.get(row[r]));
                Matcher m = p.matcher( (String)(v.get(col[c])) );
                if (m.find()) { a[r][c] = m.group(1); }
            }
        }        
        return a;
    }

    
    /* save and load */
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
    
    public static MyTableModel2 load(String path) {
        try {
            ObjectInputStream ois =
            new ObjectInputStream(
            new BufferedInputStream(
            new FileInputStream(path)));
            MyTableModel2 mtm = (MyTableModel2)ois.readObject();
            ois.close();
            return mtm;
        } catch (Throwable t) {
            t.printStackTrace();
            System.out.println("Error occured in MyTableModel.load()");
            return null;
        }
    }
    
    
    /* undo and redo */
    public void addUndo(String strdo, int row) {
        int[] rows = {row};
        addUndo(strdo,rows);
    }
    
    public void addUndo(String strdo,int[] rows) {
        Vector<Vector> vv = new Vector<Vector>();
        for (int r=0; r<rows.length; r++) {
            //vv.add((Vector) dataVector.get(rows[r])); }
            vv.add((Vector)(((Vector) dataVector.get(rows[r])).clone())); }
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
        Vector<Vector> vv = new Vector<Vector>();
        for (int r=0; r<rows.length; r++) {
            vv.add((Vector)(((Vector) dataVector.get(rows[r])).clone())); }
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
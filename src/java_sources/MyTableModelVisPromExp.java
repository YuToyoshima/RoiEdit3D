import javax.swing.table.DefaultTableModel;
import java.util.Vector;
import java.util.regex.Pattern;

public class MyTableModelVisPromExp extends DefaultTableModel{
   
    public MyTableModelVisPromExp(){
        super();
    }
   
    public MyTableModelVisPromExp(Object[][] data, Object[] columnNames){
        super(data, columnNames);
    }
    
    public MyTableModelVisPromExp(Vector data, Vector columnNames){
        super(data, columnNames);
    }

    @Override
    public Class getColumnClass(int col){
       return getValueAt(0, col).getClass();     
    }

    @Override
    public boolean isCellEditable (int row, int col) {
//    if(col==0) {
            return true;
//        }
//        else {
            //return false;
        //}        
    }

    public void setValuesAt(Object[] obj, int[] row, int[] column) {
        for (int i=0; i<obj.length; i++) {
            ((Vector) dataVector.get(row[i])).set(column[i],obj[i]);
        }
        fireTableDataChanged();
    }

    public boolean[] strcmp(String str, int column) {        
        boolean[] flag = new boolean[dataVector.size()];
        if (str!=null) {
            for (int i=0; i<dataVector.size(); i++) {
                flag[i] = str.equals(((Vector) dataVector.get(i)).get(column));
            }
        }
        return flag;
    }

    public int[] strcmp (String[] strs, int column) {
        int[] idx = new int[dataVector.size()];
        for (int j=0; j<dataVector.size(); j++) {
            idx[j] = -1; // initialize
        }
        for (int i=0; i<strs.length; i++) {
            boolean[] flag = this.strcmp(strs[i],column);
            for (int j=0; j<flag.length; j++){
                if (flag[j]) {
                    idx[j] = i;
                    //continue;
                }
            }
        }
        return idx;
    }

    
    public int[] strcmp (Vector strs, int column) {
        int[] idx = new int[dataVector.size()];
        for (int j=0; j<dataVector.size(); j++) {
            idx[j] = -1; // initialize
        }
        for (int i=0; i<strs.size(); i++) {
            boolean[] flag = this.strcmp(String.valueOf(strs.get(i)),column);            
            for (int j=0; j<flag.length; j++){            	
                if (flag[j]) {
                    idx[j] = i;
                    //continue;
                }
            }
        }
        return idx;        
    }


    public boolean[] regexp(String str, int column) {        
        boolean[] flag = new boolean[dataVector.size()];
        if (str!=null) {
            for (int i=0; i<dataVector.size(); i++) {
                flag[i] = Pattern.matches(str,String.valueOf(((Vector)dataVector.get(i)).get(column)));
            }
        }
        return flag;
    }

    public int[] regexp (String[] strs, int column) {
        int[] idx = new int[dataVector.size()];
        for (int j=0; j<dataVector.size(); j++) {
            idx[j] = -1; // initialize
        }
        for (int i=0; i<strs.length; i++) {
            boolean[] flag = this.regexp(strs[i],column);
            for (int j=0; j<flag.length; j++){
                if (flag[j]) {
                    idx[j] = i;
                    //continue;
                }
            }
        }
        return idx;
    }

    
    public int[] regexp (Vector strs, int column) {
        int[] idx = new int[dataVector.size()];
        for (int j=0; j<dataVector.size(); j++) {
            idx[j] = -1; // initialize
        }
        for (int i=0; i<strs.size(); i++) {
            boolean[] flag = this.regexp(String.valueOf(strs.get(i)),column);            
            for (int j=0; j<flag.length; j++){            	
                if (flag[j]) {
                    idx[j] = i;
                    //continue;
                }
            }
        }
        return idx;        
    }


    public Vector getDataColumn(int column) {
    	Vector vec = new Vector(dataVector.size());
    	for (int i=0; i<dataVector.size(); i++) {
    	    vec.add(((Vector)(dataVector.get(i))).get(column));    	              
    	}
    	return vec;
    }

}
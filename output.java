package micro_array;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;

public class output {
	
	public static String[] getProbes(String[][] data){
		
		String[] probes = new String[data.length];

	    for (int i = 1; i < data.length; i++){
	    	probes[i-1] = data[i][0];
	    }
	    return(probes);
	}
	
	public static void writeTXT(String[] probes, double[] values,String file) throws FileNotFoundException{
	
	      PrintStream out = new PrintStream(new FileOutputStream(file));
	      for (int i = 0; i < values.length-1; i++){
	        out.println(probes[i]+"			"+values[i]);
	      }
	      out.close();	      
	    }

}

///////////////////////////////////////////
// Einlesen und Analyse von .CEL-Dateien //
//                 von                   //
//       Nadine, Felix und Philipp       //
//               Gruppe 2                //
///////////////////////////////////////////
//				  OUTPUT				 //
///////////////////////////////////////////
package micro_array;

////////////
//Import  //
////////////
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;

public class output {
	
	/////////////////////////////
	// Isolieren der Probe-IDs //
	/////////////////////////////	
	
	public static String[] getProbes(String[][] data){
		
		String[] probes = new String[data.length];

	    for (int i = 1; i < data.length; i++){
	    	probes[i-1] = data[i][0];
	    }
	    return(probes);
	}

	//////////////////////////
	// Output in .txt-Datei //
	//////////////////////////
	public static void writeTXT(String[] probes,String[] express, double[] values,String file) throws FileNotFoundException{
	
	      PrintStream out = new PrintStream(new FileOutputStream(file));
	        out.println("probes p-values expressed");
	      for (int i = 0; i < values.length-1; i++){
	        out.println(probes[i]+" "+values[i]+" "+express[i]);
	      }
	      out.close();	      
	    }

}

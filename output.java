///////////////////////////////////////////
// Einlesen und Analyse von .CEL-Dateien //
//                 von                   //
//       Nadine, Felix und Philipp       //
//               Gruppe 2                //
///////////////////////////////////////////
//				  OUTPUT				 //
///////////////////////////////////////////


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
	public static void writeTXT(String[] probes,String[] express, double[] values, double[] slr,String[] geneSymbols, double[][] mas5, String[] Celnames, String file) throws FileNotFoundException{
		PrintStream out = new PrintStream(new FileOutputStream(file));
	    out.print("Probes geneSymbols p-Values Expressed SLR");
		for (int n = 1; n < Celnames.length-1; n++){
			out.print(" "+Celnames[n].substring(0, Celnames[n].lastIndexOf('.')));
		}
		out.println(" "+Celnames[Celnames.length-1].substring(0, Celnames[Celnames.length-1].lastIndexOf('.')));
	    for (int i = 0; i < values.length; i++){
	    	out.print(probes[i]+" "+geneSymbols[i]+" "+values[i]+" "+express[i]+" "+slr[i]);
		    for (int j = 0; j < mas5[0].length-1; j++){
		    	out.print(" "+mas5[i][j]);
		    }
		    out.println(" "+mas5[i][mas5[0].length-1]);
	    }
	    out.close();	      
	}



	//////////////////////////////////////
	// Output correlation in .txt-Datei //
	//////////////////////////////////////
	public static void writeCorrelationToTXT( double[][] data,String[] names,String file) throws FileNotFoundException{
		
		PrintStream out = new PrintStream(new FileOutputStream(file));
		out.println("Gene1 Gene1 Correlation");
		for (int i = 0; i < data.length; i++){
			for (int j = 0; j < data.length; j++){
				out.println(names[i]+" "+names[j]+" "+data[i][j]);
			}
		}
		out.close();	      
    }

	//////////////////////////////////////////
	// Output correlation all in .txt-Datei //
	//////////////////////////////////////////
	public static void writeCorrelationAllToTXT( double[] data,String[] names,String file, int Gene) throws FileNotFoundException{
		
		PrintStream out = new PrintStream(new FileOutputStream(file));
		out.println("Gene1 Gene1 Correlation");
		for (int i = 0; i < data.length; i++){
			out.println(names[Gene]+" "+names[i]+" "+data[i]);
		}
		out.close();	      
    }

	
	////////////////////////////
	// Write double[] to .txt //
	////////////////////////////
	
	public static void write1DDoubleToTXT (double[] data, String file ) throws FileNotFoundException{
		
		PrintStream out = new PrintStream(new FileOutputStream(file));
		for (int i = 0; i < data.length-1; i++){
			
			out.println(data[i]);
			
		}
		out.close();
	}
	
	
}


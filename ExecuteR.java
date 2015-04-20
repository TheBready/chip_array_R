///////////////////////////////////////////
// Einlesen und Analyse von .CEL-Dateien //
//                 von                   //
//       Nadine, Felix und Philipp       //
//               Gruppe 2                //
///////////////////////////////////////////
//				EXECUTE R				 //
///////////////////////////////////////////
package micro_array;

////////////
// Import //
////////////
import java.io.IOException;


public class ExecuteR {
	
	////////////////////////////////
	// Ausführen eines R-Skriptes //
	////////////////////////////////
	public static void runIt()throws IOException, InterruptedException{
		Process p = Runtime.getRuntime().exec("cmd runas /profile /user:Administrator /savecred /c start /wait Rscript micro_array_R\\readCel.R");
		System.out.println("Waiting for R ...");
	    p.waitFor();
	    System.out.println("R done.");	
	}
	
}

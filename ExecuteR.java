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


public class ExecuteR extends Thread{
	
	//////////////////////////////////
	// Konstruktor eines R-Skriptes //
	//////////////////////////////////
	ExecuteR(String path){
	setName(path);
	}
	
	
	////////////////////////////////
	// Ausführen eines R-Skriptes //
	////////////////////////////////
	public void run(){
		try {
			String rscriptpath = getName();
			Process p;
			p = Runtime.getRuntime().exec("cmd runas /profile /user:Administrator /savecred /c start /wait Rscript "+rscriptpath);
			System.out.println("Waiting for R with "+rscriptpath+" ...");
			p.waitFor();
	    	System.out.println("R done with "+rscriptpath);	
		} 
		catch (IOException | InterruptedException e) {
			e.printStackTrace();

		}
	}
	
}

package micro_array;


import java.io.IOException;


public class ExecuteR {
	
	////////////////////////////////
	// Ausf�hren eines R-Skriptes //
	////////////////////////////////
	public static void runIt()throws IOException, InterruptedException{
		Process p = Runtime.getRuntime().exec("cmd /c start /wait Rscript micro_array_R\\readCel.R");
		System.out.println("Waiting for R ...");
	    p.waitFor();
	    System.out.println("R done.");	
	}
	
}
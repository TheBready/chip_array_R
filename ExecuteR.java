package micro_array;


import java.io.*;


public class ExecuteR {
	
	
	
	public static void runIt()throws IOException, InterruptedException{
		Process p = Runtime.getRuntime().exec("cmd /c start /wait Rscript micro_array_R\\readCel.R");
		System.out.println("Waiting for R ...");
	    p.waitFor();
	    System.out.println("R done.");	
	}
	
	public static void main(String[] args) throws IOException, InterruptedException{
		
		runIt();
	}
}

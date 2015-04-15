package micro_array;


import java.io.*;


public class ExecuteR {
	
	public static void main(String[] args) throws IOException, java.io.FileNotFoundException{
		
		//System.setProperty("user.dir","C:\\Users\\Felix\\OwnCloud\\Studium\\7. Fachsemester\\Software-Praktikum\\R\\micro_array_R");
		//System.getProperty("user.dir");
		File currentDirectory = new File(new File(".").getAbsolutePath());
		System.out.println(currentDirectory.getCanonicalPath());
		System.out.println(currentDirectory.getAbsolutePath());
		
		Runtime.getRuntime().exec("cmd /c Rscript micro_array_R\\test.R");
		System.out.print("fertig");
		

	
	}
}

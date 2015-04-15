package micro_array;

////////////
// Import //
////////////
import java.io.*;
import java.util.*;



public class input {

	////////////////////////////
	// Einlesen der raw Daten //
	////////////////////////////
	public static String[][] readRaw(String table)throws java.io.FileNotFoundException,IOException{
		File file = new File(table);
        System.out.println("Import raw data...");
        Scanner input = new Scanner(file);
        final int maxLines = 604259;
        String[][] resultArray = new String[maxLines][];
        int linesCounter = 0;
        while (input.hasNextLine() && linesCounter < maxLines) {
            resultArray[linesCounter] = input.nextLine().split(" "); //Hier ist es ein Leerzeichen
            linesCounter++;
        }
        input.close();
        return(resultArray);
	}
	
	
	////////////////////////////
	// Einlesen der pm Daten //
	////////////////////////////
	public static String[][] readPm(String table)throws java.io.FileNotFoundException,IOException{
		File file = new File(table);
        System.out.println("Import pm data...");
        Scanner input = new Scanner(file);
        final int maxLines = 459426;
        String[][] resultArray = new String[maxLines][];
        int linesCounter = 0;
        while (input.hasNextLine() && linesCounter < maxLines) {
            resultArray[linesCounter] = input.nextLine().split(" "); //Hier ist es ein Leerzeichen
            linesCounter++;
        }
        input.close();
        return(resultArray);
	}
	
	////////////////////////////
	// Einlesen der mm Daten //
	////////////////////////////
	public static String[][] readMm(String table)throws java.io.FileNotFoundException,IOException{
		File file = new File(table);
        System.out.println("Import mm data...");
        Scanner input = new Scanner(file);
        final int maxLines = 459426;
        String[][] resultArray = new String[maxLines][];
        int linesCounter = 0;
        while (input.hasNextLine() && linesCounter < maxLines) {
            resultArray[linesCounter] = input.nextLine().split(" "); //Hier ist es ein Leerzeichen
            linesCounter++;
        }
        input.close();
        return(resultArray);
	}
	
	
	////////////////////////////////
	// Einlesen der MAS 5.0 Daten //
	////////////////////////////////
	public static String[][] readMAS5(String table)throws java.io.FileNotFoundException,IOException{
		File file = new File(table);
        System.out.println("Import MAS 5.0 data...");
        Scanner input = new Scanner(file);
        final int maxLines = 54676;
        String[][] resultArray = new String[maxLines][];
        int linesCounter = 0;
        while (input.hasNextLine() && linesCounter < maxLines) {
            resultArray[linesCounter] = input.nextLine().split("	"); //Hier ist es ein Tab
            linesCounter++;
        }
        input.close();
        return(resultArray);
	}
	
	//////////////////
	// Main-Methode //
	//////////////////
	public static void main(String[] args) throws IOException, java.io.FileNotFoundException{

		// Getwd()
		File currentDirectory = new File(new File(".").getAbsolutePath());
		System.out.println(currentDirectory.getCanonicalPath());
		System.out.println(currentDirectory.getAbsolutePath());
		
		// Aufrufen der Funktionen
		String[][] raw = readRaw("test.txt");
		String[][] mas5 = readMAS5("mas5.txt");
		String[][] pm = readPm("PM.txt");
		String[][] mm = readMm("MM.txt");
		
		
		// Test-Ausgaben
        System.out.println(raw.length);
        System.out.println(raw[1][1]);
        System.out.println(mas5.length);
        System.out.println(mas5[0][0]);
		
	}
	
	
	
	

}


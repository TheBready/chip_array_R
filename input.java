package micro_array;

////////////
// Import //
////////////
import java.io.File;
import java.io.IOException;
import java.util.Scanner;



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
	
	

}


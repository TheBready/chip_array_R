///////////////////////////////////////////
// Einlesen und Analyse von .CEL-Dateien //
//                 von                   //
//       Nadine, Felix und Philipp       //
//               Gruppe 2                //
///////////////////////////////////////////
//				  INPUT					 //
///////////////////////////////////////////
package micro_array;

////////////
// Import //
////////////
import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Scanner;


public class input  {
	
	///////////////////////////////////
	// Zählen der Zeilen einer Datei //
	///////////////////////////////////
	
	public static int countLines(String filename) throws IOException {
	    InputStream is = new BufferedInputStream(new FileInputStream(filename));
	    try {
	        byte[] c = new byte[1024];
	        int count = 0;
	        int readChars = 0;
	        boolean empty = true;
	        while ((readChars = is.read(c)) != -1) {
	            empty = false;
	            for (int i = 0; i < readChars; ++i) {
	                if (c[i] == '\n') {
	                    ++count;
	                }
	            }
	        }
	        return (count == 0 && !empty) ? 1 : count;
	    } finally {
	        is.close();
	    }
	}
	
	////////////////////////////
	// Einlesen der von Daten //
	////////////////////////////
	public static String[][] readFile(String table)throws java.io.FileNotFoundException,IOException{
		File file = new File(table);
        System.out.println("Import data...");
        Scanner input = new Scanner(file);
        int maxLines = countLines(table); 
        System.out.println("Datei hat "+maxLines+" Zeilen.");
        String[][] resultArray = new String[maxLines][];
        int linesCounter = 0;
        while (input.hasNextLine() && linesCounter < maxLines) {
            resultArray[linesCounter] = input.nextLine().split(" "); //Hier ist es ein Leerzeichen
            linesCounter++;
        }
        input.close();
        return(resultArray);
	}
	

}


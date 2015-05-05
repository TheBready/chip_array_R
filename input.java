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
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;


public class input extends Thread  {
	
    public String[][] inputString2D;

    
	//////////////////////////////
	// Konstruktor eines Inputs //
	//////////////////////////////
	input(String file){
	setName(file);
	}
    
	/////////////////////////////
	// Ausführen des Einlesens //
	/////////////////////////////
    
    public void run(){
		try {
			inputString2D = readFile(getName());
		} 
		catch (IOException e) {
			e.printStackTrace();

		}
	}

	
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
        System.out.println("Import data...");
        FileReader fileReader = new FileReader(table);
        BufferedReader bufferedReader = new BufferedReader(fileReader);
        int maxLines = countLines(table); 
        System.out.println("Datei hat "+maxLines+" Zeilen.");
        String[][] resultArray = new String[maxLines][];
        String line = null;
        int linesCounter = 0;
        while ((line = bufferedReader.readLine()) != null) {
            resultArray[linesCounter] = line.split(" ");
            linesCounter++;
        }
        bufferedReader.close();
        System.out.println("Daten importiert "+table);
        return(resultArray);
	}
	
	////////////////////////////////////////////
	// Einlesen von files über BufferedReader //
	////////////////////////////////////////////
	
	// returns BufferedReader object
	public static BufferedReader BfReader(String dir) {
		//just a dummy
		BufferedReader reader = null;	
		// read file 
		try{
			BufferedReader r = new BufferedReader(new FileReader(dir));
			return r;
				
		} catch (FileNotFoundException e) {
			System.out.println("Wrong file or directory.");
			e.printStackTrace();

		// dummy return to avoid error (function must always return a value)
		return reader;
	}
	
	
	}
	
}



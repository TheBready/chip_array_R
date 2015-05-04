///////////////////////////////////////////
// Einlesen und Analyse von .CEL-Dateien //
//                 von                   //
//       Nadine, Felix und Philipp       //
//               Gruppe 2  				 //
///////////////////////////////////////////
//				  MAIN					 //
///////////////////////////////////////////
package micro_array;

////////////
// Import //
////////////
import java.io.File;
import java.io.IOException;
import java.util.Arrays;

public class main {
	
	         
	//////////////////
	// Main-Methode //
	//////////////////
	public static void main(String[] args) throws IOException, InterruptedException{		
		
		// Getwd()
		File currentDirectory = new File(new File(".").getAbsolutePath());
		System.out.println(currentDirectory.getCanonicalPath());
		System.out.println(currentDirectory.getAbsolutePath());
		
		
		// Erstellen der Objekte der R-Skripte
        ExecuteR images = new ExecuteR("micro_array_R\\images.R");
        ExecuteR MVA = new ExecuteR("micro_array_R\\MVA.R");
        ExecuteR main_analyse = new ExecuteR("micro_array_R\\main-analyse.R");
        
		// Starten der R-Skripte
        main_analyse.start();
        MVA.start();
        images.start();
        

        // Warte auf Ende der R-Skripte
        MVA.join();
        main_analyse.join();
        images.join();

        
		
		// Laden der Dateien
		System.out.println("Laden der Daten");
	 	try{
	 		//String[][] raw =input.readRaw("output/ND_Group2_133Plus_2/exprs/ND_Group2_133Plus_2_signals.txt");
	 		String[][] mas5 = input.readFile("output/ND_Group2_133Plus_2/MAS5/ND_Group2_133Plus_2_MAS5_500.txt");
	 		//String[][] pm = input.readPm("output/ND_Group2_133Plus_2/pm/ND_Group2_133Plus_2_signals_PM.txt");
	 		//String[][] mm = input.readMm("output/ND_Group2_133Plus_2/mm/ND_Group2_133Plus_2_signals_MM.txt");
	 		String[][] pma = input.readFile("output/ND_Group2_133Plus_2/PMA/PMA_Calls.txt");
	 		
		
	 		// Namen der Chips
	 		System.out.println("Lese Chip-Namen");
	 		String[] Celnames = Arrays.copyOfRange(mas5[0], 0, mas5[0].length);
		
	 		// Namen der Probes 
	 		System.out.println("Lese Probeset-Namen");
	 		String[] mas5Names = output.getProbes(mas5);
		
	 		// Konvertiert in Matrix ohne erste Zeile und ohne erste Reihe
	 		System.out.println("Erstelle Matrix der Daten");
	 		double[][] mas5Double = micro_math.makeDouble(mas5);
	 		
	 		// Überprüft ob alle Present in Gruppe
	 		System.out.println("Erstelle Matrix der Daten");
	 		boolean[][] present = micro_math.isPresent(pma);
	 		
	 		// Filtert alle raus die nicht exprimiert werden
	 		System.out.println("Daten werden gefiltert");
	 		double[][] mas5_filtered = micro_math.filterIt(present,mas5Double);
	 		String[] probes_filtered = micro_math.filterProbes(present, mas5Names);
	 		
	 		
	 		//Erstellt Liste ob hoch oder niedrig exprimiert
	 		System.out.println("Teste ob höher exprimiert in Gruppe 1 als Gruppe 2");
	 		String[] express = micro_math.highOrLow(mas5_filtered);

	 		// T-Test
	 		System.out.println("t-Test wird durchgefürt");
	 		double[] mas5test = micro_math.studT(mas5_filtered);	 		
	 		System.out.println("Schreibe p-values in p-values.txt");
			File tdic = new File("output/ND_Group2_133Plus_2/t-Test/");
	        tdic.mkdir();
	 		output.writeTXT(probes_filtered,express,mas5test,"output/ND_Group2_133Plus_2/t-Test/p-values.txt");
	 		
	 		//Merge and Sort
	 		System.out.println("Bubble-Sort für p-values");
	 		micro_math.sortIt(mas5test,probes_filtered,express);
			System.out.println("Schreibe sortierte p-values in p-values_sorted.txt");
			output.writeTXT(probes_filtered,express,mas5test,"output/ND_Group2_133Plus_2/t-Test/p-values_sorted.txt");	
			
			micro_math.SLR("C:/Users/Felix/OwnCloud/Studium/7. Fachsemester/Software-Praktikum/workspace/micro_array");
			
			for(int i = 0; i < Celnames.length; i++){
				System.out.println(Celnames[i]);
			}

	 	}
		catch(IOException ex) {
			System.err.println("Kein Output-Ordner... Bitte lassen Sie erst readCel.R laufen!");
			ex.printStackTrace();
		}	
	}
}

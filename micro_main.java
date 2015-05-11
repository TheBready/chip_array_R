///////////////////////////////////////////
// Einlesen und Analyse von .CEL-Dateien //
//                 von                   //
//       Nadine, Felix und Philipp       //
//               Gruppe 2  				 //
///////////////////////////////////////////
//				  MAIN					 //
///////////////////////////////////////////


////////////
// Import //
////////////
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Date;

public class micro_main {
	
	         
	//////////////////
	// Main-Methode //
	//////////////////
	public static void main(String[] args) throws IOException, InterruptedException{		
		
		// Start der Laufzeit-Messung
		long start = new Date().getTime();

		// Getwd()
		File currentDirectory = new File(new File(".").getAbsolutePath());
		System.out.println(currentDirectory.getCanonicalPath());
		System.out.println(currentDirectory.getAbsolutePath());
		
		// Erstellen der Objekte der R-Skripte
        ExecuteR images = new ExecuteR("micro_array_R\\images.R");
        ExecuteR normalisation = new ExecuteR("micro_array_R\\normalisation.R");
        ExecuteR raw_analysis = new ExecuteR("micro_array_R\\raw-analysis.R");
        
		// Starten der R-Skripte
        raw_analysis.start();
        normalisation.start();
        images.start();
        
        // Warte auf Ende der R-Skripte
        normalisation.join();
        raw_analysis.join();
        images.join();        
		
		// Laden der Dateien
		System.out.println("Laden der Daten");
	 	try{

			// Erstellen der Objekte zum Input
	        input raw = new input("output/ND_Group2_133Plus_2/exprs/ND_Group2_133Plus_2_signals.txt");
	        input mas5 = new input("output/ND_Group2_133Plus_2/MAS5/ND_Group2_133Plus_2_MAS5_500.txt");
	        input pm = new input("output/ND_Group2_133Plus_2/pm/ND_Group2_133Plus_2_signals_PM.txt");
	        input mm = new input("output/ND_Group2_133Plus_2/mm/ND_Group2_133Plus_2_signals_MM.txt");
	        input pma = new input("output/ND_Group2_133Plus_2/PMA/PMA_Calls.txt");

	        // Starten des Einlesens
	        raw.start();
	        mas5.start();
	        pm.start();
	        mm.start();
	        pma.start();
	        

	        // Warte auf Ende des Einlesens
	        raw.join();
	        mas5.join();
	        pm.join();
	        mm.join();
	        pma.join();
	        	
	        
	 		// Namen der Chips
	 		System.out.println("Lese Chip-Namen");
	 		String[] Celnames = Arrays.copyOfRange(mas5.inputString2D[0], 0, mas5.inputString2D[0].length);

		
	 		// Namen der Probes 
	 		System.out.println("Lese Probeset-Namen");
	 		String[] mas5Names = output.getProbes(mas5.inputString2D);
		
	 		// Konvertiert in Matrix ohne erste Zeile und ohne erste Reihe
	 		System.out.println("Erstelle Matrix der Daten");
	 		double[][] mas5Double = micro_math.makeDouble(mas5.inputString2D);
	 		
	 		// Überprüft ob alle Present in Gruppe
	 		System.out.println("Erstelle Matrix der Daten");
	 		double treshold = 0.8; 
	 		boolean[][] present = micro_math.isPresent(pma.inputString2D,treshold);
	 		
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
			
			//Erstellen der SLR Datei aus MAS5 (Filter auf 0.2 gesetzt)
			String filePath = currentDirectory.getCanonicalPath();
			// Teststring
			double[] SlrResults = new double[micro_math.countLines(filePath + "\\output\\ND_Group2_133Plus_2\\MAS5\\ND_Group2_133Plus_2_MAS5_500.txt")];
	        SlrResults = micro_math.SLR(filePath, 0.2);
			System.out.println("Schreibe berechnete SLR-Werte in SLR_Values.txt");
			output.write1DDoubleToTXT(SlrResults,"output/ND_Group2_133Plus_2/SLR/pure_SLR_VALUES.txt");
			//output.writeTXT(probes_filtered,express,mas5test,"output/ND_Group2_133Plus_2/SLR/SLR_Values.txt");	
		
			// test for coexpressed Genes 
			System.out.println("Teste auf coexpremierte Gene");
			double[][] mas5correlation = micro_math.spearCorrelation(mas5_filtered,probes_filtered);
			System.out.println("Schreibe berechnete Correlation in correlation.txt");
			File dic = new File("output/ND_Group2_133Plus_2/coexpressed/");
	        dic.mkdir();
			output.writeCorrelationToTXT(mas5correlation,probes_filtered,"output/ND_Group2_133Plus_2/coexpressed/correlation.txt");
			
			for(int i = 0; i < Celnames.length; i++){
				System.out.println(Celnames[i]);
			}
			
			// Ende der Laufzeit-Messung
			long runningTime = new Date().getTime() - start; 
	        System.out.println("Laufzeit: "+(runningTime/1000)+" Sekunden");

		}
		catch(IOException ex) {
			System.err.println("Kein Output-Ordner... Bitte lassen Sie erst readCel.R laufen!");
			ex.printStackTrace();
		}	
	}
}

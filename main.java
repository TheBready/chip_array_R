package micro_array;

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
		
		
		//starten der R-Skriptes 
		ExecuteR.runIt();
				
	
		// Laden der Dateien
		System.out.println("Laden der Daten");
	 	try{
	 		String[][] raw =input.readRaw("output/ND_Group2_133Plus_2/exprs/ND_Group2_133Plus_2_signals.txt");
	 		String[][] mas5 = input.readMAS5("output/ND_Group2_133Plus_2/MAS5/ND_Group2_133Plus_2_MAS5_500.txt");
	 		String[][] pm = input.readPm("output/ND_Group2_133Plus_2/pm/ND_Group2_133Plus_2_signals_PM.txt");
	 		String[][] mm = input.readMm("output/ND_Group2_133Plus_2/mm/ND_Group2_133Plus_2_signals_MM.txt");
		
		
	 		// Namen der Chips
	 		System.out.println("Lese Chip-Namen");
	 		String[] Celnames = Arrays.copyOfRange(mas5[0], 0, mas5.length);
		
	 		// Namen der Probes 
	 		System.out.println("Lese Probeset-Namen");
	 		String[] mas5Names = output.getProbes(mas5);
		
	 		// Konvertiert in Matrix ohne erste Zeile und ohne erste Reihe
	 		System.out.println("Erstelle Matrix der Daten");
	 		double[][] mas5Double = micro_math.makeDouble(mas5);

	 		// T-Test
	 		System.out.println("t-Test wird durchgefürt");
	 		double[] mas5test = micro_math.studT(mas5Double);
	 		System.out.println("Schreibe p-values in p-values.txt");
	 		output.writeTXT(mas5Names,mas5test,"output/t-Test/p-values.txt");
	 		System.out.println("Bubble-Sort für p-values");
	 		micro_math.sortIt(mas5test,mas5Names);
			System.out.println("Schreibe sortierte p-values in p-values_sorted.txt");
			output.writeTXT(mas5Names,mas5test,"output/t-Test/p-values_sorted.txt");		
	 	}
		catch(IOException ex) {
			System.err.println("Kein Output-Ordner... Bitte lassen Sie erst readCel.R laufen!");
			ex.printStackTrace();
		}	
	}
}

package micro_array;

import java.io.File;
import java.io.IOException;
import java.util.*;
import micro_array.ExecuteR;
import micro_array.input;
import micro_array.micro_math;

public class main {
	//////////////////
	// Main-Methode //
	//////////////////
	public static void main (String[] args) throws IOException, InterruptedException {		
		
		// Getwd()
		File currentDirectory = new File(new File(".").getAbsolutePath());
		System.out.println(currentDirectory.getCanonicalPath());
		System.out.println(currentDirectory.getAbsolutePath());
		
		
		//starten der R-Skriptes 
		//ExecuteR.runIt();
					
		// Laden der Dateien
		//String[][] raw =input.readRaw("output/ND_Group2_133Plus_2/exprs/ND_Group2_133Plus_2_signals.txt");
		String[][] mas5 = input.readMAS5("output/ND_Group2_133Plus_2/MAS5/ND_Group2_133Plus_2_MAS5_500.txt");
		//String[][] pm = input.readPm("output/ND_Group2_133Plus_2/pm/ND_Group2_133Plus_2_signals_PM.txt");
		//String[][] mm = input.readMm("output/ND_Group2_133Plus_2/mm/ND_Group2_133Plus_2_signals_MM.txt");
			
		// Namen der Chips
		String[] Celnames = Arrays.copyOfRange(mas5[0], 0, mas5.length);
		
		// Konvertiert in Matrix ohne erste Zeile und ohne erste Reihe
		double[][] mas5Double = micro_math.makeDouble(mas5);

		
		double[] mas5test = micro_math.studT(mas5Double);
		
		
		// Test-Ausgaben
		//System.out.println(raw.length);
		//System.out.println(raw[1][1]);
		//System.out.println(mas5.length);
		//System.out.println(mas5[0][0]);
	}

}

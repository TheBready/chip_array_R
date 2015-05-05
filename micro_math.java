///////////////////////////////////////////
// Einlesen und Analyse von .CEL-Dateien //
//                 von                   //
//       Nadine, Felix und Philipp       //
//               Gruppe 2                //
///////////////////////////////////////////
//				MICRO_MATH				 //
///////////////////////////////////////////
package micro_array;

////////////
// Import //
////////////
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.LineNumberReader;
import java.util.Arrays;

import org.apache.commons.math3.stat.inference.TestUtils;


public class micro_math {
	
	////////////////////////////////////////////
	// convert string-array into double-array //
	////////////////////////////////////////////
	public static double[][] makeDouble(String[][] stringArray){
		
		double[][] doubleArray = new double[stringArray.length][6];
		
		for(int j=1;  j < stringArray[0].length; j++){
			for(int i=1;  i < stringArray.length; i++){	
				doubleArray[i-1][j-1] = Double.parseDouble(stringArray[i][j]);	
			}
		}
		
		return(doubleArray);
	}
	
	////////////
	// t-test //
	////////////
	
	public static double[] studT(double[][] data){
		
		double[][] sample1 = new double[data.length][3];
		double[][] sample2 = new double[data.length][3];
		double[] result = new double[data.length];

		for(int i=0;  i < data.length-1; i++){	
			sample1[i] = Arrays.copyOfRange(data[i], 0, 3);
			sample2[i] =  Arrays.copyOfRange(data[i], 3, 6);
			result[i] = TestUtils.pairedTTest(sample1[i], sample2[i]);
		}

		return(result);
	}
	
	/////////////////
	// Sort t-test //
	/////////////////
	public static void sortIt(double[] mas5test, String[] mas5Names, String[] express) {
		for (int n = 0; n < mas5test.length-1; n++) {
	        for (int m = 0; m < mas5test.length-2 - n; m++) {
	            if ((mas5test[m]-mas5test[m + 1]) > 0) {
	                double swapDouble = mas5test[m];
	                mas5test[m] = mas5test[m + 1];
	                mas5test[m + 1] = swapDouble;
	                String swapString = mas5Names[m];
	                mas5Names[m] = mas5Names[m + 1];
	                mas5Names[m + 1] = swapString;
	                String swapExpr = express[m];
	                express[m] = express[m + 1];
	                express[m + 1] = swapExpr;
	            }
	        }
	    }	  
	}
	
	///////////////////////////
	// High or low expressed //
	///////////////////////////
	public static String[] highOrLow(double[][] mas5Double) {
		String[] express = new String[mas5Double.length];
		
		for (int n = 0; n < mas5Double.length-1; n++) {
			double group1 = (mas5Double[n][0]+mas5Double[n][1]+mas5Double[n][2])/3;
			double group2 = (mas5Double[n][3]+mas5Double[n][4]+mas5Double[n][5])/3;
			if (group1 > group2){
				express[n] = "high";
			}
			else{
				express[n] = "low";

			}
	    }
		return(express);
	}
	
	////////////////////////
	// is a probe present //
	////////////////////////
	public static boolean[][] isPresent(String[][] pma, double treshold) {
		boolean[][] present = new boolean[pma.length][2];
		for (int n = 1; n < pma.length-1; n++) {
			double group1count = 0;
			double group2count = 0;
			boolean group1 = false;
			boolean group2 = false;

			// Gruppe 1
			if(pma[n][1].charAt(0)=='P'){
				group1count++;
			}
			if(pma[n][2].charAt(0)=='P'){
				group1count++;			
			}
			if(pma[n][3].charAt(0)=='P'){
				group1count++;
			}	
			
			// Gruppe 2
			if(pma[n][4].charAt(0)=='P'){
				group2count++;
			}
			if(pma[n][5].charAt(0)=='P'){
				group2count++;
			}
			if(pma[n][6].charAt(0)=='P'){
				group2count++;		
			}
			
			//Auswertung
			if((group1count/3)>=treshold){
				group1 = true;
			}
			if((group2count/3)>=treshold){
				group2 = true;
			}
			present[n-1][0] = group1;
			present[n-1][1] = group2;
	    }
		return(present);
	}
	
	
	///////////////////////
	// Filter input Data //
	///////////////////////
	public static double[][] filterIt(boolean[][] present, double[][] mas5Double) {
		int size = 0; 
		for (int i = 0; i < present.length; i++) {
			if(present[i][0]||present[i][1]){
				size++;
			}
		}
		int counter = 0;
		double[][] mas5_filtered = new double[size][2];
		for (int n = 0; n < present.length; n++) {
			if(present[n][0]||present[n][1]){
				mas5_filtered[counter] = Arrays.copyOfRange(mas5Double[n], 0, 6);
				counter++;
			}
	    }
		return(mas5_filtered);
	}
	
	//////////////////////
	// Filter Probe-IDs //
	//////////////////////
	public static String[] filterProbes(boolean[][] present, String[] probes) {
		int size = 0; 
		for (int i = 0; i < probes.length; i++) {
			if(present[i][0]||present[i][1]){
				size++;			
			}
		}
		int counter = 0;
		String[] probes_filtered = new String[size];
		for (int n = 0; n < present.length; n++) {
			if(present[n][0]||present[n][1]){
				probes_filtered[counter] = probes[n];
				counter++;
			}
	    }
		return(probes_filtered);
	}
    ///////////////////////////////
	// Log to base 2 calculation //
	///////////////////////////////
	
	// needed for SLR
	public static double logBase2(double var){
		
		return Math.log(var) / Math.log(2);
	}
	
	
	//////////////////////////////
	// countLines of a Document //		(online Source : http://stackoverflow.com/questions/453018/number-of-lines-in-a-file-in-java)
	//////////////////////////////
	
	public static int countLines(String dir) throws IOException {
		
		LineNumberReader  lnr = new LineNumberReader(new FileReader(new File(dir)));
		lnr.skip(Long.MAX_VALUE);
		lnr.close();
		return lnr.getLineNumber()+1;
	}
	
	////////////////////////////////
	// convert String[] to double //
	////////////////////////////////
	
	// gets a line of MAS5 and converts it's values to double values + returns
	public static double convMAS5toDouble(String[] line){
		
		// save converted values in double array
		double[] doubleParser = new double[line.length];
		
		// convert
		for (int i = 1; i < line.length; i++) {
			doubleParser[i] = Double.parseDouble(line[i]);
		}
		
		// return value calculated by slrMAS5Value function
		return slrMAS5Value(doubleParser);
	}
	
	/////////////////////////////
	// calc SLR (MAS5 - 2sets) //
	/////////////////////////////
	public static double slrMAS5Value(double[] convLine){
		
		// where to devide the sets
		int setBoarder = (convLine.length/2) + (convLine.length % 2); // odd length -> name : value-value | value-value, even length -> name: value | value-value
		
		double set1 = 0.0;
		double set2 = 0.0;
		
		for (int i = 1; i < setBoarder ; i++) {
			set1 = set1+convLine[i];
			set2 = set2+convLine[setBoarder+i-1];
		}
		
		double medianSet1 = set1/(setBoarder-1);
		double medianSet2 = set2/(setBoarder-1);
		
		// final log2 based calculation
		double slr = logBase2(medianSet1)/logBase2(medianSet2);
		
		return slr;
	}
	
	///////////////////////
	// SLR - Calculation //
	// (Signal Log Ratio)//
	///////////////////////
	

	public static void SLR(String folder, double thresholdFilter){
			
		// get directories ( works just for our path-system)
		String MAS5dir = (folder + "\\output\\ND_Group2_133Plus_2\\MAS5\\ND_Group2_133Plus_2_MAS5_500.txt");
		String PMAdir = (folder + "\\output\\ND_Group2_133Plus_2\\PMA\\PMA_Calls.txt"); // PMA file
	
		//create your BufferedReaders to get specific content (see code below ;))
		BufferedReader rPMA;
		BufferedReader rMAS5;		

		// create output Folder (Windows)
		// if used for other OS, change output
		File f = null;
		f = new File(folder + "\\output\\ND_Group2_133Plus_2\\SLR");
		f.mkdirs(); // creates the folder
				
		////////////////////////////
		// get chip names & count //
		////////////////////////////
		
		int size 	= 0;
		String chipNames[] = null;
		// save chip names (from PMA_calls.txt)

		// read PMA file (created by R_script) to get the used amount of chips
		try{
			// get new FileReader object to get every different line
			rPMA = input.BfReader(PMAdir);		
			// get columns count
			chipNames = rPMA.readLine().split(" ");
			size = chipNames.length;
			// close stream
			rPMA.close();
			
		// catch exceptions			
		} catch (FileNotFoundException e) {
			System.out.println("Wrong file or directory.");
			e.printStackTrace();
		} catch (IOException e) {
			System.out.println("Couldnt read file.");
			e.printStackTrace();			
		}

		/////////////////////////////////
		// count absent genes in chips //
		/////////////////////////////////
		
		int countedLines = 0;
		// get Lines from document (needed for correct count calculations -> counter needs to be an Array of countedLines size)
		try {
			countedLines = countLines(PMAdir);
		} catch (IOException e) {
			e.printStackTrace();
		} 
		// counts for every line, how many values where absent
		int[] counter = new int[countedLines];
		
		
		try{
			// read PMA.txt file content into buffer
			rPMA = input.BfReader(PMAdir);
			String countStr;
			// to get the current position in the file
			int currentLine = 0;
			// save every line separate in countStr until there is no more 
			while((countStr = rPMA.readLine())!=null) {
				// count all absents in PMA.txt
				for (int i=0; i < countStr.length()-1; ++i){
					// look at PMA file to understand (typical PMA file line : "1552275_s_at A M M P P P") -> every letter represents one chiparray
					if(countStr.charAt(i)== ' ' && countStr.charAt(i+1)=='A'){
						counter[currentLine]++;
					} // close if
				// close for (int i=0; i < countStr.length()-1; ++i){
				}
				currentLine++;
			}// close while((countStr = rPMA.readLine())!=null) {
			// close stream
			rPMA.close();
		
		// catch exceptions
		} catch (FileNotFoundException e) {
			System.out.println("Wrong file or directory.");
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		
		//////////////////////////////////
		// create Filtered_MAS5_500.txt //
		//////////////////////////////////
				
		
		int currentLine = 0;
		// get every line, which is higher than the threshold and store it in Filtered_MAS5_500.txt
		try {
			// create and write to file (windows)
			rMAS5 = input.BfReader(MAS5dir);
			//FileWriter PMA = null;
			FileWriter MAS5 = null;
			// create file Filtered_PMA.txt to get just the fitting values
			MAS5 = new FileWriter(""+f+"\\Filtered_MAS5_500.txt");
			//to store every line for one loop
			String readLine;
			// get every line from MAS5_500.txt and filter it based on the wanted value (parameter thresholdFilter)
			while(( readLine = rMAS5.readLine())!= null){
					//if more than 80% not absent 
		            if ( (1- ((double)counter[currentLine] / (double)(size))) > thresholdFilter){
		            	
		            	// write lines to file
						MAS5.write(readLine);
						// to  save every separate line
						MAS5.append( System.lineSeparator() );
													
		            }	// end if ( (1- ((double)counter[currentLine] / (double)(size))) > thresholdFilter){
				currentLine++;
				// writes buffer data to file but doesn't close it (for next while loop)
				MAS5.flush();				
			}	// end while			
			// close stream
			MAS5.close();
			rMAS5.close();
		// catch exceptions	
		} catch (FileNotFoundException e) {
			System.out.println("Wrong file or directory.");
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
				
		//////////////////////////////////
		// finally calculate SLR Values //
		//////////////////////////////////		
		// slr value
		double slrValue;		
		// reset currentLine
		currentLine = 0;		
		// new file writer
		FileWriter SLR = null;
			

		try {
			
			rMAS5 = input.BfReader(""+f+"\\Filtered_MAS5_500.txt");
			SLR = new FileWriter(""+f+"\\SLR_Values.txt");
			//save all slr-values in SLR_Values.txt
			
			// store all 
			String[] splitLine = null;
			String readLine;

			// get every line from Filtered_PMA.txt
			while((readLine = rMAS5.readLine())!=null){
				
				// ignore first line
				if (currentLine != 0) {

				// MAS5 variant:	
				// get 2 sets of chips (example : 6 chips total ->  set1 = 1-3, set2 = 4-6)
				// use the mean value of both sets for each gene to logBase2 (siehe functions convMAS5toDouble and slrMAS5Value)
					
					splitLine = readLine.split(" ");
					//System.out.println(splitLine[1]);					
					slrValue = convMAS5toDouble(splitLine);
					// write gene name + slr value to file
					SLR.write(splitLine[0] + "\t" + slrValue);
					SLR.append( System.lineSeparator() );
				
				} // end if (currentLine != 0) {
				currentLine++;
			} // end while loop (line 311)
			// close writing stream
			SLR.close();
			
		// catch exceptions	
		}catch (FileNotFoundException e) {
			System.out.println("Wrong file or directory.");
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}		
		
	} // end SLR
	
	
	
}


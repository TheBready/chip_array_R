///////////////////////////////////////////
// Einlesen und Analyse von .CEL-Dateien //
//                 von                   //
//       Nadine, Felix und Philipp       //
//               Gruppe 2                //
///////////////////////////////////////////
//				MICRO_MATH				 //
///////////////////////////////////////////


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
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

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

	/////////////////////////////
	// Check coexpressed Genes //
	/////////////////////////////
	public static double[][] pearCorrelationFiltered(double[][] data,String[] names){
		int rows = data.length;
		double[][] correlation = new double[rows][rows];
		for(int i=0;  i < rows; i++){
			for(int j=0;  j < rows; j++){
				correlation[i][j] = pearCorrelation(data[i],data[j]);
			}
		}	
		return(correlation);
	}
	
	/////////////////////////////////
	// Check all coexpressed Genes //
	/////////////////////////////////
	public static void pearCorrelationAll(double[][] data,String[] names){	
		// Number of CPU-Cores
		int cores = Runtime.getRuntime().availableProcessors();
		// Creates the threadpool
		ExecutorService threadPool = Executors.newFixedThreadPool(cores);
		//Anzahl der Gene
		int rows = data.length;	
		
		for (int i = 0; i < rows; i++) {
			int threadNumber = i;
			threadPool.submit(new Runnable() {
				public void run() {
					System.out.println("Thread "+threadNumber+" gestartet für "+names[threadNumber]+" ...");
					double[] correlation = new double[rows];
					for(int j=0;  j < rows; j++){
						correlation[j] = micro_math.pearCorrelation(data[threadNumber],data[j]);
					}
					try {
						output.writeCorrelationAllToTXT(correlation,names,"output/ND_Group2_133Plus_2/coexpressed/"+names[threadNumber]+".txt",threadNumber);
					} catch (FileNotFoundException e) {
						e.printStackTrace();
					}
				}
			});
		}
		
		
		// Schließe den threadpool
		threadPool.shutdown();
		
		// Warte auf alle threads
        while (!threadPool.isTerminated()) {
        }
        
        System.out.println("Finished all threads");
	}
	
	/*
	 for(int j=0;  j < rows; j++){
			correlation[j] = micro_math.pearCorrelation(data[i],data[j]);
	 }
	try {
		System.out.println("Schreibe berechnete Correlation von Gen "+names[i]+" in coexpressed/"+names[i]+".txt");
		output.writeCorrelationAllToTXT(correlation,names,"output/ND_Group2_133Plus_2/coexpressed/"+names[i]+".txt",i);
	} catch (FileNotFoundException e) {
		e.printStackTrace();
	}
	
	*/
	/////////////////////////
	// Pearson Correlation //
	/////////////////////////
	public static double pearCorrelation(double[] sample1, double[] sample2){
		//cor(X, Y) = (Sum[(xi - E(X))(yi - E(Y))] / [s(X)s(Y)])/n
		double correlation = 0;
		double sd1 = sd(sample1);
		double sd2 = sd(sample2);
		double mean1 = mean(sample1);
		double mean2 = mean(sample2);
		int elements = sample1.length;
		for(int i=0;  i < elements; i++){
			correlation = correlation + ((sample1[i]-mean1)*(sample2[i]-mean2));
		}
		
		correlation = ((correlation) / (sd1*sd2))/elements;
		return(correlation);
	}
	
	//////////
	// Mean //
	//////////
	public static double mean(double[] list){
		double mean = 0;
		int elements = list.length;
		for(int i=0;  i < elements; i++){
			mean = mean + list[i];
		}
		
		mean = mean/elements;
		return(mean);
	}
	
	//////////////
	// Variance //
	//////////////
	public static double var(double[] list){
		double var = 0;
		int elements = list.length;
		double meanList = mean(list);
		for(int i=0;  i < elements; i++){
			double step = list[i]-meanList;
			var = var + step*step;  
		}
		
		var = var/elements;
		return(var);
	}

	/////////////////////////
	// Standard-derivation //
	/////////////////////////
	public static double sd(double[] list){
		double var = var(list);
		double sd = Math.sqrt(var);
		return(sd);
	}
	
	////////////
	// t-test //
	////////////
	
	public static double[] studT(double[][] data){
		
		double[][] sample1 = new double[data.length][3];
		double[][] sample2 = new double[data.length][3];
		double[] result = new double[data.length];

		for(int i=0;  i < data.length; i++){	
			sample1[i] = Arrays.copyOfRange(data[i], 0, 3);
			sample2[i] =  Arrays.copyOfRange(data[i], 3, 6);
			result[i] = TestUtils.pairedTTest(sample1[i], sample2[i]);
		}

		return(result);
	}
	
	/////////////////
	// Sort t-test //
	/////////////////
	public static void sortIt(double[] mas5test, String[] mas5Names, String[] express, double[] slr, String[] symbols, double[][] mas5_filtered) {
		for (int n = 0; n < mas5test.length; n++) {
	        for (int m = 0; m < mas5test.length-1 - n; m++) {
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
	                double swapSLR = slr[m];
	                slr[m] = slr[m + 1];
	                slr[m + 1] = swapSLR;
	                String swapSymbols = symbols[m];
	                symbols[m] = symbols[m + 1];
	                symbols[m + 1] = swapSymbols;
	                double[] swapmas5 = mas5_filtered[m];
	                mas5_filtered[m] = mas5_filtered[m + 1];
	                mas5_filtered[m + 1] = swapmas5;
	            }
	        }
	    }	  
	}
	
	///////////////////////////
	// High or low expressed //
	///////////////////////////
	public static String[] highOrLow(double[][] mas5Double) {
		String[] express = new String[mas5Double.length];
		
		for (int n = 0; n < mas5Double.length; n++) {
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
		for (int n = 1; n < pma.length; n++) {
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
	public static double[][] filterIt(boolean[][] present, double[][] mas5Double, double[] slr) {
		int size = 0; 
		for (int i = 0; i < present.length; i++) {
			if((present[i][0]||present[i][1])&&(slr[i]>2.0)){
				size++;
			}
		}
		int counter = 0;
		double[][] mas5_filtered = new double[size][2];
		for (int n = 0; n < present.length; n++) {
			if((present[n][0]||present[n][1])&&(slr[n]>2.0)){
				mas5_filtered[counter] = Arrays.copyOfRange(mas5Double[n], 0, 6);
				counter++;				
			}
	    }
		return(mas5_filtered);
	}
	
	////////////////
	// Filter SLR //
	////////////////
	public static double[] filterItSLR(boolean[][] present, double[] slr) {
		int size = 0; 
		for (int i = 0; i < present.length; i++) {
			if((present[i][0]||present[i][1])&&(slr[i]>2.0)){
				size++;
			}
		}
		int counter = 0;
		double[] slr_filtered = new double[size];
		for (int n = 0; n < present.length; n++) {
			if((present[n][0]||present[n][1])&&(slr[n]>2.0)){
				slr_filtered[counter] = slr[n];
				counter++;				
			}
	    }
		return(slr_filtered);
	}

	
	
	
	//////////////////////
	// Filter Probe-IDs //
	//////////////////////
	public static String[] filterProbes(boolean[][] present, String[] probes, double[] slr) {
		int size = 0; 
		for (int i = 0; i < probes.length; i++) {
			if((present[i][0]||present[i][1])&&(slr[i]>2.0)){
				size++;			
			}
		}
		int counter = 0;
		String[] probes_filtered = new String[size];
		for (int n = 0; n < present.length; n++) {
			if((present[n][0]||present[n][1])&&(slr[n]>2.0)){
				probes_filtered[counter] = probes[n];
				counter++;
			}
	    }
		return(probes_filtered);
	}
	
	//////////////////////
	// decode Probe-IDs //
	//////////////////////
	public static String[] decodeProbes(String[][] symbols, String[] probes) {
		String[] probes_decoded = new String[probes.length];
		for (int i = 0; i < probes.length; i++) {
			for (int j = 1; j < symbols.length; j++) {
				if(symbols[j][1].equals(probes[i])){
					probes_decoded[i] = symbols[j][2];
					break;
				}
				else if (j+1 == symbols.length){
					probes_decoded[i] = "unknown";
					break;
		        }
			}
		}
		return(probes_decoded);
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
	

	public static double[] SLR(String folder){
			
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
		// Get Line-Count in PMA-file  //
		/////////////////////////////////
		
		int countedLines = 0;
		// get Lines from document (needed for correct count calculations -> counter needs to be an Array of countedLines size)
		try {
			countedLines = countLines(PMAdir);
		} catch (IOException e) {
			e.printStackTrace();
		} 
		
		
		int currentLine = 0;
	
		//////////////////////////////////
		// finally calculate SLR Values //
		//////////////////////////////////		
		// slr value
		double slrValue;		
		// reset currentLine
		currentLine = 0;		
		// new file writer
		FileWriter SLR = null;
			
		double[] output = new double[countedLines];

		try {
			
			//rMAS5 = input.BfReader(""+f+"\\Filtered_MAS5_500.txt");
			rMAS5 = input.BfReader(MAS5dir);
			//SLR = new FileWriter(""+f+"\\SLR_Values.txt");
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
				
					output[currentLine-1] = slrValue;	
					
					
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
		
		return output;
	} // end SLR
	
}


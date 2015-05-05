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
import java.util.ArrayList;
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
	public static boolean[][] isPresent(String[][] pma) {
		boolean[][] present = new boolean[pma.length][2];
		for (int n = 1; n < pma.length-1; n++) {
			boolean group1 = (pma[n][1].charAt(0)=='P'&&pma[n][2].charAt(0)=='P'&&pma[n][3].charAt(0)=='P');
			boolean group2 = (pma[n][4].charAt(0)=='P'&&pma[n][5].charAt(0)=='P'&&pma[n][6].charAt(0)=='P');
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
	
	// gets the basic folder directory
	// a Signal Log Ratio of 1.0 indicates an increase of the transcript level by 2 fold and -1.0 indicates a decrease by 2 fold
	// creates folder SLR and files with all SLR values (calculated out of MAS5 and PMA)
	// create table with R (call SLR_script)
	public static void SLR(String folder){
		
		
		// get directories ( works just for our path-system)
		String MAS5dir = (folder + "\\output\\ND_Group2_133Plus_2\\MAS5\\ND_Group2_133Plus_2_MAS5_500.txt");
		String PMAdir = (folder + "\\output\\ND_Group2_133Plus_2\\PMA\\PMA_Calls.txt"); // PMA file
		
		//create your BufferedReaders to get specific content (see later code)
		BufferedReader r;
		BufferedReader r1;
		BufferedReader r2;
		BufferedReader r3 = null;
		

		// create output Folder (Windows)
		// if used for other OS, change output
		File f = null;
		f = new File(folder + "\\output\\ND_Group2_133Plus_2\\SLR");
		f.mkdirs(); // creates the folder
		
		
		int counter = 0;
		// get the amount of columns in our files
		int size 	= 0;
		
		// read MAS5 file (created by R_script)	
		try{
			// get new FileReader object to get every different line
			r1 = new BufferedReader(new FileReader(MAS5dir));
			
			// get columns count
			size = r1.readLine().split("\t").length;
			
			// close stream
			r1.close();
		// catch exceptions			
		} catch (FileNotFoundException e) {
			System.out.println("Wrong file or directory.");
			e.printStackTrace();
		} catch (IOException e) {
			System.out.println("Couldnt read file.");
			e.printStackTrace();
					
		}

		// read PMA file (created by R_script)
		
		try{
			// read PMA.txt file content into buffer
			r2 = new BufferedReader(new FileReader(PMAdir));
			String countStr;
			// save every line separate in countStr until there is no more 
			while((countStr = r2.readLine())!=null) {
				// count all absents in PMA.txt
				//System.out.println(counter); -> 229683 complete
				for (int i=0; i < countStr.length()-1; ++i){
					// look at PMA file to understand (typical PMA file line : "1552275_s_at A M M P P P") -> every letter represents one chiparray
					if(countStr.charAt(i)== ' ' && countStr.charAt(i+1)=='A'){
						counter++;
					} // close if
				
				}	// close for
			}// close while	
			// close stream
			r2.close();
		
			
		// catch exceptions
		} catch (FileNotFoundException e) {
			System.out.println("Wrong file or directory.");
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		

		try {
			
			// create and write to file (windows)
			r = new BufferedReader(new FileReader(MAS5dir));

			FileWriter fw = null;
			
			// create file SLR_MAS5_500.txt -> we use MAS5 with 500
			fw = new FileWriter( ""+f+"\\SLR_MAS5_500.txt" );
		
			String line;
			
			// get every line from SLR_MAS5_500.txt
			while(( line =r.readLine())!= null){
					//if more than 80% not absent
		            if((counter / (size-1))>0.2){
		            	//System.out.println(line); //testprint
						fw.write(line);
						// to  keep the right format
						fw.append( System.getProperty("line.separator") );
													
		            }	// end if
					
					// writes buffer data to file but doesn't close it (for next while loop)
					fw.flush();
			}	// end while
			
			// close stream
			fw.close();
			
		// catch exceptions	
		} catch (FileNotFoundException e) {
			System.out.println("Wrong file or directory.");
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		
		// slr value
		Double slrValue; 
		
		// new file writer
		FileWriter w = null;
		
		// used in while loop (line 311) to check in which line we are in
		int countLine = 0;
		
		try{
			// read SLR_MAS5_500 and load into buffer
			r3 = new BufferedReader(new FileReader(""+f+"\\SLR_MAS5_500.txt"));
			
			String[] line;
			String readLine;
			
			w = new FileWriter(""+f+"\\SLR_MAS5_500_values.txt");
			//save all slr-values in SLR_MAS5_500_values.txt
			
			// 
			ArrayList <String> chipNames = new ArrayList<String>();
			

			// get every line from SLR_MAS5_500.txt
			while((readLine=r3.readLine())!=null){
				

				// saves in line every sub string of the line
				line = readLine.split("\t");
				// for the first line, where our chip-Names are located, we save them in our chipNames string

				if (countLine == 0){
					for(int l=0; l<line.length; ++l){
						chipNames.add(line[l]);
						// System.out.println(line[i]); // would print our chipNames
					}
				}
				// not first line anymore
				countLine++;

				// use logBase2(double var) to get the log2 values

				if (countLine > 1){
					for (int i = 1; i < line.length; ++i){
						for( int j = 1; j < line.length; ++j){

							if(i!=j){

								if(Double.parseDouble(line[i])<50&&Double.parseDouble(line[j])>500){
									double d1 =Double.parseDouble(line[i]);

									double d2 =Double.parseDouble(line[j]);

									//calculate slr value
									slrValue =(logBase2(d1))-(logBase2(d2));
									// convert to string to save in file
									String SLRstring = slrValue.toString();
									// write slr in file
									w.write(line[0]+"\t" + chipNames.get(i-1) + "\t" + chipNames.get(j-1) + ": " + SLRstring);
									w.append( System.getProperty("line.separator") );

								} // end inner if
								if(Double.parseDouble(line[i])>500&&Double.parseDouble(line[j])<50){
									//calculate slr value
									slrValue=logBase2(Double.parseDouble(line[i])-logBase2(Double.parseDouble(line[j])));
									// convert to string to save in file
									String slrstring = slrValue.toString();
									// write slr in file
									w.write(line[0]+"\t"+chipNames.get(i-1)+"\t"+chipNames.get(j-1)+": "+slrstring);									
									w.append( System.getProperty("line.separator") );

								} // end 2. inner if
							} // end if
						} // end inner for loop
					} // end outer for loop
	
				} // end outer if 
			
			
			} // end while loop (line 311)
			// close writing stream
			w.close();
			
		// catch exceptions	
		}catch (FileNotFoundException e) {
			System.out.println("Wrong file or directory.");
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		
		
		
		
		
	} // end SLR
}


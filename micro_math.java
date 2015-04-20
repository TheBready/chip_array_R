package micro_array;

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
			System.out.print(sample1[i][2]);
			System.out.println("   "+sample2[i][2]+"    "+result[i]);
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
}

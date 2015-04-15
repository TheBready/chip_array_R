package micro_array;

import java.util.Arrays;
import org.apache.commons.math3.stat.inference.TestUtils;

public class micro_math {
	
	public static double[][] makeDouble(String[][] stringArray){
		
		double[][] doubleArray = new double[stringArray.length][6];
		
		for(int j=1;  j < stringArray[0].length; j++){
			for(int i=1;  i < stringArray.length; i++){	
				doubleArray[i-1][j-1] = Double.parseDouble(stringArray[i][j]);	
			}
		}
		
		return(doubleArray);
	}
	
	
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
	
	
	public static void sortIt(double[] mas5test, String[] mas5Names) {
		for (int n = 0; n < mas5test.length-1; n++) {
	        for (int m = 0; m < mas5test.length-2 - n; m++) {
	            if ((mas5test[m]-mas5test[m + 1]) > 0) {
	                double swapDouble = mas5test[m];
	                mas5test[m] = mas5test[m + 1];
	                mas5test[m + 1] = swapDouble;
	                String swapString = mas5Names[m];
	                mas5Names[m] = mas5Names[m + 1];
	                mas5Names[m + 1] = swapString;
	            }
	        }
	    }	  
	}
	

}

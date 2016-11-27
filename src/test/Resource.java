//: test/Resource.java
package test;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.ArrayList;

/**
 * @author chn
 * 
 * Resource: class for load source data from file
 */
public class Resource {
	
	// Reference: George Karypis et al.. Chameleon: Hierarchical Clustering Using Dynamic Modeling.
	// Reference: http://archive.ics.uci.edu/ml/index.html
	private static final String[] DATA_PATHS = {
		// Chameleon
		"/resource/chameleon/t4.8k.txt", 
		"/resource/chameleon/t5.8k.txt", 
		"/resource/chameleon/t7.10k.txt", 
		"/resource/chameleon/t8.8k.txt", 
		// 3D Road Network
		"/resource/roadnetwork/3D_spatial_network.txt", 
		// Cloud
		"/resource/cloud/cloud.data1.txt", 
		"/resource/cloud/cloud.data2.txt", 
		// Spambase
		"/resource/spambase/spambase.data.txt"
	};
	
	private static final String[] DATA_NAMES = {
		// Chameleon
		"Chameleon(t4)", 
		"Chameleon(t5)", 
		"Chameleon(t7)", 
		"Chameleon(t8)", 
		// 3D Road Network
		"3D Road Network", 
		// Cloud
		"Cloud(1)", 
		"Cloud(2)", 
		// Spambase
		"Spambase"
	};
	
	private static final int[] DIMENSION = {
		// Chameleon
		2, 2, 2, 2, 
		// 3D Road Network
		3, 
		// Cloud
		10, 10, 
		// Spambase
		58
	};
	
	/**
	 * private constructor: prevent construction
	 */
	private Resource() {}
	
	/**
	 * get the dimension of data according to index
	 * @param index index
	 * @return the dimension
	 */
	public static int getDimension(int index) {
		return DIMENSION[index];
	}
	
	/**
	 * get the name of source data according to index 
	 * @param index index
	 * @return the name of source data
	 */
	public static String getName(int index) {
		return DATA_NAMES[index];
	}
	
	/**
	 * load data from file according to index
	 * @param index index
	 * @return data in 2-dimensional double array
	 */
	public static double[][] loadData(int index) {
		ArrayList<double[]> list = new ArrayList<double[]>();
		
		try {
			String line = null;
			String[] tuple = null;
			BufferedReader reader = new BufferedReader(new InputStreamReader(
				Resource.class.getResourceAsStream(DATA_PATHS[index])
			));
			
			while ((line = reader.readLine()) != null) {
				tuple = line.split("[\t ]+");
				double[] point = new double[DIMENSION[index]];
				for (int i = 0; i < DIMENSION[index]; i++) point[i] = Double.parseDouble(tuple[i]);
				list.add(point);
			}
			
			reader.close();
		} catch (Exception e) {
			System.err.println(e.getMessage());
			e.printStackTrace();
		}
		
		return list.isEmpty() ? null : list.toArray(new double[0][]);
	}
	
}
///:~
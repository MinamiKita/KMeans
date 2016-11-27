//: main/EuclideanSpace.java
package main;

import java.util.Arrays;
import java.util.List;

/**
 * @author chn
 * 
 * EuclideanSpace: class for providing specified method for k-means clustering algorithm
 */
public class EuclideanSpace {
	
	private int dimension = 0;
	private double squaredPrecision = 0.0;
	
	/**
	 * constructor
	 * @param dimension dimension
	 * @param precision precision for checking whether 2 point is approximately equal
	 */
	public EuclideanSpace(int dimension, double precision) {
		if (dimension < 1) 
			throw new IllegalArgumentException("Invalid parameter dimension: "+dimension);
		if (precision < 0.0) 
			throw new IllegalArgumentException("Invalid parameter precision: "+precision);
		
		this.dimension = dimension;
		this.squaredPrecision = precision * precision;
	}
	
	/**
	 * @return the dimension
	 */
	public int getDimension() {
		return dimension;
	}
	
	/**
	 * transform a point to String
	 * @param x a point in double array
	 * @return the point in String
	 */
	public String print(double[] x) {
		StringBuffer buffer = new StringBuffer("(").append(x[0]);
		for (int d = 1; d < dimension; d++) 
			buffer.append(", ").append(x[d]);
		return buffer.append(")").toString();
	}
	
	/**
	 * check whether source data is valid
	 * @param data source data in 2-dimensional double array
	 */
	public void checkData(double[][] data) {
		if (data.length < 1) 
			throw new IllegalArgumentException("Invalid length of parameter data: "+data.length);
		for (int i = 0; i < data.length; i++) {
			if (data[i].length != dimension) 
				throw new IllegalArgumentException("Invalid length of parameter data["+i+"]: "+data[i].length);
			for (int d = 0; d < dimension; d++) 
				if (Double.isNaN(data[i][d])) 
					throw new IllegalArgumentException("Invalid parameter data["+i+"]["+d+"]: "+data[i][d]);
		}
	}
	
	/**
	 * get the square of distance between 2 points
	 * @param x a point in double array
	 * @param y other point in double array
	 * @return the square of distance
	 */
	public double getSquaredDistance(double[] x, double[] y) {
		double result = 0.0, difference = 0.0;
		
		for (int d = 0; d < dimension; d++) {
			difference = x[d] - y[d];
			result += difference * difference;
		}
		
		return result;
	}
	
	/**
	 * get the centers of given points in different clusters
	 * @param k the number of clusters
	 * @param data all points
	 * @param labels clusters' labels
	 * @return the centers
	 */
	public double[][] getCenters(int k, double[][] data, int[] labels) {
		int[] size = new int[k];
		Arrays.fill(size, 0);
		double[][] centers = new double[k][dimension];
		for (int i = 0; i < k; i++) Arrays.fill(centers[i], 0.0);
		
		// adding up
		for (int i = 0; i < data.length; i++) {
			size[labels[i]]++;
			for (int d = 0; d < dimension; d++) 
				centers[labels[i]][d] += data[i][d];
		}
		
		// averaging
		for (int i = 0; i < k; i++) 
			for (int d = 0; d < dimension; d++) 
				centers[i][d] /= size[i];
		
		return centers;
	}
	
	/**
	 * get the centers of given information in buffer
	 * @param buffers buffer recording essential information from parallel Lloyd's algorithm procedure
	 * @return the centers
	 */
	public double[][] getCenters(List<ParallelKMeansProcedure.ParallelLloydsAlgorithmBuffer> buffers) {
		int[] size = new int[ParallelKMeansProcedure.ParallelLloydsAlgorithmBuffer.centers.length];
		Arrays.fill(size, 0);
		double[][] centers = new double[ParallelKMeansProcedure.ParallelLloydsAlgorithmBuffer.centers.length][dimension];
		for (int i = 0; i < ParallelKMeansProcedure.ParallelLloydsAlgorithmBuffer.centers.length; i++) 
			Arrays.fill(centers[i], 0.0);
		
		// adding up
		for (ParallelKMeansProcedure.ParallelLloydsAlgorithmBuffer buffer: buffers) {
			for (int i = 0; i < ParallelKMeansProcedure.ParallelLloydsAlgorithmBuffer.centers.length; i++) {
				size[i] += buffer.size[i];
				for (int d = 0; d < dimension; d++) 
					centers[i][d] += buffer.centroids[i][d];
			}
		}
		
		// averaging
		for (int i = 0; i < ParallelKMeansProcedure.ParallelLloydsAlgorithmBuffer.centers.length; i++) 
			for (int d = 0; d < dimension; d++) 
				centers[i][d] /= size[i];
		
		return centers;
	}
	
	/**
	 * add up assigned data and size
	 * @param buffer buffer recording assigned data and size
	 */
	public void addUp(ParallelKMeansProcedure.ParallelLloydsAlgorithmBuffer buffer) {
		buffer.size = new int[ParallelKMeansProcedure.ParallelLloydsAlgorithmBuffer.centers.length];
		Arrays.fill(buffer.size, 0);
		buffer.centroids = new double[ParallelKMeansProcedure.ParallelLloydsAlgorithmBuffer.centers.length][dimension];
		for (int i = 0; i < ParallelKMeansProcedure.ParallelLloydsAlgorithmBuffer.centers.length; i++) 
			Arrays.fill(buffer.centroids[i], 0.0);
		
		for (int i = buffer.begin; i < buffer.end; i++) {
			buffer.size[ParallelKMeansProcedure.ParallelLloydsAlgorithmBuffer.labels[i]]++;
			for (int d = 0; d < dimension; d++) 
				buffer.centroids[ParallelKMeansProcedure.ParallelLloydsAlgorithmBuffer.labels[i]][d] += 
					ParallelKMeansProcedure.ParallelLloydsAlgorithmBuffer.data[i][d];
		}
	}
	
	/**
	 * check whether a point is equal to the other
	 * @param x a point in double array
	 * @param y other point in double array
	 * @return a boolean flag
	 */
	public boolean isEqual(double[] x, double[] y) {
		for (int d = 0; d < dimension; d++) 
			if (x[d] != y[d]) 
				return false;
		return true;
	}
	
	/**
	 * check whether a point is approximately equal to the other
	 * @param x a point in double array
	 * @param y other point in double array
	 * @return a boolean flag
	 */
	public boolean isApproximatelyEqual(double[] x, double[] y) {
		if (getSquaredDistance(x, y) > squaredPrecision) return false;
		else return true;
	}
	
}
///:~
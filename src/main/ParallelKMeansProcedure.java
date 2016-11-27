//: main/ParallelKMeansProcedure.java
package main;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.Callable;

/**
 * @author chn
 * 
 * ParallelKMeansProcedure: static class for parallel k-means procedure
 */
public class ParallelKMeansProcedure {
	
	/**
	 * private constructor: prevent construction
	 */
	private ParallelKMeansProcedure() {}
	
	/**
	 * DoubleArray: class for double array
	 */
	private static class DoubleArray {
		
		private double[] array;
		
		/**
		 * constructor
		 * @param array a double array
		 */
		public DoubleArray(double[] array) {
			this.array = array;
		}
		
		/**
		 * @return the array
		 */
		public double[] getArray() {
			return array;
		}

		@Override
		public int hashCode() {
			int result = 1;
			long longBits = 0L;
			
			for (double element: array) {
				longBits = Double.doubleToLongBits(element);
				result = 31 * result + (int)(longBits ^ (longBits >>> 32));
			}
			
			return result;
		}
		
		@Override
		public boolean equals(Object object) {
			DoubleArray other = (DoubleArray) object;
			if (array.length == other.array.length) {
				for (int i = 0; i < array.length; i++) 
					if (array[i] != other.array[i]) 
						return false;
				return true;
			} else return false;
		}
		
	}
	
	/**
	 * execute parallel Lloyd's algorithm procedure
	 * @param number the number of threads
	 * @param labels clusters' labels
	 * @param data source data
	 * @param centers current centers
	 * @param space Euclidean space
	 * @return candidate centers
	 */
	public static double[][] parallelLloydsAlgorithmProcedure(
		int number, 
		int[] labels, 
		double[][] data, 
		double[][] centers, 
		EuclideanSpace space
	) {
		final int stride = data.length / number + ((data.length % number > 0) ? 1 : 0);
		List<ParallelLloydsAlgorithmThread> threads = new ArrayList<ParallelLloydsAlgorithmThread>(number);
		
		ParallelLloydsAlgorithmBuffer.labels = labels;
		ParallelLloydsAlgorithmBuffer.data = data;
		ParallelLloydsAlgorithmBuffer.centers = centers;
		ParallelLloydsAlgorithmBuffer.space = space;
		
		int begin = 0;
		final int range = data.length - stride;
		for (; begin < range; begin += stride) 
			threads.add(new ParallelLloydsAlgorithmThread(
				new ParallelLloydsAlgorithmBuffer(begin, begin+stride)
			));
		threads.add(new ParallelLloydsAlgorithmThread(
			new ParallelLloydsAlgorithmBuffer(begin, data.length)
		));
		
		List<ParallelLloydsAlgorithmBuffer> results = null;
		try {
			results = SimpleThreadPool.execute(threads);
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		return ParallelLloydsAlgorithmBuffer.space.getCenters(results);
	}
	
	/**
	 * execute parallel k-means|| algorithm procedure
	 * Reference: Bahman Bahmani et al.. Scalable KMeans++.
	 * @param l over-sampling factor
	 * @param iterations the times of iteration
	 * @param number the number of threads
	 * @param k the number of clusters
	 * @param data source data
	 * @param space Euclidean space
	 * @param random a pseudo-random number generator
	 * @return weighted candidate centers
	 * @throws Exception if the number of parallel sampling is less than k
	 */
	public static TaggedDoubleArray[] parallelKMeansVbarVbarProcedure(
		int l, 
		int iterations, 
		int number, 
		int k, 
		double[][] data, 
		EuclideanSpace space, 
		Random random
	) throws Exception {
		final int stride = data.length / number + ((data.length % number > 0) ? 1 : 0), range = data.length - stride;
		int label = 0, begin = 0;
		double[] center = data[random.nextInt(data.length)]; // uniformly random sampling
		TaggedDoubleArray[] weightedCenters = null;
		List<ParallelKMeansVbarVbarUpdatingThread> updatingThreads = new ArrayList<ParallelKMeansVbarVbarUpdatingThread>(number);
		List<ParallelKMeansVbarVbarSamplingThread> samplingThreads = new ArrayList<ParallelKMeansVbarVbarSamplingThread>(number);
		List<ParallelKMeansVbarVbarAggregatingThread> aggregatingThreads = new ArrayList<ParallelKMeansVbarVbarAggregatingThread>(number);
		Hashtable<DoubleArray, Integer> centers = new Hashtable<DoubleArray, Integer>();
		centers.put(new DoubleArray(center), Integer.valueOf(label));
		
		ParallelKMeansVbarVbarUpdatingBuffer.labels = new int[data.length];
		Arrays.fill(ParallelKMeansVbarVbarUpdatingBuffer.labels, 0);
		ParallelKMeansVbarVbarUpdatingBuffer.minimumSquaredDistance = new double[data.length];
		Arrays.fill(ParallelKMeansVbarVbarUpdatingBuffer.minimumSquaredDistance, Double.MAX_VALUE);
		ParallelKMeansVbarVbarUpdatingBuffer.data = data;
		ParallelKMeansVbarVbarUpdatingBuffer.space = space;
		ParallelKMeansVbarVbarUpdatingBuffer.newCenters = new ArrayList<TaggedDoubleArray>();
		ParallelKMeansVbarVbarUpdatingBuffer.newCenters.add(new TaggedDoubleArray(label++, center));
		
		ParallelKMeansVbarVbarSamplingBuffer.l = l;
		ParallelKMeansVbarVbarSamplingBuffer.minimumSquaredDistance = ParallelKMeansVbarVbarUpdatingBuffer.minimumSquaredDistance;
		ParallelKMeansVbarVbarSamplingBuffer.data = data;
		ParallelKMeansVbarVbarSamplingBuffer.random = random;
		
		ParallelKMeansVbarVbarAggregatingBuffer.labels = ParallelKMeansVbarVbarUpdatingBuffer.labels;
		ParallelKMeansVbarVbarAggregatingBuffer.minimumSquaredDistance = ParallelKMeansVbarVbarUpdatingBuffer.minimumSquaredDistance;
		ParallelKMeansVbarVbarAggregatingBuffer.data = ParallelKMeansVbarVbarUpdatingBuffer.data;
		ParallelKMeansVbarVbarAggregatingBuffer.space = ParallelKMeansVbarVbarUpdatingBuffer.space;
		ParallelKMeansVbarVbarAggregatingBuffer.newCenters = ParallelKMeansVbarVbarUpdatingBuffer.newCenters;
		
		for (begin = 0; begin < range; begin += stride) 
			updatingThreads.add(new ParallelKMeansVbarVbarUpdatingThread(
				new ParallelKMeansVbarVbarUpdatingBuffer(begin, begin+stride)
			));
		updatingThreads.add(new ParallelKMeansVbarVbarUpdatingThread(
			new ParallelKMeansVbarVbarUpdatingBuffer(begin, data.length)
		));
		
		for (begin = 0; begin < range; begin += stride) 
			samplingThreads.add(new ParallelKMeansVbarVbarSamplingThread(
				new ParallelKMeansVbarVbarSamplingBuffer(begin, begin+stride)
			));
		samplingThreads.add(new ParallelKMeansVbarVbarSamplingThread(
			new ParallelKMeansVbarVbarSamplingBuffer(begin, data.length)
		));
		
		for (int iteration = 0; iteration < iterations; iteration++) {
			// updating
			List<ParallelKMeansVbarVbarUpdatingBuffer> updatingResults = null;
			try {
				updatingResults = SimpleThreadPool.execute(updatingThreads);
			} catch (Exception e) {
				e.printStackTrace();
			}
			
			ParallelKMeansVbarVbarSamplingBuffer.potential = 0.0;
			for (ParallelKMeansVbarVbarUpdatingBuffer result: updatingResults) 
				ParallelKMeansVbarVbarSamplingBuffer.potential += result.cost;
			
			// sampling
			List<ParallelKMeansVbarVbarSamplingBuffer> samplingResults = null;
			try {
				samplingResults = SimpleThreadPool.execute(samplingThreads);
			} catch (Exception e) {
				e.printStackTrace();
			}
			
			ParallelKMeansVbarVbarUpdatingBuffer.newCenters.clear();
			for (ParallelKMeansVbarVbarSamplingBuffer result: samplingResults) {
				for (double[] centroid: result.centroids) {
					DoubleArray array = new DoubleArray(centroid);
					if (!centers.containsKey(array)) {
						centers.put(array, new Integer(label));
						ParallelKMeansVbarVbarUpdatingBuffer.newCenters.add(
							new TaggedDoubleArray(label++, centroid)
						);
					}
				}
			}
		}
		
		if (label < k) 
			throw new Exception("the number of samples ("+label+") is less than "+k);
		
		for (begin = 0; begin < range; begin += stride) 
			aggregatingThreads.add(new ParallelKMeansVbarVbarAggregatingThread(
				new ParallelKMeansVbarVbarAggregatingBuffer(begin, begin+stride, label)
			));
		aggregatingThreads.add(new ParallelKMeansVbarVbarAggregatingThread(
			new ParallelKMeansVbarVbarAggregatingBuffer(begin, data.length, label)
		));
		
		// aggregating
		List<ParallelKMeansVbarVbarAggregatingBuffer> aggregatingResults = null;
		try {
			aggregatingResults = SimpleThreadPool.execute(aggregatingThreads);
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		int[] sum = new int[label];
		Arrays.fill(sum, 0);
		for (ParallelKMeansVbarVbarAggregatingBuffer result: aggregatingResults) 
			for (int i = 0; i < label; i++) 
				sum[i] += result.sum[i];

		weightedCenters = new TaggedDoubleArray[label];
		for (Map.Entry<DoubleArray, Integer> entry: centers.entrySet()) {
			int index = entry.getValue().intValue();
			weightedCenters[index] = new TaggedDoubleArray(sum[index], entry.getKey().getArray());
		}
			
		return weightedCenters;
	}
	
	/**
	 * TaggedDoubleArray: class for double array and an integer tag
	 */
	public static class TaggedDoubleArray {
		
		private int tag;
		private double[] array;
		
		/**
		 * constructor
		 * @param tag an integer
		 * @param array a double array
		 */
		public TaggedDoubleArray(int tag, double[] array) {
			this.tag = tag;
			this.array = array;
		}

		/**
		 * @return the tag
		 */
		public int getTag() {
			return tag;
		}

		/**
		 * @return the array
		 */
		public double[] getArray() {
			return array;
		}
		
	}
	
	/**
	 * ParallelLloydsAlgorithmBuffer: class for buffer for ParallelLloydsAlgorithmThread
	 */
	public static class ParallelLloydsAlgorithmBuffer {
		
		public static int[] labels;
		public static double[][] data;
		public static double[][] centers;
		public static EuclideanSpace space;
		
		public int begin;
		public int end;
		public int[] size;
		public double[][] centroids;
		
		/**
		 * constructor
		 * @param begin the beginning of slice
		 * @param end the end of slice
		 */
		public ParallelLloydsAlgorithmBuffer(int begin, int end) {
			this.begin = begin;
			this.end = end;
		}
		
	}
	
	/**
	 * ParallelKMeansVbarVbarSamplingBuffer: class for buffer for ParallelKMeansVbarVbarSamplingThread
	 */
	public static class ParallelKMeansVbarVbarSamplingBuffer {
		
		public static double l;
		public static double potential;
		public static double[] minimumSquaredDistance;
		public static double[][] data;
		public static Random random;
		
		public int begin;
		public int end;
		public List<double[]> centroids;
		
		/**
		 * constructor
		 * @param begin the beginning of slice
		 * @param end the end of slice
		 */
		public ParallelKMeansVbarVbarSamplingBuffer(int begin, int end) {
			this.begin = begin;
			this.end = end;
			this.centroids = new ArrayList<double[]>();
		}
		
	}
	
	/**
	 * ParallelKMeansVbarVbarUpdatingBuffer: class for buffer for ParallelKMeansVbarVbarUpdatingThread
	 */
	public static class ParallelKMeansVbarVbarUpdatingBuffer {
		
		public static int[] labels;
		public static double[] minimumSquaredDistance;
		public static double[][] data;
		public static EuclideanSpace space;
		public static List<TaggedDoubleArray> newCenters;
		
		public int begin;
		public int end;
		public double cost;
		
		/**
		 * constructor
		 * @param begin the beginning of slice
		 * @param end the end of slice
		 */
		public ParallelKMeansVbarVbarUpdatingBuffer(int begin, int end) {
			this.begin = begin;
			this.end = end;
		}
		
	}
	
	/**
	 * ParallelKMeansVbarVbarAggregatingBuffer: class for buffer for ParallelKMeansVbarVbarAggregatingThread
	 */
	public static class ParallelKMeansVbarVbarAggregatingBuffer {
		
		public static int[] labels;
		public static double[] minimumSquaredDistance;
		public static double[][] data;
		public static EuclideanSpace space;
		public static List<TaggedDoubleArray> newCenters;
		
		public int begin;
		public int end;
		public int[] sum;
		
		/**
		 * constructor
		 * @param begin the beginning of slice
		 * @param end the end of slice
		 * @param samplesNumber the number of samples
		 */
		public ParallelKMeansVbarVbarAggregatingBuffer(int begin, int end, int samplesNumber) {
			this.begin = begin;
			this.end = end;
			this.sum = new int[samplesNumber];
		}
		
	}
	
	/**
	 * ParallelLloydsAlgorithmThread: class for thread of parallel Lloyd's algorithm
	 */
	public static class ParallelLloydsAlgorithmThread 
	implements Callable<ParallelLloydsAlgorithmBuffer> {
		
		private ParallelLloydsAlgorithmBuffer buffer;
		
		/**
		 * constructor
		 * @param buffer buffer
		 */
		public ParallelLloydsAlgorithmThread(ParallelLloydsAlgorithmBuffer buffer) {
			this.buffer = buffer;
		}
		
		@Override
		public ParallelLloydsAlgorithmBuffer call() throws Exception {
			// assignment
			for (int i = buffer.begin; i < buffer.end; i++) {
				double squaredDistance = 0.0, minimum = Double.MAX_VALUE;
				
				for (int j = 0; j < ParallelLloydsAlgorithmBuffer.centers.length; j++) {
					squaredDistance = ParallelLloydsAlgorithmBuffer.space.getSquaredDistance(
						ParallelLloydsAlgorithmBuffer.data[i], 
						ParallelLloydsAlgorithmBuffer.centers[j]
					);
					if (squaredDistance < minimum) {
						minimum = squaredDistance;
						ParallelLloydsAlgorithmBuffer.labels[i] = j;
					}
				}
			}
			
			// adding up assigned data and size
			ParallelLloydsAlgorithmBuffer.space.addUp(buffer);
			
			return buffer;
		}
		
	}
	
	/**
	 * ParallelKMeansVbarVbarSamplingThread: class for thread of the sampling phase of k-means|| algorithm
	 */
	public static class ParallelKMeansVbarVbarSamplingThread 
	implements Callable<ParallelKMeansVbarVbarSamplingBuffer> {
		
		private ParallelKMeansVbarVbarSamplingBuffer buffer;
		
		/**
		 * constructor
		 * @param buffer buffer
		 */
		public ParallelKMeansVbarVbarSamplingThread(ParallelKMeansVbarVbarSamplingBuffer buffer) {
			this.buffer = buffer;
		}
		
		@Override
		public ParallelKMeansVbarVbarSamplingBuffer call() throws Exception {
			buffer.centroids.clear();
			
			// independent sampling
			for (int i = buffer.begin; i < buffer.end; i++) 
				if (ParallelKMeansVbarVbarSamplingBuffer.random.nextDouble() * 
					ParallelKMeansVbarVbarSamplingBuffer.potential < 
					ParallelKMeansVbarVbarSamplingBuffer.l * 
					ParallelKMeansVbarVbarSamplingBuffer.minimumSquaredDistance[i]) 
					buffer.centroids.add(ParallelKMeansVbarVbarSamplingBuffer.data[i]);
			
			return buffer;
		}
		
	}
	
	/**
	 * ParallelKMeansVbarVbarUpdatingThread: class for thread of the updating phase of k-means|| algorithm
	 */
	public static class ParallelKMeansVbarVbarUpdatingThread 
	implements Callable<ParallelKMeansVbarVbarUpdatingBuffer> {
		
		private ParallelKMeansVbarVbarUpdatingBuffer buffer;
		
		/**
		 * constructor
		 * @param buffer buffer
		 */
		public ParallelKMeansVbarVbarUpdatingThread(ParallelKMeansVbarVbarUpdatingBuffer buffer) {
			this.buffer = buffer;
		}
		
		@Override
		public ParallelKMeansVbarVbarUpdatingBuffer call() throws Exception {
			double squaredDistance = 0.0;
			buffer.cost = 0.0;
			
			// (re)assignment and updating partial potential
			for (int i = buffer.begin; i < buffer.end; i++) {
				for (TaggedDoubleArray center: ParallelKMeansVbarVbarUpdatingBuffer.newCenters) {
					squaredDistance = ParallelKMeansVbarVbarUpdatingBuffer.space.getSquaredDistance(
						ParallelKMeansVbarVbarUpdatingBuffer.data[i], center.getArray()
					);
					if (squaredDistance < ParallelKMeansVbarVbarUpdatingBuffer.minimumSquaredDistance[i]) {
						ParallelKMeansVbarVbarUpdatingBuffer.minimumSquaredDistance[i] = squaredDistance;
						ParallelKMeansVbarVbarUpdatingBuffer.labels[i] = center.getTag();
					}
				}
				buffer.cost += ParallelKMeansVbarVbarUpdatingBuffer.minimumSquaredDistance[i];
			}
			
			return buffer;
		}
		
	}
	
	/**
	 * ParallelKMeansVbarVbarAggregatingThread: class for thread of the aggregating phase of k-means|| algorithm
	 */
	public static class ParallelKMeansVbarVbarAggregatingThread 
	implements Callable<ParallelKMeansVbarVbarAggregatingBuffer> {
		
		private ParallelKMeansVbarVbarAggregatingBuffer buffer;
		
		/**
		 * constructor
		 * @param buffer buffer
		 */
		public ParallelKMeansVbarVbarAggregatingThread(ParallelKMeansVbarVbarAggregatingBuffer buffer) {
			this.buffer = buffer;
		}

		@Override
		public ParallelKMeansVbarVbarAggregatingBuffer call() throws Exception {
			double squaredDistance = 0.0;
			Arrays.fill(buffer.sum, 0);
			
			// reassignment and adding up partial number of each cluster
			for (int i = buffer.begin; i < buffer.end; i++) {
				for (TaggedDoubleArray center: ParallelKMeansVbarVbarUpdatingBuffer.newCenters) {
					squaredDistance = ParallelKMeansVbarVbarUpdatingBuffer.space.getSquaredDistance(
						ParallelKMeansVbarVbarUpdatingBuffer.data[i], center.getArray()
					);
					if (squaredDistance < ParallelKMeansVbarVbarUpdatingBuffer.minimumSquaredDistance[i]) {
						ParallelKMeansVbarVbarUpdatingBuffer.minimumSquaredDistance[i] = squaredDistance;
						ParallelKMeansVbarVbarUpdatingBuffer.labels[i] = center.getTag();
					}
				}
				buffer.sum[ParallelKMeansVbarVbarUpdatingBuffer.labels[i]]++;
			}
			
			return buffer;
		}
		
	}
	
}
///:~
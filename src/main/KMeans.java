//: main/KMeans.java
package main;

import java.util.Arrays;
import java.util.Random;

/**
 * @author chn
 * 
 * KMeans: class for k-means algorithm and its variant
 */
public class KMeans {
	
	private int k;
	private int[] labels;
	private double[][] data;
	private double[][] centers;
	private EuclideanSpace space;
	private Random random;
	
	/**
	 * uniformly random sampling without dup
	 * @param range the number of selected centers
	 * @return the index of sample
	 */
	private int noDupRandomSampling(int range) {
		int index = 0;
		boolean isDup = false;
		
		do {
			index = random.nextInt(data.length);
			isDup = false;
			for (int i = 0; i < range; i++) {
				if (space.isEqual(data[index], centers[i])) {
					isDup = true;
					break;
				}
			}
		} while (isDup);
		
		return index;
	}
	
	/**
	 * D^2 weighting sampling (a kind of Roulette Wheel Selection)
	 * Reference: David Arthur, Sergei Vassilvitskii. k-means++: The Advantages of Careful Seeding.
	 * @param minimumSquaredDistance all minimum squared distance in double array
	 * @param potential the sum of all minimum squared distance
	 * @return the index of sample
	 */
	private int d2WeightingSampling(double[] minimumSquaredDistance, double potential) {
		double threshold = random.nextDouble() * potential, sum = 0.0;
		
		for (int index = 0; index < data.length; index++) {
			sum += minimumSquaredDistance[index];
			if (sum >= threshold) return index;
		}
		
		return data.length - 1;
	}
	
	/**
	 * sampling based on Markov Chain Monte Carlo (Metropolis–Hastings algorithm)
	 * Reference: Olivier Bachem et al.. Approximate K-Means++ in Sublinear Time.
	 * @param length the length of Markov chain
	 * @param range the number of selected centers
	 * @param counter all updated times in integer array
	 * @param minimumSquaredDistance all minimum squared distance in double array
	 * @return the index of sample
	 */
	private int markovChainMonteCarloSampling(
		int length, 
		int range, 
		int[] counter, 
		double[] minimumSquaredDistance
	) {
		int index = noDupRandomSampling(range), candidate = 0; // initialization
		double squaredDistance = 0.0;
		for (; counter[index] < range; counter[index]++) {
			squaredDistance = space.getSquaredDistance(data[index], centers[counter[index]]);
			if (squaredDistance < minimumSquaredDistance[index]) 
				minimumSquaredDistance[index] = squaredDistance;
		}
		
		for (int i = 1; i < length; i++) {
			candidate = random.nextInt(data.length); // random sampling
			for (; counter[candidate] < range; counter[candidate]++) {
				squaredDistance = space.getSquaredDistance(data[candidate], centers[counter[candidate]]);
				if (squaredDistance < minimumSquaredDistance[candidate]) 
					minimumSquaredDistance[candidate] = squaredDistance;
			}
			if (random.nextDouble() * minimumSquaredDistance[index] < // acceptance
				minimumSquaredDistance[candidate]) index = candidate;
		}
		
		return index;
	}
	
	/**
	 * sampling based on Markov Chain Monte Carlo (Multiple-try Metropolis)
	 * Reference: Jun S. Liu et al.. The Multiple-Try Method and Local Optimization in Metropolis Sampling.
	 * @param length the length of Markov chain
	 * @param tryNumber the number of samples in one iteration
	 * @param range the number of selected centers
	 * @param counter all updated times in integer array
	 * @param minimumSquaredDistance all minimum squared distance in double array
	 * @return the index of sample
	 */
	private int multipleTryMetropolisSampling(
		int length, 
		int tryNumber, 
		int range, 
		int[] counter, 
		double[] minimumSquaredDistance
	) {
		final int randomIndicesNumber = tryNumber - 1;
		int index = noDupRandomSampling(range); // initialization
		double squaredDistance = 0.0, candidatesSum = 0.0, indicesSum = 0.0;
		int[] candidates = new int[tryNumber], indices = new int[randomIndicesNumber];
		
		for (; counter[index] < range; counter[index]++) {
			squaredDistance = space.getSquaredDistance(data[index], centers[counter[index]]);
			if (squaredDistance < minimumSquaredDistance[index]) 
				minimumSquaredDistance[index] = squaredDistance;
		}
		
		for (int i = 1; i < length; i++) {
			candidatesSum = 0.0;
			for (int j = 0; j < tryNumber; j++) {
				candidates[j] = random.nextInt(data.length); // random sampling
				for (; counter[candidates[j]] < range; counter[candidates[j]]++) {
					squaredDistance = space.getSquaredDistance(data[candidates[j]], centers[counter[candidates[j]]]);
					if (squaredDistance < minimumSquaredDistance[candidates[j]]) 
						minimumSquaredDistance[candidates[j]] = squaredDistance;
				}
				candidatesSum += minimumSquaredDistance[candidates[j]];
			}
			
			indicesSum = minimumSquaredDistance[index];
			for (int j = 0; j < randomIndicesNumber; j++) {
				indices[j] = random.nextInt(data.length); // random history
				for (; counter[indices[j]] < range; counter[indices[j]]++) {
					squaredDistance = space.getSquaredDistance(data[indices[j]], centers[counter[indices[j]]]);
					if (squaredDistance < minimumSquaredDistance[indices[j]]) 
						minimumSquaredDistance[indices[j]] = squaredDistance;
				}
				indicesSum += minimumSquaredDistance[indices[j]];
			}
			
			if (random.nextDouble() * indicesSum < candidatesSum) { // acceptance
				candidatesSum *= random.nextDouble();
				indicesSum = 0.0;
				for (index = 0; index < tryNumber; index++) {
					indicesSum += minimumSquaredDistance[candidates[index]];
					if (indicesSum >= candidatesSum) break;
				}
				index = candidates[index];
			}
		}
		
		return index;
	}
	
	/**
	 * k-means++ algorithm for weighted points
	 * @param data weighted points in an array of TaggedDoubleArray
	 */
	private void kMeansPlusPlus(ParallelKMeansProcedure.TaggedDoubleArray[] data) {
		double potential = 0.0;
		double[] minimumSquaredDistance = new double[data.length];
		ParallelKMeansProcedure.TaggedDoubleArray[] centers = new ParallelKMeansProcedure.TaggedDoubleArray[k];
		
		// uniformly random sampling and assignment
		centers[0] = data[random.nextInt(data.length)];
		for (int i = 0; i < data.length; i++) {
			minimumSquaredDistance[i] = data[i].getTag() * 
				space.getSquaredDistance(data[i].getArray(), centers[0].getArray());
			potential += minimumSquaredDistance[i];
		}
		
		final int range = k - 1;
		for (int j = 1; j < range; j++) {
			// D^2 weighting sampling
			centers[j] = data[d2WeightingSampling(minimumSquaredDistance, potential)];
			
			// reassignment
			for (int i = 0; i < data.length; i++) {
				double squaredDistance = data[i].getTag() * 
					space.getSquaredDistance(data[i].getArray(), centers[j].getArray());
				if (squaredDistance < minimumSquaredDistance[i]) {
					potential -= minimumSquaredDistance[i] - squaredDistance;
					minimumSquaredDistance[i] = squaredDistance;
				}
			}
		}
		
		// D^2 weighting sampling
		centers[k-1] = data[d2WeightingSampling(minimumSquaredDistance, potential)];
		
		for (int j = 0; j < k; j++) 
			this.centers[j] = centers[j].getArray();
	}
	
	/**
	 * constructor
	 * @param data source data
	 * @param k the number of clusters
	 * @param space an instance of Euclidean space
	 */
	public KMeans(double[][] data, int k, EuclideanSpace space) {
		if (k < 1) throw new IllegalArgumentException("Invalid parameter k: "+k);
		space.checkData(data);
		
		this.k = k;
		this.data = data;
		this.space = space;
		
		this.labels = new int[this.data.length];
		this.centers = new double[this.k][];
		this.random = new Random();
	}
	
	/**
	 * random seeding
	 */
	public void random() {
		// uniformly random sampling
		centers[0] = data[random.nextInt(data.length)];
		
		// uniformly random sampling without dup
		for (int i = 1; i < k; i++) 
			centers[i] = data[noDupRandomSampling(i)];
	}
	
	/**
	 * seeding using k-means++ algorithm
	 * Reference: David Arthur, Sergei Vassilvitskii. k-means++: The Advantages of Careful Seeding.
	 */
	public void kMeansPlusPlus() {
		double potential = 0.0;
		double[] minimumSquaredDistance = new double[data.length];
		
		// uniformly random sampling and assignment
		centers[0] = data[random.nextInt(data.length)];
		for (int i = 0; i < data.length; i++) {
			minimumSquaredDistance[i] = space.getSquaredDistance(data[i], centers[0]);
			potential += minimumSquaredDistance[i];
		}
		
		final int range = k - 1;
		for (int j = 1; j < range; j++) {
			// D^2 weighting sampling
			centers[j] = data[d2WeightingSampling(minimumSquaredDistance, potential)];
			
			// reassignment
			for (int i = 0; i < data.length; i++) {
				double squaredDistance = space.getSquaredDistance(data[i], centers[j]);
				if (squaredDistance < minimumSquaredDistance[i]) {
					potential -= minimumSquaredDistance[i] - squaredDistance;
					minimumSquaredDistance[i] = squaredDistance;
				}
			}
		}
		
		// D^2 weighting sampling
		centers[k-1] = data[d2WeightingSampling(minimumSquaredDistance, potential)];
	}
	
	/**
	 * seeding using k-MC^2 algorithm
	 * Reference: Olivier Bachem et al.. Approximate K-Means++ in Sublinear Time.
	 * @param length the length of Markov chain
	 */
	public void kMCMC(int length) {
		/*if (length < 1) 
			throw new IllegalArgumentException("Invalid parameter length: "+length);*/
		
		int[] counter = new int[data.length];
		double[] minimumSquaredDistance = new double[data.length];
		Arrays.fill(counter, 0);
		Arrays.fill(minimumSquaredDistance, Double.MAX_VALUE);
		
		// uniformly random sampling
		centers[0] = data[random.nextInt(data.length)];
		
		// Markov Chain Monte Carlo sampling
		for (int i = 1; i < k; i++) 
			centers[i] = data[markovChainMonteCarloSampling(length, i, counter, minimumSquaredDistance)];
	}
	
	/**
	 * seeding using k-MTM algorithm
	 * @param length the length of Markov chain
	 * @param tryNumber the number of samples in one iteration
	 */
	public void kMTM(int length, int tryNumber) {
		/*if (length < 1) 
			throw new IllegalArgumentException("Invalid parameter length: "+length);
		if (tryNumber < 1) 
			throw new IllegalArgumentException("Invalid parameter tryNumber: "+tryNumber);*/
		
		int[] counter = new int[data.length];
		double[] minimumSquaredDistance = new double[data.length];
		Arrays.fill(counter, 0);
		Arrays.fill(minimumSquaredDistance, Double.MAX_VALUE);
		
		// uniformly random sampling
		centers[0] = data[random.nextInt(data.length)];
		
		// Markov Chain Monte Carlo(Multiple-try Metropolis) sampling
		for (int i = 1; i < k; i++) 
			centers[i] = data[multipleTryMetropolisSampling(length, tryNumber, i, counter, minimumSquaredDistance)];
	}
	
	/**
	 * seeding using k-means|| algorithm
	 * Reference: Bahman Bahmani et al.. Scalable KMeans++.
	 * @param l over-sampling factor
	 * @param iterations the times of iteration
	 * @param threadsNumber the number of threads
	 * @throws Exception if the number of parallel sampling is less than k
	 */
	public void kMeansVbarVbar(int l, int iterations, int threadsNumber) throws Exception {
		/*if (l < 1) 
			throw new IllegalArgumentException("Invalid parameter l: "+l);
		if (iterations < 1) 
			throw new IllegalArgumentException("Invalid parameter iterations: "+iterations);
		if (threadsNumber < 1) 
			throw new IllegalArgumentException("Invalid parameter threadsNumber: "+threadsNumber);*/
		
		kMeansPlusPlus(ParallelKMeansProcedure.parallelKMeansVbarVbarProcedure(
			l, iterations, threadsNumber, k, data, space, random
		));
	}
	
	/**
	 * clustering using Lloyd's algorithm
	 * Reference: https://en.wikipedia.org/wiki/Lloyd%27s_algorithm
	 * @param iterations the times of iteration
	 * @return actual number of iterations
	 */
	public int lloydsAlgorithm(int iterations) {
		/*if (iterations < 1) 
			throw new IllegalArgumentException("Invalid parameter iterations: "+iterations);*/
		
		boolean isNotConverged = true;
		int iteration = 0;
		
		do {
			iteration++;
			
			// assignment
			assign();
			
			// centering
			double[][] centroids = space.getCenters(k, data, labels);
			isNotConverged = false;
			for (int i = 0; i < k; i++) {
				if (!space.isApproximatelyEqual(centers[i], centroids[i])) {
					isNotConverged = true;
					break;
				}
			}
			System.arraycopy(centroids, 0, centers, 0, k);
		} while (isNotConverged && iteration < iterations);
		
		return iteration;
	}
	
	/**
	 * clustering using parallel Lloyd's algorithm
	 * @param iterations the times of iteration
	 * @param threadsNumber the number of threads
	 * @return actual number of iterations
	 */
	public int parallelLloydsAlgorithm(int iterations, int threadsNumber) {
		/*if (iterations < 1) 
			throw new IllegalArgumentException("Invalid parameter iterations: "+iterations);
		if (threadsNumber < 1) 
			throw new IllegalArgumentException("Invalid parameter threadsNumber: "+threadsNumber);*/
		
		boolean isNotConverged = true;
		int iteration = 0;
		
		do {
			iteration++;
			
			double[][] centroids = ParallelKMeansProcedure.parallelLloydsAlgorithmProcedure(
				threadsNumber, labels, data, centers, space
			);
			
			isNotConverged = false;
			for (int i = 0; i < k; i++) {
				if (!space.isApproximatelyEqual(centers[i], centroids[i])) {
					isNotConverged = true;
					break;
				}
			}
			System.arraycopy(centroids, 0, centers, 0, k);
		} while (isNotConverged && iteration < iterations);
		
		return iteration;
	}
	
	/**
	 * assign each data to its closest current centers and get labels
	 */
	public void assign() {
		for (int i = 0; i < data.length; i++) {
			double squaredDistance = 0.0, minimum = Double.MAX_VALUE;
			
			for (int j = 0; j < k; j++) {
				squaredDistance = space.getSquaredDistance(data[i], centers[j]);
				if (squaredDistance < minimum) {
					minimum = squaredDistance;
					labels[i] = j;
				}
			}
		}
	}
	
	/**
	 * get the potential (or cost, or the sum of all minimum squared distance) of data
	 * @return the potential
	 */
	public double getPotential() {
		double potential = 0.0;
		
		for (int i = 0; i < data.length; i++) 
			potential += space.getSquaredDistance(data[i], centers[labels[i]]);
		
		return potential;
	}
	
	/**
	 * get the size of the set of source data
	 * @return the length of data
	 */
	public int getDataLength() {
		return data.length;
	}
	
	@Override
	public String toString() {
		StringBuilder builder = new StringBuilder();
        
		//builder.append("k-means: divided "+data.length+" data into "+k+" clusters\n");
		builder.append("centers: \n");
		for (int i = 0; i < k; i++) 
			builder.append(space.print(centers[i])).append('\n');
		/*
		builder.append("labels: \n");
		for (int i = 0; i < data.length; i++) 
			builder.append(space.print(data[i])).append(" -> ").append(labels[i]).append('\n');
        */
        return builder.toString();
	}
	
	/**
	 * for test
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		KMeans kmeans = new KMeans(
			new double[][]{{1.0, 1.0}, {1.1, 1.0}, {0.9, 1.0}, {1.0, 1.1}, {1.1, 1.1}, {1.0, 0.9}, {2.0, 2.0}, {2.1, 2.0}, {2.0, 2.1}, {2.1, 2.1}}, 
			2, 
			new EuclideanSpace(2, 1e-3)
		);
		/*
		for (int i = 0; i < 10; i++) {
			kmeans.kMeansVbarVbar(2, 2, SimpleThreadPool.getAvailableProcessorsNumber());
			System.out.println(kmeans.toString());
		}
		*/
		
		kmeans.kMeansVbarVbar(2, 2, SimpleThreadPool.getAvailableProcessorsNumber());
		kmeans.assign();
		System.out.println(kmeans.getPotential());
		System.out.println(kmeans.toString());
		kmeans.lloydsAlgorithm(100);
		System.out.println(kmeans.getPotential());
		System.out.println(kmeans.toString());
		
	}

}
///:~
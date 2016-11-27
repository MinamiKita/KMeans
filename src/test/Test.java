//: test/Test.java
package test;

import java.io.FileWriter;
import java.io.IOException;

import main.EuclideanSpace;
import main.KMeans;
import main.SimpleThreadPool;

/**
 * @author chn
 * 
 * Test: class for test
 */
public class Test {
	
	/**
	 * private constructor: prevent construction
	 */
	private Test() {}
	
	/**
	 * get the mean of a set of value in int array
	 * @param array a set of value
	 * @return the mean
	 */
	private static double getMean(int[] array) {
		int mean = 0;
		
		for (int value: array) mean += value;
		
		return (double) mean / array.length;
	}
	
	/**
	 * get the mean of a set of value in double array
	 * @param array a set of value
	 * @return the mean
	 */
	private static double getMean(double[] array) {
		double mean = 0.0;
		
		for (double value: array) mean += value;
		
		return mean / array.length;
	}
	
	/**
	 * get the mean of a set of value in long array
	 * @param array a set of value
	 * @return the mean
	 */
	private static double getMean(long[] array) {
		long mean = 0L;
		
		for (long value: array) mean += value;
		
		return (double) mean / array.length;
	}
	
	/**
	 * get the standard error of the mean(SEM) of a set of value in double array
	 * Reference: https://en.wikipedia.org/wiki/Standard_error
	 * @param array a set of value
	 * @param mean the mean
	 * @return the standard error of the mean(SEM)
	 */
	private static double getSEM(double[] array, double mean) {
		double sem = 0.0, difference = 0.0;
		
		for (double value: array) {
			difference = value - mean;
			sem += difference * difference;
		}
		
		return Math.sqrt(sem/array.length/(array.length-1));
	}
	
	/**
	 * write information into csv file
	 * @param writer FileWriter instance
	 * @param objects information
	 * @throws IOException
	 */
	private static void write(FileWriter writer, Object ... objects) throws IOException {
		if (objects.length > 0) {
			StringBuffer buffer = new StringBuffer();
			
			if (objects[0] instanceof Double) buffer.append(String.format("%.4f", objects[0]));
			else buffer.append(objects[0].toString());
			for (int i = 1; i < objects.length; i++) {
				if (objects[i] instanceof Double) buffer.append(String.format(",%.4f", objects[i]));
				else buffer.append(',').append(objects[i].toString());
			}
			
			writer.write(buffer.toString());
		}
	}
	
	/**
	 * write a row of information into csv file
	 * @param writer FileWriter instance
	 * @param objects information
	 * @throws IOException
	 */
	private static void writeLine(FileWriter writer, Object ... objects) throws IOException {
		if (objects.length > 0) {
			StringBuffer buffer = new StringBuffer();
			
			if (objects[0] instanceof Double) buffer.append(String.format("%.4f", objects[0]));
			else buffer.append(objects[0].toString());
			for (int i = 1; i < objects.length; i++) {
				if (objects[i] instanceof Double) buffer.append(String.format(",%.4f", objects[i]));
				else buffer.append(',').append(objects[i].toString());
			}
			
			writer.write(buffer.append('\n').toString());
		} else writer.write('\n');
	}
	
	/**
	 * append information into csv file
	 * @param writer FileWriter instance
	 * @param objects information
	 * @throws IOException
	 */
	private static void append(FileWriter writer, Object ... objects) throws IOException {
		if (objects.length > 0) {
			StringBuffer buffer = new StringBuffer();
			
			for (Object object: objects) {
				if (object instanceof Double) buffer.append(String.format(",%.4f", object));
				else buffer.append(',').append(object.toString());
			}
			
			writer.write(buffer.toString());
		}
	}
	
	/**
	 * append a row of information into csv file
	 * @param writer FileWriter instance
	 * @param objects information
	 * @throws IOException
	 */
	private static void appendLine(FileWriter writer, Object ... objects) throws IOException {
		if (objects.length > 0) {
			StringBuffer buffer = new StringBuffer();
			
			for (Object object: objects) {
				if (object instanceof Double) buffer.append(String.format(",%.4f", object));
				else buffer.append(',').append(object.toString());
			}
			
			writer.write(buffer.append('\n').toString());
		} else writer.write('\n');
	}
	
	/**
	 * test (seeding and main procedure; potential and time)
	 * @param writer a FileWriter instance
	 * @param index data's index
	 * @param times repeated times
	 * @param k the number of clusters
	 * @param precision precision for checking whether 2 point is approximately equal
	 * @param length the length of Markov chain (k-MC^2, k-MTM)
	 * @param tryNumber the number of samples in one iteration (k-MTM)
	 * @param l over-sampling factor (k-means||)
	 * @param iterations the times of iteration (k-means||)
	 * @param threadsNumber the number of threads (k-means||, parallel Lloyd's algorithm)
	 * @throws IOException
	 */
	private static void test(
		FileWriter writer, 
		int index, 
		int times, 
		int k, 
		double precision, 
		int[] length, 
		int[] tryNumber, 
		int[] l, 
		int[] iterations, 
		int threadsNumber
	) throws IOException {
		double mean = 0.0, halfInterval = 0.0;
		int[] iteration = new int[times];
		double[] potential = new double[times];
		long[] time = new long[times];
		KMeans kmeans = new KMeans(
			Resource.loadData(index), 
			k, 
			new EuclideanSpace(Resource.getDimension(index), precision)
		);
		
		// ==================== random ====================
		System.out.println("random ... ");
		// test seeding
		for (int i = 0; i < times; i++) {
			time[i] = System.currentTimeMillis();
			kmeans.random();
			time[i] = System.currentTimeMillis() - time[i];
			kmeans.assign();
			potential[i] = kmeans.getPotential();
		}
		
		mean = getMean(potential);
		halfInterval = 1.96 * getSEM(potential, mean);
		write(writer, "random", 0, getMean(time), mean, halfInterval);
		
		// test main procedure (serial)
		for (int i = 0; i < times; i++) {
			kmeans.random();
			time[i] = System.currentTimeMillis();
			iteration[i] = kmeans.lloydsAlgorithm(Integer.MAX_VALUE);
			time[i] = System.currentTimeMillis() - time[i];
			potential[i] = kmeans.getPotential();
		}
		
		mean = getMean(potential);
		halfInterval = 1.96 * getSEM(potential, mean);
		append(writer, getMean(iteration), getMean(time), mean, halfInterval);
		
		// test main procedure (parallel)
		for (int i = 0; i < times; i++) {
			kmeans.random();
			time[i] = System.currentTimeMillis();
			iteration[i] = kmeans.parallelLloydsAlgorithm(Integer.MAX_VALUE, threadsNumber);
			time[i] = System.currentTimeMillis() - time[i];
			potential[i] = kmeans.getPotential();
		}
		
		mean = getMean(potential);
		halfInterval = 1.96 * getSEM(potential, mean);
		appendLine(writer, getMean(iteration), getMean(time), mean, halfInterval);
		// ==================== random ====================
		
		// ==================== k-means++ ====================
		System.out.println("k-means++ ... ");
		// test seeding
		for (int i = 0; i < times; i++) {
			time[i] = System.currentTimeMillis();
			kmeans.kMeansPlusPlus();
			time[i] = System.currentTimeMillis() - time[i];
			kmeans.assign();
			potential[i] = kmeans.getPotential();
		}
		
		mean = getMean(potential);
		halfInterval = 1.96 * getSEM(potential, mean);
		write(writer, "k-means++", k*kmeans.getDataLength(), getMean(time), mean, halfInterval);
		
		// test main procedure (serial)
		for (int i = 0; i < times; i++) {
			kmeans.kMeansPlusPlus();
			time[i] = System.currentTimeMillis();
			iteration[i] = kmeans.lloydsAlgorithm(Integer.MAX_VALUE);
			time[i] = System.currentTimeMillis() - time[i];
			potential[i] = kmeans.getPotential();
		}
		
		mean = getMean(potential);
		halfInterval = 1.96 * getSEM(potential, mean);
		append(writer, getMean(iteration), getMean(time), mean, halfInterval);
		
		// test main procedure (parallel)
		for (int i = 0; i < times; i++) {
			kmeans.kMeansPlusPlus();
			time[i] = System.currentTimeMillis();
			iteration[i] = kmeans.parallelLloydsAlgorithm(Integer.MAX_VALUE, threadsNumber);
			time[i] = System.currentTimeMillis() - time[i];
			potential[i] = kmeans.getPotential();
		}
		
		mean = getMean(potential);
		halfInterval = 1.96 * getSEM(potential, mean);
		appendLine(writer, getMean(iteration), getMean(time), mean, halfInterval);
		// ==================== k-means++ ====================
		
		// ==================== k-MC^2 ====================
		for (int j = 0; j < length.length; j++) {
			System.out.println("k-MC^2(m="+length[j]+") ... ");
			// test seeding
			for (int i = 0; i < times; i++) {
				time[i] = System.currentTimeMillis();
				kmeans.kMCMC(length[j]);
				time[i] = System.currentTimeMillis() - time[i];
				kmeans.assign();
				potential[i] = kmeans.getPotential();
			}
			
			mean = getMean(potential);
			halfInterval = 1.96 * getSEM(potential, mean);
			write(writer, "k-MC^2(m="+length[j]+")", length[j]*k*k, getMean(time), mean, halfInterval);
			
			// test main procedure (serial)
			for (int i = 0; i < times; i++) {
				kmeans.kMCMC(length[j]);
				time[i] = System.currentTimeMillis();
				iteration[i] = kmeans.lloydsAlgorithm(Integer.MAX_VALUE);
				time[i] = System.currentTimeMillis() - time[i];
				potential[i] = kmeans.getPotential();
			}
			
			mean = getMean(potential);
			halfInterval = 1.96 * getSEM(potential, mean);
			append(writer, getMean(iteration), getMean(time), mean, halfInterval);
			
			// test main procedure (parallel)
			for (int i = 0; i < times; i++) {
				kmeans.kMCMC(length[j]);
				time[i] = System.currentTimeMillis();
				iteration[i] = kmeans.parallelLloydsAlgorithm(Integer.MAX_VALUE, threadsNumber);
				time[i] = System.currentTimeMillis() - time[i];
				potential[i] = kmeans.getPotential();
			}
			
			mean = getMean(potential);
			halfInterval = 1.96 * getSEM(potential, mean);
			appendLine(writer, getMean(iteration), getMean(time), mean, halfInterval);
		}
		// ==================== k-MC^2 ====================
		
		// ==================== k-MTM ====================
		for (int m = 0; m < tryNumber.length; m++) {
			for (int j = 0; j < length.length; j++) {
				System.out.println("k-MTM(m="+length[j]+";k'="+tryNumber[m]+") ... ");
				// test seeding
				for (int i = 0; i < times; i++) {
					time[i] = System.currentTimeMillis();
					kmeans.kMTM(length[j], tryNumber[m]);
					time[i] = System.currentTimeMillis() - time[i];
					kmeans.assign();
					potential[i] = kmeans.getPotential();
				}
				
				mean = getMean(potential);
				halfInterval = 1.96 * getSEM(potential, mean);
				write(
					writer, "k-MTM(m="+length[j]+";k'="+tryNumber[m]+")", 
					length[j]*tryNumber[m]*k*k, getMean(time), mean, halfInterval
				);
				
				// test main procedure (serial)
				for (int i = 0; i < times; i++) {
					kmeans.kMTM(length[j], tryNumber[m]);
					time[i] = System.currentTimeMillis();
					iteration[i] = kmeans.lloydsAlgorithm(Integer.MAX_VALUE);
					time[i] = System.currentTimeMillis() - time[i];
					potential[i] = kmeans.getPotential();
				}
				
				mean = getMean(potential);
				halfInterval = 1.96 * getSEM(potential, mean);
				append(writer, getMean(iteration), getMean(time), mean, halfInterval);
				
				// test main procedure (parallel)
				for (int i = 0; i < times; i++) {
					kmeans.kMTM(length[j], tryNumber[m]);
					time[i] = System.currentTimeMillis();
					iteration[i] = kmeans.parallelLloydsAlgorithm(Integer.MAX_VALUE, threadsNumber);
					time[i] = System.currentTimeMillis() - time[i];
					potential[i] = kmeans.getPotential();
				}
				
				mean = getMean(potential);
				halfInterval = 1.96 * getSEM(potential, mean);
				appendLine(writer, getMean(iteration), getMean(time), mean, halfInterval);
			}
		}
		// ==================== k-MTM ====================
		
		// ==================== k-means|| ====================
		for (int m = 0; m < l.length; m++) {
			for (int j = 0; j < iterations.length; j++) {
				System.out.println("k-means||(l="+l[m]+";r="+iterations[j]+") ... ");
				// test seeding
				int i = 0, failure = 0;
				while (i < times) {
					try {
						time[i] = System.currentTimeMillis();
						kmeans.kMeansVbarVbar(l[m], iterations[j], threadsNumber);
						time[i] = System.currentTimeMillis() - time[i];
						kmeans.assign();
						potential[i] = kmeans.getPotential();
					} catch (Exception e) {
						failure++;
						continue;
					}
					i++;
				}
				
				mean = getMean(potential);
				halfInterval = 1.96 * getSEM(potential, mean);
				write(
					writer, "k-means||(l="+l[m]+";r="+iterations[j]+")", 
					l[m]*iterations[j]*kmeans.getDataLength(), 
					getMean(time), mean, halfInterval
				);
				
				// test main procedure (serial)
				i = 0;
				while (i < times) {
					try {
						kmeans.kMeansVbarVbar(l[m], iterations[j], threadsNumber);
						time[i] = System.currentTimeMillis();
						iteration[i] = kmeans.lloydsAlgorithm(Integer.MAX_VALUE);
						time[i] = System.currentTimeMillis() - time[i];
						potential[i] = kmeans.getPotential();
					} catch (Exception e) {
						failure++;
						continue;
					}
					i++;
				}
				
				mean = getMean(potential);
				halfInterval = 1.96 * getSEM(potential, mean);
				append(writer, getMean(iteration), getMean(time), mean, halfInterval);
				
				// test main procedure (parallel)
				i = 0;
				while (i < times) {
					try {
						kmeans.kMeansVbarVbar(l[m], iterations[j], threadsNumber);
						time[i] = System.currentTimeMillis();
						iteration[i] = kmeans.parallelLloydsAlgorithm(Integer.MAX_VALUE, threadsNumber);
						time[i] = System.currentTimeMillis() - time[i];
						potential[i] = kmeans.getPotential();
					} catch (Exception e) {
						failure++;
						continue;
					}
					i++;
				}
				
				mean = getMean(potential);
				halfInterval = 1.96 * getSEM(potential, mean);
				appendLine(writer, getMean(iteration), getMean(time), mean, halfInterval, failure/3.0);
			}
		}
		// ==================== k-means|| ====================
	}
	
	/**
	 * entry function for test
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		final int times = 1000, threadsNumber = SimpleThreadPool.getAvailableProcessorsNumber();
		final double precision = 1e-4;
		final int[] indices = {2, 4, 5, 7}, k = {100, 10, 10, 40};
		final int[] length = {1, 2, 4, 8, 16, 32, 64, 128}, tryNumber = {1, 2, 4, 8, 16};
		final int[] iterations = {2, 3, 4, 5, 6, 7, 8};
		final int[][] l = { // 0.4k, 0.5k, 0.7k, k, 1.5k
			{40, 50, 70, 100, 150}, 
			{4, 5, 7, 10, 15}, 
			{4, 5, 7, 10, 15}, 
			{16, 20, 28, 40, 60}
		};
		FileWriter writer = new FileWriter("test_result.csv");
		
		writeLine(writer, "global parameters");
		writeLine(writer, "repeated times", "threads' number", "precision");
		writeLine(writer, times, threadsNumber, precision);
		writeLine(writer);
		writeLine(
			writer, "method", 
			"#distance evaluation", "avg. elapsed time(ms)", "potential", "half confidence interval", 
			"avg. #iteration", "avg. elapsed time(ms)", "potential", "half confidence interval", 
			"avg. #iteration", "avg. elapsed time(ms)", "potential", "half confidence interval", 
			"avg. #failure"
		);
		
		for (int i = 0; i < indices.length; i++) {
			System.out.println((i+1)+") "+Resource.getName(indices[i]));
			writeLine(writer, Resource.getName(indices[i]), "dimension="+Resource.getDimension(indices[i]), "k="+k[i]);
			test(writer, indices[i], times, k[i], precision, length, tryNumber, l[i], iterations, threadsNumber);
			writer.flush();
		}
		
		writer.close();
	}
	
}
///:~
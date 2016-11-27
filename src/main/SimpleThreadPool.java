//: main/SimpleThreadPool.java
package main;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadFactory;
import java.util.concurrent.ThreadPoolExecutor;

/**
 * @author chn
 * 
 * SimpleThreadPool: static class for a simple fixed-size thread pool
 */
public class SimpleThreadPool {
	
	private static final int processorsNumber;
	private static ThreadPoolExecutor executor = null;
	
	static {
		processorsNumber = Runtime.getRuntime().availableProcessors();
		if (processorsNumber > 1) 
			executor = (ThreadPoolExecutor) Executors.newFixedThreadPool(
				processorsNumber, 
				new ThreadFactory() {
					
					@Override
					public Thread newThread(Runnable r) {
						Thread thread = new Thread(r);
						thread.setDaemon(true); // !
						return thread;
					}
					
				}
			);
	}
	
	/**
	 * private constructor: prevent construction
	 */
	private SimpleThreadPool() {}
	
	/**
	 * get the number of available processors
	 * @return the number of available processors
	 */
	public static int getAvailableProcessorsNumber() {
		return processorsNumber;
	}
	
	/**
	 * execute a collection of tasks concurrently if the local machine has available multi-processors
	 * if not, run these tasks serially
	 * @param tasks a collection of tasks
	 * @return a list of results of tasks
	 * @throws Exception if executing tasks or getting results has something wrong
	 */
	public static <Type> List<Type> execute(Collection<? extends Callable<Type>> tasks) 
		throws Exception {
		List<Type> results = new ArrayList<Type>();
		
		if (executor != null && executor.getActiveCount() < processorsNumber) {
			List<Future<Type>> futures = executor.invokeAll(tasks);
			for (Future<Type> future: futures) 
				results.add(future.get());
		} else {
			for (Callable<Type> task: tasks) 
				results.add(task.call());
		}
		
		return results;
	}
	
	/**
	 * shutdown the thread pool
	 */
	public static void shutdown() {
		if (executor != null) executor.shutdown();
	}
	
}
///:~
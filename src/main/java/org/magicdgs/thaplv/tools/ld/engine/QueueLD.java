/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2015 Daniel Gómez-Sánchez
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

package org.magicdgs.thaplv.tools.ld.engine;

import org.magicdgs.thaplv.haplotypes.light.LightGenotype;
import org.magicdgs.thaplv.haplotypes.light.SNPpair;
import org.magicdgs.thaplv.utils.concurrent.OtherExecutors;
import org.magicdgs.thaplv.utils.stats.popgen.LDfunctions;

import htsjdk.samtools.util.Log;
import htsjdk.variant.variantcontext.VariantContext;

import java.io.Closeable;
import java.io.IOException;
import java.util.Deque;
import java.util.LinkedList;
import java.util.concurrent.ThreadPoolExecutor;

/**
 * Class to use multi-thread processing of the pairs
 * Only works for binning now
 * TODO: create an LD writer that sort on the fly the SNP pairs to allow the multi-processing to
 * TODO: write the results
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class QueueLD implements Closeable {

    private final static int BUFFER_SIZE = 500_000;

    // Executor for the thread pool
    private final ThreadPoolExecutor executor;
    // queue with the light genotypes
    private final Deque<LightGenotype> snpQueue;
    // logger for this class
    private final Log logger;
    // logger for the variants
    private final VariantLogger progress;

    // current contig
    private String currentContig = null;
    // bins
    private final LDbinning ldBins;
    // bin Writer
    private BinLDwriter output;

    // PARAMETERS
    final int minimumDistance, maximumDistance, minimumSamples;
    final boolean rmSingletons;
    final double chiSqrQuantile;

    // TODO: check the parameters here instead of in the tool class
    public QueueLD(String prefix, int binLength, int minimumDistance, int maximumDistance,
            int minimumSamples, boolean rmSingletons, double chiSqrQuantile, int nThreads,
            Log logger) {
        if (binLength < 1) {
            throw new IllegalArgumentException("Bin length cannot be smaller than 1");
        }
        if (maximumDistance > 0 && minimumDistance > maximumDistance) {
            throw new IllegalArgumentException(
                    "Minimum distance between markers cannot be bigger than maximum distance");
        } else if (minimumDistance < 0) {
            throw new IllegalArgumentException(
                    "Minimum distance between markers cannot be negative");
        }
        if (minimumSamples < 1) {
            throw new IllegalArgumentException("Minimum samples cannot be smaller than 1");
        }
        if (chiSqrQuantile < 0 || chiSqrQuantile >= 1) {
            throw new IllegalArgumentException("Chi-square quantile should be in the range (0, 1)");
        }
        if (nThreads < 1) {
            throw new IllegalArgumentException("Threads cannot be less than 1");
        }
        this.minimumDistance = minimumDistance;
        this.maximumDistance = maximumDistance;
        this.chiSqrQuantile = chiSqrQuantile;
        this.minimumSamples = minimumSamples;
        this.rmSingletons = rmSingletons;
        // get a logger for this class
        this.logger = Log.getInstance(this.getClass());
        // start with a certain number of threads
        if (nThreads < 2) {
            this.logger.debug("Less than 2 threads were provided; using a single thread QueueLD");
            this.executor = null;
        } else {
            // don't waste resources when generating threads; if the user provides more than the buffer size it's his own risk
            int bufferSize = (BUFFER_SIZE < nThreads) ? nThreads : BUFFER_SIZE;
            // create an executor that is bounded in the buffer size
            this.logger.debug("Using ", nThreads, " threads with buffer-size=", bufferSize);
            this.executor = (ThreadPoolExecutor) OtherExecutors
                    .newFixedThreadPoolWithBoundedQueue(nThreads, bufferSize);
        }
        // start the logger
        // TODO: change how many times the logger is triggered
        this.progress = new VariantLogger(logger, 1_000_000, "Computed LD on", "pairs");
        // I don't think that I need a ConcurrentList with the current implemented idea
        this.snpQueue = new LinkedList<>();
        // start the bins
        this.ldBins = new LDbinning(binLength);
        try {
            this.output = new BinLDwriter(prefix, ldBins);
        } catch (IOException e) {
            throw new IllegalArgumentException(e);
        }
    }

    /**
     * Add the variant to the queue
     *
     * @param variant the variant to add
     */
    public boolean add(VariantContext variant) throws IOException {
        // get the contig
        String contig = variant.getContig();
        if (currentContig == null) {
            currentContig = contig;
        }
        // check if the queue is empty
        if (snpQueue.isEmpty() ||
                // or if it is n the contig and in the first variant, add directly
                (currentContig.equals(contig) && (maximumDistance < 1 || (variant.getStart()
                        <= snpQueue.peek().getPosition() + maximumDistance)))) {
            // directly added
            return snpQueue.add(new LightGenotype(variant));
        }
        // if it is not, check if it change the chromosome
        if (snpQueue.peekLast().getContig().equals(contig)) {
            // compute only for the firs variant
            computeQueueLD();
            // recursive call
            return add(variant);
        } else {
            // empty the queue
            while (!snpQueue.isEmpty()) {
                computeQueueLD();
            }
            // write the bins
            writeBins();
            // we changed the contig
            currentContig = contig;
            // add variant to the queue
            return snpQueue.add(new LightGenotype(variant));
        }
    }

    /**
     * Compute LD for the first element in the queue against all the other
     */
    public void computeQueueLD() {
        if (executor == null) {
            singleThreadComputation();
        } else {
            multiThreadComputation();
        }
    }

    /**
     * Process the queue in a multi-thread way
     */
    private void multiThreadComputation() {
        // extract the first variant
        LightGenotype firstVariant = snpQueue.remove();
        // for each of the other variants
        for (LightGenotype second : snpQueue) {
            // run the thread to compute LD in the ThreadPool
            executor.execute(new LDcomputation(firstVariant, second));
        }
    }

    /**
     * Process the queue in a sinlge thread
     */
    private void singleThreadComputation() {
        // extract the first variant
        LightGenotype firstVariant = snpQueue.remove();
        // construct the SNP pair for the rest of the variants
        // TODO: filter by maf
        for (LightGenotype second : snpQueue) {
            LDcomputation toCompute = new LDcomputation(firstVariant, second);
            toCompute.run();
        }
    }

    /**
     * Write the bins and reset them
     */
    private void writeBins() throws IOException {
        // wait till completion
        if (executor != null) {
            waitUntilFinish();
        }
        output.write(currentContig, ldBins);
        ldBins.clear();
    }

    /**
     * Wait for the executor to finish all the loaded task
     */
    private void waitUntilFinish() {
        // TODO: check if this works properly
        while (!executor.getQueue().isEmpty()) {

        }
        monitorLogging();
    }

    /**
     * Return the number of records in memory
     */
    public int size() {
        return snpQueue.size();
    }

    /**
     * Log the class with a monitor
     */
    private synchronized void monitorLogging() {
        if (executor != null) {
            logger.debug(String
                    .format("[thread-pool] [%d/%d] Active: %d, Completed: %d, Task: %d, isShutdown: %s, isTerminated: %s",
                            executor.getPoolSize(), executor.getCorePoolSize(),
                            executor.getActiveCount(),
                            executor.getCompletedTaskCount(), executor.getTaskCount(),
                            executor.isShutdown(),
                            executor.isTerminated()));
        }
    }

    /**
     * Close the queue and finish computation
     */
    @Override
    public void close() throws IOException {
        logger.debug("Computing QueueLD when close");
        // empty the queue
        while (!snpQueue.isEmpty()) {
            computeQueueLD();
        }
        if (executor != null) {
            logger.debug("Shutdown the executor");
            executor.shutdown();
        }
        writeBins();
        progress.logNumberOfVariantsProcessed();
        output.close();
        monitorLogging();
    }

    /**
     * Runnable class for compute SNP pair
     */
    private class LDcomputation implements Runnable {

        private final LightGenotype geno1, geno2;

        public LDcomputation(LightGenotype genotype1, LightGenotype genotype2) {
            this.geno1 = genotype1;
            this.geno2 = genotype2;
        }

        @Override
        public void run() {
            // only generate the pair if the minimum distance threshold is hold; may be faster than generate it
            if (minimumDistance == 0 || !(Math.abs(geno1.getPosition() - geno2.getPosition())
                    < minimumDistance)) {
                // generate the SNP pair
                SNPpair pair = new SNPpair(geno1, geno2);
                // check if the conditions are met
                if (!(pair.getTotalCounts() < minimumSamples
                        // total counts are lower than the minimum samples
                        || pair.isInvariantA()                        // A is invariant
                        || pair.isInvariantB()                        // B is invariant
                        || (rmSingletons && pair
                        .oneIsSingleton())))    // rmSingletons is set and the pair contains a singleton
                {
                    // compute LD statistic for the pair
                    // logger.debug("Computing LD for pair in a distance of ", pair.getDistance());
                    double[] LD = LDfunctions.rStatistics(pair, chiSqrQuantile);
                    // add to the binning if computed
                    if (LD != null) {
                        ldBins.add(pair, LD);
                        // log the pair if it is computed
                        progress.variants(pair.getReferenceA(), pair.getPositionA(),
                                "first variant");
                    }
                }
            }
        }
    }

}

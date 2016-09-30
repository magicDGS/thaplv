/*
 * Copyright (c) 2016, Daniel Gomez-Sanchez <daniel.gomez.sanchez@hotmail> All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted
 * provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this list of conditions
 * and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions
 * and the following disclaimer in the documentation and/or other materials provided with the
 * distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse
 * or promote products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
 * FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
 * IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
 * THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

package org.magicdgs.thaplv.tools.ld.engine;

import org.magicdgs.thaplv.haplotypes.light.LightGenotype;
import org.magicdgs.thaplv.haplotypes.light.SNPpair;
import org.magicdgs.thaplv.utils.concurrent.ExecutorsFactory;
import org.magicdgs.thaplv.utils.stats.popgen.LDfunctions;

import htsjdk.samtools.util.AbstractProgressLogger;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

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

    // logger for this class
    private static final Logger logger = LogManager.getLogger(QueueLD.class);


    // Executor for the thread pool
    private final ThreadPoolExecutor executor;
    // queue with the light genotypes
    private final Deque<LightGenotype> snpQueue;
    // logger for the variants
    // TODO: probably change for ProgressMeter
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

    public QueueLD(String prefix, int binLength, int minimumDistance, int maximumDistance,
            int minimumSamples, boolean rmSingletons, double chiSqrQuantile, int nThreads,
            Logger logger) {
        // TODO: make argument validation here too?
        this.minimumDistance = minimumDistance;
        this.maximumDistance = maximumDistance;
        this.chiSqrQuantile = chiSqrQuantile;
        this.minimumSamples = minimumSamples;
        this.rmSingletons = rmSingletons;
        // start with a certain number of threads
        if (nThreads < 2) {
            this.logger.debug("Less than 2 threads were provided; using a single thread QueueLD");
            this.executor = null;
        } else {
            // don't waste resources when generating threads; if the user provides more than the buffer size it's his own risk
            int bufferSize = (BUFFER_SIZE < nThreads) ? nThreads : BUFFER_SIZE;
            // create an executor that is bounded in the buffer size
            this.logger.debug("Using {} threads with buffer-size={}", nThreads, bufferSize);
            this.executor = (ThreadPoolExecutor) ExecutorsFactory
                    .newFixedThreadPoolWithBoundedQueue(nThreads, bufferSize);
        }
        // start the logger
        this.progress = new VariantLogger(logger);
        // TODO: I don't think that I need a ConcurrentList with the current implemented idea
        this.snpQueue = new LinkedList<>();
        // start the bins
        this.ldBins = new LDbinning(binLength);
        try {
            this.output = new BinLDwriter(prefix, ldBins);
        } catch (IOException e) {
            // TODO: change this exception handling
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
            logger.debug(
                    "[thread-pool] [{}/{}] Active: {}, Completed: {}, Task: {}, isShutdown: {}, isTerminated: {}",
                    () -> new Object[] {
                            executor.getPoolSize(), executor.getCorePoolSize(),
                            executor.getActiveCount(),
                            executor.getCompletedTaskCount(), executor.getTaskCount(),
                            executor.isShutdown(),
                            executor.isTerminated()
                    });
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
        // TODO: logging all pairs processed by this queue
        // progress.logNumberOfVariantsProcessed();
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
                        logger.debug("First pair: {}", pair);
                        progress.record(pair.getReferenceA(), pair.getPositionA());
                    }
                }
            }
        }
    }

    private static class VariantLogger extends AbstractProgressLogger {

        private final Logger variantLog;

        // TODO: javadoc
        protected VariantLogger(final Logger variantLog) {
            super("pairs", "Computed LD on", 1_000_000);
            this.variantLog = variantLog;
        }

        @Override
        protected void log(String... message) {
            StringBuilder stringBuilder = new StringBuilder();
            for (String s : message) {
                stringBuilder.append(s);
            }
            variantLog.info(stringBuilder.toString());
        }
    }

}

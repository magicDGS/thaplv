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

import org.magicdgs.thaplv.cmd.argumentcollections.LengthBinningArgumentCollection;
import org.magicdgs.thaplv.cmd.argumentcollections.MultiThreadComputationArgumentCollection;
import org.magicdgs.thaplv.haplotypes.light.LightGenotype;
import org.magicdgs.thaplv.haplotypes.light.SNPpair;
import org.magicdgs.thaplv.utils.stats.popgen.LDfunctions;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.Closeable;
import java.util.Deque;
import java.util.LinkedList;
import java.util.Optional;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.Consumer;

/**
 * Class for computing LD using a queue of variants, keep the necessary ones.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 * @deprecated the engine for ld will be re-implemented from scratch
 */
@Deprecated
public class QueueLD implements Closeable {

    // logger for this class
    private static final Logger logger = LogManager.getLogger(QueueLD.class);

    private final AtomicInteger computedPairs;
    private final AtomicInteger addedPairs;

    // Executor for the thread pool
    private final Optional<ThreadPoolExecutor> executor;
    // this is a simple way of encapsulating both multi-thread and no-multi thread processing
    private final Consumer<LDcomputation> runner;

    // queue with the light genotypes
    private final Deque<LightGenotype> snpQueue;

    // current contig
    private String currentContig = null;
    // bin Writer
    private final LDdecayOutput output;

    // PARAMETERS
    private final LengthBinningArgumentCollection binningParams;
    private final int minimumSamples;
    private final boolean rmSingletons;
    private final double chiSqrQuantile;

    /**
     * Initialize the queue for LD computation.
     *
     * @param output            the output for bin the results.
     * @param binningParams     the parameters for binning.
     * @param minimumSamples    the minimum number of samples to allow computation.
     * @param rmSingletons      if {@code true}, filter-out singletons.
     * @param chiSqrQuantile    the chi-square quantile to use for compute if max. correlation is
     *                          significant.
     * @param multiThreadParams parameters for multi-thread computation.
     */
    public QueueLD(final LDdecayOutput output, final LengthBinningArgumentCollection binningParams,
            final int minimumSamples, final boolean rmSingletons, final double chiSqrQuantile,
            final MultiThreadComputationArgumentCollection multiThreadParams) {
        this.binningParams = binningParams;
        this.minimumSamples = minimumSamples;
        this.chiSqrQuantile = chiSqrQuantile;
        this.rmSingletons = rmSingletons;
        this.executor = multiThreadParams.getFixedThreadPoolWithBoundedQueue();
        if (executor.isPresent()) {
            this.runner = executor.get()::execute;
        } else {
            logger.debug("Using QueueLD in the current thread.");
            this.runner = LDcomputation::run;
        }
        // linked list could handle up to 1.5 million variants (stored as LightGenotype)
        this.snpQueue = new LinkedList<>();
        this.computedPairs = new AtomicInteger(0);
        this.addedPairs = new AtomicInteger(0);
        this.output = output;

    }

    /** Returns the number of records in memory. */
    public int variantsInRam() {
        return snpQueue.size();
    }

    /** Gets the number of computed pairs. Note: this are the ones which passed all the filters. */
    public int computedPairs() {
        return computedPairs.intValue();
    }

    /** Gets the number of added pairs. Note: the unfiltered ones. */
    public int addedPairs() {
        return addedPairs.intValue();
    }

    /** Adds the variant to the queue and computes the queue if necessary. */
    public boolean add(final VariantContext variant) {
        // get the contig
        final String contig = variant.getContig();
        if (currentContig == null) {
            currentContig = contig;
        }
        // check if the queue is empty
        if (snpQueue.isEmpty() ||
                // or if it is in the contig and in the first variant, add directly
                (currentContig.equals(contig) && binningParams.bellowMaximumDistance(
                        variant.getStart(), snpQueue.peek().getPosition()))) {
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
            finalizeQueue();
            // we changed the contig
            currentContig = contig;
            // add variant to the queue
            return snpQueue.add(new LightGenotype(variant));
        }
    }

    /** Computes LD for the first element in the queue against all the other. */
    private void computeQueueLD() {
        // extract the first variant
        final LightGenotype firstVariant = snpQueue.remove();
        for (final LightGenotype second : snpQueue) {
            final LDcomputation toCompute = new LDcomputation(firstVariant, second);
            // run this job, either in this thread or in the thread pool
            runner.accept(toCompute);
        }
    }

    /** Writes the bins and reset them. */
    private void writeBins() {
        // wait till completion
        waitUntilFinish();
        output.write(currentContig);
        output.clear();
    }

    /** Waits for the executor to finish all the loaded task. */
    private void waitUntilFinish() {
        if (executor.isPresent()) {
            while (!executor.get().getQueue().isEmpty()) {
                try {
                    monitorLogging();
                    Thread.sleep(1000);
                } catch (InterruptedException e) {
                    Thread.currentThread().interrupt();
                }
            }
            monitorLogging();
        }
    }

    /** Logs the executor status with a monitor. */
    private synchronized void monitorLogging() {
        if (executor.isPresent()) {
            logger.debug(() -> String.format(
                    "[thread-pool] [%s/%s] Active: %s/%s | Completed: %s/%s | Shutdown/Terminated: %s/%s",
                    executor.get().getPoolSize(), executor.get().getCorePoolSize(),
                    executor.get().getActiveCount(), executor.get().getQueue().size(),
                    executor.get().getCompletedTaskCount(), executor.get().getTaskCount(),
                    executor.get().isShutdown(), executor.get().isTerminated()));
        }
    }

    /** Finalizes all the pending variants in the queue. */
    public void finalizeQueue() {
        logger.debug("Finalizing queue");
        // empty the queue
        while (!snpQueue.isEmpty()) {
            computeQueueLD();
        }
        // write the bins
        writeBins();
    }

    /** Closes the queue and finishes computation. */
    @Override
    public void close() {
        logger.debug("Closing queue");
        // finalize the current state
        finalizeQueue();
        if (executor.isPresent()) {
            logger.debug("Shutdown the executor");
            executor.get().shutdown();
        }
        output.close();
        monitorLogging();
    }

    /** Runnable class for compute SNP pairs and LD statistics. */
    private class LDcomputation implements Runnable {

        private final LightGenotype geno1, geno2;

        public LDcomputation(final LightGenotype genotype1, final LightGenotype genotype2) {
            this.geno1 = genotype1;
            this.geno2 = genotype2;
        }

        @Override
        public void run() {
            // only generate the pair if the minimum distance threshold is hold; may be faster than generate it
            if (binningParams.exceedMinimumDistance(geno1.getPosition(), geno2.getPosition())) {
                // generate the SNP pair
                final SNPpair pair = new SNPpair(geno1, geno2);
                addedPairs.incrementAndGet();
                // check if the conditions are met
                if (!(pair.getTotalCounts() < minimumSamples
                        // total counts are lower than the minimum samples
                        || pair.isInvariantA()                        // A is invariant
                        || pair.isInvariantB()                        // B is invariant
                        || (rmSingletons && pair.oneIsSingleton())))  // remove singletons
                {
                    // compute LD statistic for the pair
                    final double[] LD = LDfunctions.rStatistics(pair, chiSqrQuantile);
                    // add to the binning if computed
                    if (LD != null) {
                        output.add(pair, LD);
                        computedPairs.incrementAndGet();
                    }
                }
            }
        }
    }

}

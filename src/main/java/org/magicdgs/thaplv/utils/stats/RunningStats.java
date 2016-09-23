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

package org.magicdgs.thaplv.utils.stats;

import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.math3.stat.descriptive.rank.PSquarePercentile;

import java.util.Map;
import java.util.TreeMap;

/**
 * Class for computing running statistics including mean, standard deviation and quantiles (using
 * {@link PSquarePercentile}).
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
// TODO: synchronization could be done in a better/safer way
// TODO: change the implementation to use the exact implementations if the number of values added is low
public class RunningStats {
    // TODO: change the mean/stdev implementation for apache.commons.math3 methods?
    // TODO: will be much better for consistency
    private double newMean, newStdDev;
    private final PSquarePercentile median;
    private final Variance sampleVar;
    // the quantiles always includes the median
    private final Map<Double, PSquarePercentile> quantiles;

    /**
     * Constructs a new running statistic without quantiles. Computes the mean, variance, standard
     * deviation and median.
     */
    public RunningStats() {
        median = new PSquarePercentile(50.0D);
        sampleVar = new Variance();
        quantiles = new TreeMap<>();
    }

    /**
     * Constructs a new running statistics with the provided quantiles.
     *
     * @param quantiles quantiles to add in the range (1, 100).
     */
    public RunningStats(final double... quantiles) {
        this();
        for (final double quan : quantiles) {
            if (quan != 50.0D) {
                this.quantiles.put(quan, new PSquarePercentile(quan));
            }
        }
    }

    /**
     * Adds a new value to the statistics. {@link Double#NaN} will be ignored.
     *
     * @param x the value to add to the running statistics.
     *
     * @return {@code true} if the value is pushed; {@code true} otherwise
     */
    // TODO: change the name to add()
    public synchronized boolean push(double x) {
        if (!Double.isNaN(x)) {
            median.increment(x);
            // See Knuth TAOCP vol 2, 3rd edition, page 232
            if (median.getN() == 1) {
                newMean = x;
                newStdDev = 0.0;
            } else {
                final double oldMean = newMean;
                newMean = oldMean + (x - oldMean) / median.getN();
                newStdDev = newStdDev + (x - oldMean) * (x - newMean);
            }
            sampleVar.increment(x);
            for (final PSquarePercentile quan : quantiles.values()) {
                quan.increment(x);
            }
            return true;
        }
        return false;
    }

    /** Clear all the data in this instance. */
    public final synchronized void clear() {
        median.clear();
        for (final PSquarePercentile q : quantiles.values()) {
            q.clear();
        }
        newMean = 0;
        newStdDev = 0;
    }

    /** Get the number of data values. */
    public final synchronized long numDataValues() {
        return median.getN();
    }

    /**
     * Gets the actual mean.
     *
     * @return the current mean; {@link Double#NaN} if no value was added.
     */
    public final synchronized double mean() {
        return (median.getN() > 0) ? newMean : Double.NaN;
    }

    /**
     * Gets the actual variance.
     *
     * @return the current variance; {@link Double#NaN} if no value was added.
     */
    public final synchronized double variance() {
        return ((median.getN() > 1) ? newStdDev / (median.getN() - 1) : Double.NaN);
    }

    /**
     * Gets the actual sample variancel.
     *
     * @return the current sample variance; {@link Double#NaN} if no value was added.
     *
     * @deprecated use {@link #variance()} instead.
     */
    @Deprecated
    public final synchronized double sampleVariance() {
        return sampleVar.getResult();
    }

    /**
     * Gets the actual median.
     *
     * @return the current median; {@link Double#NaN} if no value was added.
     */
    public final synchronized double median() {
        return median.getResult();
    }

    /**
     * Returns {@code true} if the quantile is computed; {@code false} otherwise.
     */
    public final synchronized boolean hasQuantile(final double quantile) {
        return quantiles.containsKey(quantile);
    }

    /**
     * Gets the actual quantile value.
     *
     * @return the actual quantile; {@link Double#NaN} if no value was added.
     *
     * @throws IllegalArgumentException if the quantile is not computed.
     */
    public final synchronized double getQuantile(final double quantile) {
        if (quantile == 50.0D) {
            return median.getResult();
        }
        if (hasQuantile(quantile)) {
            return quantiles.get(quantile).getResult();
        }
        throw new IllegalArgumentException(
                "Running stats must be initialized with the requested quantile");
    }

    /**
     * Gets the current standard deviation.
     *
     * @return the current standard deviation;{@link Double#NaN} if no value was added.
     */
    public final synchronized double standardDeviation() {
        return Math.sqrt(variance());
    }

    /**
     * Gets the current sample standard deviation.
     *
     * @return the current sample standard deviation; {@link Double#NaN} if no value was added.
     *
     * @deprecated use {@link #standardDeviation()} instead.
     */
    @Deprecated
    public final synchronized double sampleStandardDeviation() {
        return Math.sqrt(sampleVariance());
    }

    /**
     * Gets all the quantiles except the median.
     *
     * @return a sorted map with the quantiles and its value.
     */
    public final synchronized Map<Double, Double> getAllQuantiles() {
        final TreeMap<Double, Double> toReturn = new TreeMap<>();
        quantiles.entrySet().stream()
                .forEach(q -> toReturn.put(q.getKey(), q.getValue().getResult()));
        return toReturn;
    }

    // for testing
    @Override
    public String toString() {
        final StringBuilder builder = new StringBuilder();
        builder.append("Stats for ");
        builder.append(numDataValues());
        builder.append(" datapoints. ");
        builder.append("Mean: ");
        builder.append(newMean);
        builder.append("; ");
        builder.append("Variance: ");
        builder.append(variance());
        builder.append("; ");
        builder.append("Median: ");
        builder.append(median.getResult());
        builder.append("; ");
        builder.append("Quantiles computed: ");
        builder.append(quantiles.keySet());
        return builder.toString();
    }
}

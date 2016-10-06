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

import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.stream.DoubleStream;

/**
 * Class that bin a double based on the length, computing mean, median and quantiles.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class LengthBinning {

    // value that contains the bins
    private final SortedMap<Integer, RunningStats> binStat;

    // value that contains the bin size
    private final int binSize;

    // which quantiles should be computed?
    private final double[] quantiles;

    private boolean empty;

    /**
     * Constructs a binning by length using the default quantiles (1, 5, 95, 99).
     *
     * @param binSize size of the bin.
     */
    public LengthBinning(final int binSize) {
        this(binSize, 1, 5, 95, 99);
    }

    /**
     * Constructs a binning by lenght using the provided quantiles.
     *
     * @param binSize   size of the bin.
     * @param quantiles quantiles to add in the range (1, 100).
     */
    public LengthBinning(final int binSize, final double... quantiles) {
        Utils.validateArg(binSize >= 1, "Binning cannot be performed for bins lower than 1 bp");
        this.binStat = new TreeMap<>();
        this.binSize = binSize;
        this.quantiles = quantiles;
        empty = true;
    }

    /** Returns {@code true} if no value was tracked in this object; {@code false} otherwise. */
    public synchronized boolean isEmpty() {
        return empty;
    }

    /**
     * Add a value to the associated length bin.
     *
     * @param length the length associated with this value.
     * @param value  the value to add to the computation of the statistic.
     */
    public synchronized void add(final int length, final double value) {
        final int cBin = toBin(length, binSize);
        RunningStats stats = binStat.get(cBin);
        if (stats == null) {
            stats = new RunningStats(quantiles);
            binStat.put(cBin, stats);
        }
        if (stats.add(value)) {
            empty = false;
        }
    }

    /** Gets the bin for a concrete value. */
    @VisibleForTesting
    static int toBin(final int length, final int binSize) {
        int bin = (int) Math.floor((length - 1) / binSize);
        bin *= binSize;
        bin += binSize;
        return bin;
    }

    // for testing
    @Override
    public synchronized String toString() {
        final StringBuilder builder = new StringBuilder();
        for (final Map.Entry<Integer, RunningStats> bin : binStat.entrySet()) {
            builder.append("\tBin #");
            builder.append(bin.getKey());
            builder.append("\t-> ");
            builder.append(bin.getValue());
            builder.append('\n');
        }
        return builder.toString();
    }

    /**
     * Gets the binned statistics for this object. The sorted map will contain all bins until the
     * next one.
     *
     * Note that the object is copied and modifications of the statistics will not be reflected in
     * this instance. Every time that the method is called it generates a new map.
     *
     * @return a copy of the binned statistics (including 0-element ones); emty map if {@link
     * #isEmpty()} is {@code true}.
     */
    public synchronized SortedMap<Integer, RunningStats> getBinStats() {
        final SortedMap<Integer, RunningStats> toReturn = new TreeMap<>();
        if (!isEmpty()) {
            final int lastBin = binStat.lastKey();
            for (int i = binSize; i <= lastBin; i += binSize) {
                final RunningStats stats = binStat.getOrDefault(i, new RunningStats(quantiles));
                toReturn.put(i, stats);
            }
        }
        return toReturn;
    }

    /** Returns a fresh array with the quantiles computed for the statistics. */
    public double[] getQuantiles() {
        return DoubleStream.of(quantiles).toArray();
    }

    /**
     * Clear all the statistics computed.
     *
     * After calling this method, {@link #isEmpty()} will return {@code true}.
     */
    public synchronized void clear() {
        binStat.clear();
        empty = true;
    }
}

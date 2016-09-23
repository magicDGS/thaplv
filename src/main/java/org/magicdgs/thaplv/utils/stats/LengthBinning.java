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

import java.text.DecimalFormat;
import java.util.Map;
import java.util.concurrent.ConcurrentSkipListMap;

/**
 * Class that bin a double based on the length, computing mean, median and quantiles.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
// TODO: implement unit tests
public class LengthBinning {

    // value that contains the bins
    // TODO: implement a getter?
    private final ConcurrentSkipListMap<Integer, RunningStats> binStat;

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
        // TODO: change contract to compute defaults with 0.005, 0.025, 0.975, 0.995 to have the 99% and 95%
        this(binSize, 1, 5, 95, 99);
    }

    /**
     * Constructs a binning by lenght using the provided quantiles.
     *
     * @param binSize   size of the bin.
     * @param quantiles quantiles to add in the range (1, 100).
     */
    public LengthBinning(final int binSize, final double... quantiles) {
        if (binSize < 1) {
            throw new IllegalArgumentException(
                    "Binning cannot be performed for bins lower than 1 bp");
        }
        this.binStat = new ConcurrentSkipListMap<>();
        this.binSize = binSize;
        this.quantiles = quantiles;
        empty = true;
    }

    /** Returns {@code true} if no value was tracked in this object; {@code false} otherwise. */
    public boolean isEmpty() {
        return empty;
    }

    /**
     * Add a value to the associated length bin.
     *
     * @param length the length associated with this value.
     * @param value  the value to add to the computation of the statistic.
     */
    public void add(final int length, final double value) {
        final int cBin = toBin(length);
        RunningStats stats = binStat.get(cBin);
        if (stats == null) {
            stats = new RunningStats(quantiles);
            binStat.put(cBin, stats);
        }
        if (stats.push(value)) {
            empty = false;
        }
    }

    /** Returns the string with the representation of the header for a table format of this object. */
    public String formatHeaderBins() {
        final StringBuilder builder = new StringBuilder("bin\tn\tmean\tsd\tmedian");
        for (final double q : quantiles) {
            builder.append("\tq");
            builder.append(q / 100);
        }
        return builder.toString();
    }

    /** Format the bins stored in this object. */
    public String formatBins(final DecimalFormat digits) {
        // TODO: format bins with no data???
        final StringBuilder builder = new StringBuilder();
        binStat.entrySet().stream().filter(bin -> bin.getValue().numDataValues() != 0)
                .forEach(bin -> {
                    builder.append('\t');
                    builder.append(bin.getKey());
                    builder.append('\t');
                    builder.append(bin.getValue().numDataValues());
                    builder.append('\t');
                    builder.append(digits.format(bin.getValue().mean()));
                    builder.append('\t');
                    builder.append(digits.format(bin.getValue().sampleStandardDeviation()));
                    builder.append('\t');
                    builder.append(digits.format(bin.getValue().median()));
                    builder.append('\t');
                    for (final double quantile : bin.getValue().getAllQuantiles().values()) {
                        builder.append(digits.format(quantile));
                        builder.append('\t');
                    }
                    builder.append('\n');
                });
        return builder.toString();
    }

    /** Getd the bin for a concrete value. */
    private int toBin(final int length) {
        int bin = (int) Math.floor((length - 1) / binSize);
        bin *= binSize;
        bin += binSize;
        return bin;
    }

    // for testing
    @Override
    public String toString() {
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
     * Clear all the statistics computed.
     *
     * After calling this method, {@link #isEmpty()} will return {@code true}.
     */
    public void clear() {
        binStat.clear();
        empty = true;
    }
}

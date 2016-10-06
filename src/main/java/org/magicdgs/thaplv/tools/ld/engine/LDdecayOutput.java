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

import org.magicdgs.thaplv.haplotypes.light.SNPpair;
import org.magicdgs.thaplv.utils.stats.LengthBinning;
import org.magicdgs.thaplv.utils.stats.popgen.LDfunctions;

import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.BufferedWriter;
import java.io.Closeable;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Output for LDdecay. It contains the abstraction of adding/writing the binning statistics.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 * @deprecated the engine for ld will be re-implemented from scratch
 */
@Deprecated
public class LDdecayOutput implements Closeable {

    private static final char TAB_SEPARATOR = '\t';

    // format for the values
    private static final DecimalFormat ROUNDED_FORMAT = new DecimalFormat("#.#######");

    private static final List<Double> quantiles = new ArrayList<Double>() {{
        add(1D);
        add(5D);
        add(95D);
        add(99D);
    }};

    private static List<String> stringQuantiles = quantiles.stream().map(q -> "q" + q / 100)
            .collect(Collectors.toList());

    private static final String baseHeader = String.join(String.valueOf(TAB_SEPARATOR),
            new String[] {"chr", "bin", "n", "mean", "sd", "median"});


    // list with the bins to write (for some statistics we don't need binning)
    @VisibleForTesting
    public static final List<String> statsToBin = new ArrayList<String>() {{
        add("r2");
        add("r2norm");
        add("rw_pos");
        add("rw_neg");
    }};

    // array that contains the binning in the order: positive rw, negative rw, r2 and r2'
    private final LengthBinning[] binningStats;
    private final BufferedWriter[] writers;
    private final File[] outputFiles;

    /**
     * Creates a LDdecay output.
     *
     * @param outputPrefix prefix for binned histograms.
     * @param binLength    lenght for binning.
     *
     * @throws UserException.CouldNotCreateOutputFile if an IO occurs
     */
    public LDdecayOutput(final String outputPrefix, final int binLength) {
        final String header = baseHeader + TAB_SEPARATOR
                + String.join(String.valueOf(TAB_SEPARATOR), stringQuantiles);

        // init a new Array with the length og the stats to bin
        binningStats = new LengthBinning[statsToBin.size()];
        writers = new BufferedWriter[statsToBin.size()];
        outputFiles = new File[statsToBin.size()];

        for (int i = 0; i < statsToBin.size(); i++) {
            binningStats[i] =
                    new LengthBinning(binLength, quantiles.stream().mapToDouble(d -> d).toArray());
            outputFiles[i] = new File(String.format("%s.%s", outputPrefix, statsToBin.get(i)));
            try {
                writers[i] = new BufferedWriter(new FileWriter(outputFiles[i]));
                writers[i].write(header);
                writers[i].newLine();
            } catch (IOException e) {
                throw new UserException.CouldNotCreateOutputFile(outputFiles[i], e);
            }
        }
    }

    /** Writes all the histograms to their respective files. */
    public void write(final String contig) {
        for (int i = 0; i < statsToBin.size(); i++) {
            // get the length binning
            try {
                writers[i].write(formatBins(contig, binningStats[i]));
            } catch (IOException e) {
                throw new UserException.CouldNotCreateOutputFile(outputFiles[i], e);
            }
        }

    }

    /** Converts a {@link LengthBinning} into a histogram string format, appending the contig to it. */
    private static String formatBins(final String contig, final LengthBinning bins) {
        final StringBuilder builder = new StringBuilder();
        bins.getBinStats().entrySet().stream().filter(bin -> bin.getValue().numDataValues() != 0)
                .forEach(bin -> {
                    builder.append(contig).append(TAB_SEPARATOR);
                    builder.append(bin.getKey()).append(TAB_SEPARATOR);
                    builder.append(bin.getValue().numDataValues()).append(TAB_SEPARATOR);
                    builder.append(ROUNDED_FORMAT.format(bin.getValue().mean()))
                            .append(TAB_SEPARATOR);
                    builder.append(ROUNDED_FORMAT.format(bin.getValue().sampleStandardDeviation()))
                            .append(TAB_SEPARATOR);
                    builder.append(ROUNDED_FORMAT.format(bin.getValue().median()))
                            .append(TAB_SEPARATOR);
                    for (final double quantile : bin.getValue().getAllQuantiles().values()) {
                        builder.append(ROUNDED_FORMAT.format(quantile)).append(TAB_SEPARATOR);
                    }
                    builder.append('\n');
                });
        return builder.toString();
    }


    /** Adds to a distance bin a value for the statistic idx. */
    private synchronized void addToBin(int idx, int distance, double value) {
        binningStats[idx].add(distance, value);
    }


    /** Adds a SNP pair and associated values. */
    public void add(final SNPpair pair, final double[] values) {
        final String chromosomeName = pair.getReferenceA();
        if (!chromosomeName.equals(pair.getReferenceB())) {
            throw new IllegalArgumentException(
                    "Output for pairs in different chromosomes not implemented");
        }
        add(pair.getDistance(), values);
    }

    /** Adds a new value to bin. */
    private void add(int distance, double[] values) {
        if (values.length != LDfunctions.rStatisticsOrder.length) {
            throw new IllegalArgumentException("Only " + LDfunctions.rStatisticsOrder.length
                    + " statistics could be written in the output file");
        }
        // iterate over the values
        for (int i = 0; i < LDfunctions.rStatisticsOrder.length; i++) {
            // update the binnings if this value is binned
            final int idx = statsToBin.indexOf(LDfunctions.rStatisticsOrder[i]);
            // if it have an index, added
            if (idx != -1) {
                // Only this needs to be synchronized
                addToBin(idx, distance, values[i]);
            } else {
                // check if it is asked to separate positive and negative values
                final int posIndex = statsToBin
                        .indexOf(String.format("%s_pos", LDfunctions.rStatisticsOrder[i]));
                final int negIndex = statsToBin
                        .indexOf(String.format("%s_neg", LDfunctions.rStatisticsOrder[i]));
                // bin according to the sign of the value
                // TODO: value of 0 is actually bin to the negative, but it could be added to the positive too or both/none
                // TODO: discussion in https://github.com/magicDGS/thaplv/issues/39
                if (posIndex != -1 && values[i] > 0) {
                    addToBin(posIndex, distance, values[i]);
                } else if (negIndex != -1 && values[i] <= 0) {
                    addToBin(negIndex, distance, values[i]);
                }

            }
        }
    }

    /** Clears all the already computed binned statistics. */
    public void clear() {
        for (int i = 0; i < statsToBin.size(); i++) {
            binningStats[i].clear();
        }
    }

    /** Closes all the output binned stats. */
    @Override
    public void close() {
        for (int i = 0; i < statsToBin.size(); i++) {
            // get the length binning
            try {
                writers[i].close();
            } catch (IOException e) {
                throw new UserException.CouldNotCreateOutputFile(outputFiles[i], e);
            }
        }
    }
}

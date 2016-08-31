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
package org.magicdgs.thaplv.tools.ibd.engine;

import org.magicdgs.thaplv.haplotypes.pairs.DifferencesDistancePair;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.Closeable;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Creates a new processor for IBD
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class IBDOutput implements Closeable {

    private final static Logger logger = LogManager.getLogger(IBDOutput.class);

    private final PrintStream ibdTracksWriter;

    private final PrintStream pairwiseDiffWritter;

    private final IBDcollector collector;

    private final double minDiff;

    private final LinkedHashMap<String, SimpleInterval> pairIBDmap;

    /**
     * Construct an output for IBD regions
     *
     * @param outputPrefix       the prefix for the output files
     * @param collector          collector for compute IBD regions
     * @param minimumDifferences the minimum number of differences
     * @param onlyIBD            should only the IBD-tracks be output?
     *
     * @throws UserException.CouldNotCreateOutputFile if there is an IO error
     */
    public IBDOutput(final String outputPrefix, final IBDcollector collector,
            final double minimumDifferences, final boolean onlyIBD) {
        final File ibdFile = new File(String.format("%s.ibd", outputPrefix));
        try {
            this.ibdTracksWriter = new PrintStream(ibdFile);
        } catch (FileNotFoundException e) {
            throw new UserException.CouldNotCreateOutputFile(ibdFile, e);
        }
        final File pairwiseDiffFile = new File(String.format("%s.diff", outputPrefix));
        try {
            this.pairwiseDiffWritter = (onlyIBD) ? null : new PrintStream(pairwiseDiffFile);
        } catch (FileNotFoundException e) {
            throw new UserException.CouldNotCreateOutputFile(pairwiseDiffFile, e);
        }
        this.collector = collector;
        this.minDiff = minimumDifferences;
        printHeader();
        pairIBDmap = new LinkedHashMap<>();
    }

    /**
     * Add the variant to the output. It adds the variant to the collector and outputs as much as
     * possible.
     *
     * @param variant the variant
     */
    public void addToCollector(final VariantContext variant) {
        final List<PairwiseDifferencesWindow> done = collector.addVariant(variant);
        endWindows(done, true);
    }

    /**
     * Print the windows as pair-wise differences if requested, and add to the IBD tracks.
     *
     * @param done      windows that are already finished
     * @param dontCheck add also regions with 0 variants included
     */
    private void endWindows(final Collection<PairwiseDifferencesWindow> done,
            final boolean dontCheck) {
        for (final PairwiseDifferencesWindow window : done) {
            if (dontCheck || window.getVariantsInWindow() != 0) {
                for (final DifferencesDistancePair pair : window) {
                    // first print the pairwise difference
                    printPairwiseDiff(window, pair);
                    // and later the IBD track
                    addToIBDtracks(pair, window);

                }
                flushWriters();
            }
        }
    }

    /** Add a pair to the IBD tracks and update if possible. In addition, output the IBD track. */
    private void addToIBDtracks(final DifferencesDistancePair pair,
            final PairwiseDifferencesWindow window) {
        // retrieve the window
        final SimpleInterval win = pairIBDmap.get(pair.getPairNames());
        // see if the window is bigger than the threshold and have sites
        if (differencesPerSite(pair, window) < minDiff && pair.getNumberOfSites() != 0) {
            // if there is a window present for this pair
            if (win != null) {
                // if it is intersecting
                // TODO: using overlaps() and meregeWithContiguous() breaks compatibility
                // TODO: this was checked by converting to bed file and joining intervals
                // TODO: now a test is covering this, but keeping this TODO for other tools that use this class
                if (win.overlaps(window)) {
                    // extend the window and return
                    pairIBDmap.put(pair.getPairNames(), win.mergeWithContiguous(window));
                    return;
                } else {
                    // if not, print the IBD region
                    printIBDregion(pair.getSample1(), pair.getSample2(), win);
                }
            }
            // either if the window is not extended or is null for this sample, store the new window
            pairIBDmap.put(pair.getPairNames(), window.getInterval());
        } else if (win != null && !win.overlaps(window)) {
            // if the window is not overlapping with the current, save space printing the IBD region
            printIBDregion(pair.getSample1(), pair.getSample2(), win);
            pairIBDmap.remove(pair.getPairNames());
        }
    }

    /** Print the header for the writers. */
    private void printHeader() {
        // write the header for IBD tracks
        ibdTracksWriter.println("Sample1\tSample2\tRef\tStart\tEnd");
        if (pairwiseDiffWritter == null) {
            logger.warn("Pairwise-difference file (.diff) won't be written");
        } else {
            logger.debug("Pairwise-difference file (.diff) will be written");
            pairwiseDiffWritter.println(
                    "Sample1\tSample2\tRef\tStart\t\tEnd\tWin_length_length\tN_variants\tMissing\tN_diferences\tDiff_per_site");
        }
    }

    /** Print the IBD region for a pair of samples */
    private void printIBDregion(final String sample1, final String sample2,
            final SimpleInterval window) {
        ibdTracksWriter.println(String.format("%s\t%s\t%s\t%d\t%d",
                sample1, sample2, window.getContig(), window.getStart(), window.getEnd()));
    }

    /** Print the pair-wise differences for a pair. */
    private void printPairwiseDiff(final PairwiseDifferencesWindow window,
            final DifferencesDistancePair pair) {
        if (pairwiseDiffWritter != null) {
            pairwiseDiffWritter.println(String.format("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%f",
                    pair.getSample1(), pair.getSample2(), window.getContig(), window.getStart(),
                    window.getEnd(),
                    window.getAvailableSites(), window.getVariantsInWindow(),
                    pair.getNumberOfMissing(), pair.getNumberOfDifferences(),
                    differencesPerSite(pair, window)));
        }
    }

    /** Computes the differences per site. */
    private double differencesPerSite(final DifferencesDistancePair pair,
            final PairwiseDifferencesWindow window) {
        return pair.getNumberOfDifferences() / (double) (window.getAvailableSites() - pair
                .getNumberOfMissing());
    }

    /** Flush the writers. */
    private void flushWriters() {
        ibdTracksWriter.flush();
        if (pairwiseDiffWritter != null) {
            pairwiseDiffWritter.flush();
        }
    }

    public void close() {
        logger.debug("Ending pending windows");
        // TODO: should we break compatibility using dontCheck = true?
        endWindows(collector.getWindows(), false);
        logger.debug("Output the rest of IBD tracks");
        // output the rest of IBD tracks
        if (!pairIBDmap.isEmpty()) {
            for (final Map.Entry<String, SimpleInterval> entry : pairIBDmap.entrySet()) {
                if (entry.getValue() != null) {
                    final String[] sampleNames = entry.getKey().split(",");
                    printIBDregion(sampleNames[0].substring(1),
                            sampleNames[1].substring(0, sampleNames[1].length() - 1),
                            entry.getValue());
                }
            }
        }
        flushWriters();
        // close the writers
        logger.debug("Trying to close the .ibd file");
        CloserUtil.close(ibdTracksWriter);
        if (pairwiseDiffWritter != null) {
            logger.debug("Trying to close the .diff file");
            CloserUtil.close(pairwiseDiffWritter);
        }
        logger.debug("Trying to close the collectors counter");
        CloserUtil.close(collector.nCounter);
    }
}

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
import org.magicdgs.thaplv.io.FastaNsCounter;

import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * Pair-wise difference window for compute differences
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class PairwiseDifferencesWindow implements Iterable<DifferencesDistancePair>, Locatable {

    private final static Logger logger = LogManager.getLogger(IBDOutput.class);

    private final SimpleInterval interval;

    /**
     * The ordered pairs
     */
    protected final List<DifferencesDistancePair> orderedPairsBySampleNameInHeader;

    /**
     * The number of availables sites in the window
     */
    protected final int availableSites;

    /**
     * The number of variants added to the window
     */
    protected int variantsInWindow;

    /** Is this window the last one over the genome */
    private final boolean isLast;

    /**
     * Private constructor for initialize everything except the pairs
     *
     * @param interval the interval configuring the window
     * @param counter  counter for the Ns in the FASTA file
     */
    private PairwiseDifferencesWindow(final SimpleInterval interval, final FastaNsCounter counter) {
        Utils.nonNull(interval, "null interval");
        Utils.nonNull(counter, "null counter");
        if (IntervalUtils.intervalIsOnDictionaryContig(interval, counter.getDictionary())) {
            this.interval = interval;
            this.isLast = interval.getEnd() == counter.getChromosomeLength(interval.getContig());
        } else {
            this.interval = interval.expandWithinContig(0, counter.getDictionary());
            this.isLast = true;
        }
        // TODO: this is -1 to keep consistency with the first implementation
        // TODO: we should check at some point if this is real
        this.availableSites = interval.size() - counter.countNsRegion(interval) - 1;
        this.variantsInWindow = 0;
        orderedPairsBySampleNameInHeader = new ArrayList<>();
        logger.debug("Generated PairwiseDifferencesWindow at {} with {} available sites.",
                interval, availableSites);
    }

    /**
     * Constructor for a pair-wise comparison against only one sample
     *
     * @param interval         the interval configuring the window
     * @param comparisonSample the sample to use for comparison
     * @param sampleNames      the rest of the samples to compare with {@code comparisonSample}
     * @param counter          counter for the Ns in the FASTA file
     */
    public PairwiseDifferencesWindow(final SimpleInterval interval, final String comparisonSample,
            final List<String> sampleNames, final FastaNsCounter counter) {
        this(interval, counter);
        Utils.nonNull(comparisonSample, "null comparisonSample");
        Utils.nonEmpty(sampleNames, "empty sampleNames");
        initPairs(comparisonSample, sampleNames);
    }

    /**
     * Constructor for all pair-wise comparison of the samples
     *
     * @param interval    the interval configuring the window
     * @param sampleNames the rest of the samples to compare with {@code comparisonSample}
     * @param counter     counter for the Ns in the FASTA file
     */
    public PairwiseDifferencesWindow(final SimpleInterval interval, final List<String> sampleNames,
            final FastaNsCounter counter) {
        this(interval, counter);
        Utils.nonEmpty(sampleNames, "empty sampleNames");
        initPairs(sampleNames);
    }

    public PairwiseDifferencesWindow(final String contig, final int start, final int end,
            final List<String> sampleNames, final FastaNsCounter counter) {
        this(new SimpleInterval(contig, start, end), sampleNames, counter);
    }

    public PairwiseDifferencesWindow(final String contig, final int start, final int end,
            final String comparisonSample, final List<String> sampleNames,
            final FastaNsCounter counter) {
        this(new SimpleInterval(contig, start, end), comparisonSample, sampleNames, counter);
    }

    /**
     * Get the number of variants in this window
     *
     * @return the number of variants accumulated
     */
    public int getVariantsInWindow() {
        return variantsInWindow;
    }

    /**
     * Get the number of available sites
     *
     * @return the number of available sites
     */
    public int getAvailableSites() {
        return availableSites;
    }

    /** Is this the last window in the genome? */
    public boolean isLast() {
        return isLast;
    }

    /**
     * Initialize the pairs.
     */
    private void initPairs(final String sampleName, final List<String> otherSamples) {
        logger.debug("Comparing {} against other {} samples",
                sampleName, otherSamples.size());
        for (final String sample : otherSamples) {
            // only if they are different
            if (!sampleName.equals(sample)) {
                orderedPairsBySampleNameInHeader
                        .add(new DifferencesDistancePair(sampleName, sample));
            } else {
                logger.debug("Equal samples found: ignoring this pair");
            }
        }

        logger.debug("Number of pair-wise comparisons to perform: {}",
                orderedPairsBySampleNameInHeader.size());
    }

    /**
     * Initialize the pairs
     */
    private void initPairs(final List<String> sampleNames) {
        for (int i = 0; i < sampleNames.size() - 1; i++) {
            for (int j = i + 1; j < sampleNames.size(); j++) {
                orderedPairsBySampleNameInHeader
                        .add(new DifferencesDistancePair(sampleNames.get(i),
                                sampleNames.get(j)));
            }
        }
    }

    /**
     * Get the number of pairs in the window
     *
     * @return the number of pairs
     */
    public int getNumberOfPairs() {
        return orderedPairsBySampleNameInHeader.size();
    }

    /**
     * Add the variant to the window
     *
     * @param variant the variant to add
     */
    public void addVariant(final VariantContext variant) {
        if (interval.overlaps(variant)) {
            for (final DifferencesDistancePair currentPair : orderedPairsBySampleNameInHeader) {
                currentPair.add(variant);
            }
            variantsInWindow++;
        }
    }

    /**
     * Add a variant for a reference comparison
     *
     * @param variant the variant to add
     */
    public void addVariantReferenceComparison(final VariantContext variant) {
        if (interval.overlaps(variant)) {
            for (final DifferencesDistancePair currentPair : orderedPairsBySampleNameInHeader) {
                currentPair.addReference(variant.getGenotype(currentPair.getSample2()));
            }
            variantsInWindow++;
        }
    }

    @Override
    public Iterator<DifferencesDistancePair> iterator() {
        return orderedPairsBySampleNameInHeader.iterator();
    }

    @Override
    public String getContig() {
        return interval.getContig();
    }

    @Override
    public int getStart() {
        return interval.getStart();
    }

    @Override
    public int getEnd() {
        return interval.getEnd();
    }

    public SimpleInterval getInterval() {
        return interval;
    }
}

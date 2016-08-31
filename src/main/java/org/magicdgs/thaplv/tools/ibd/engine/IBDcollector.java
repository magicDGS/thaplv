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

import org.magicdgs.thaplv.io.FastaNsCounter;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.Deque;
import java.util.LinkedList;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.function.BiFunction;
import java.util.function.Consumer;
import java.util.function.Function;

/**
 * Collector for variants to compute IBD.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class IBDcollector {

    // logger for this class
    private final static Logger logger = LogManager.getLogger(IBDcollector.class);

    // queue for store difference windows until they are not overlapping
    private final Deque<PairwiseDifferencesWindow> queue;

    // the name of the samples
    private final List<String> sampleNames;

    // a simple class for count the number of Ns in a sequence
    protected final FastaNsCounter nCounter;

    // the window size
    private final int windowSize;

    // the window step
    private final int windowStep;

    // cached maximum number of windows in RAM
    private final int maximumWindow;

    private final BiConsumer<PairwiseDifferencesWindow, VariantContext> function;

    /** Creator for the windows using the contig and the start. */
    private final BiFunction<String, Integer, PairwiseDifferencesWindow> windowsCreator;

    /**
     * New IBD collector for the samples in a VCF header
     *
     * @param header     the header where the samples are stored
     * @param nCounter   the class for count the Ns in the reference
     * @param windowSize the window size
     * @param windowStep the window step
     */
    public IBDcollector(final VCFHeader header, final FastaNsCounter nCounter, final int windowSize,
            final int windowStep) {
        Utils.nonNull(header, "null header");
        Utils.nonNull(nCounter, "null nCounter");
        Utils.nonNull(nCounter.getDictionary(), "null nCounter dictionary");
        this.sampleNames = header.getSampleNamesInOrder();
        this.nCounter = nCounter;
        this.windowSize = windowSize;
        this.windowStep = windowStep;
        // compute the cached number of windows
        this.maximumWindow = computeMaximumNumberOfWindows(windowSize, windowStep);
        this.queue = new LinkedList<>();
        this.function = PairwiseDifferencesWindow::addVariant;
        this.windowsCreator = (contig, start) -> new PairwiseDifferencesWindow(contig, start,
                start + this.windowSize, sampleNames, this.nCounter);
    }

    /**
     * New IBD collector for the reference comparison
     *
     * @param sampleNames the name of the samples to include
     * @param nCounter   the class for count the Ns in the reference
     * @param windowSize the window size
     * @param windowStep the window step
     */
    public IBDcollector(final List<String> sampleNames, final FastaNsCounter nCounter, final int windowSize,  final int windowStep) {
        Utils.nonNull(nCounter, "null nCounter");
        Utils.nonNull(nCounter.getDictionary(), "null nCounter dictionary");
        this.sampleNames = sampleNames;
        this.nCounter = nCounter;
        this.windowSize = windowSize;
        this.windowStep = windowStep;
        // compute the cached number of windows
        this.maximumWindow = computeMaximumNumberOfWindows(windowSize, windowStep);
        this.queue = new LinkedList<>();
        this.function = PairwiseDifferencesWindow::addVariantReferenceComparison;
        this.windowsCreator = (contig, start) -> new PairwiseDifferencesWindow(contig, start,
                start + this.windowSize, "Reference", sampleNames, this.nCounter);
    }

    /** Compute the cached number of windows. */
    private static int computeMaximumNumberOfWindows(final int windowSize, final int windowStep) {
        return (int) (Math.floor(windowSize / windowStep)) - 1;
    }

    /**
     * Add a variant and return all the pairwise-differences window that are not longer needed
     *
     * @param variant the variant to add
     *
     * @return windows that do not overlap with the variant
     */
    public List<PairwiseDifferencesWindow> addVariant(final VariantContext variant) {
        Utils.nonNull(variant);
        final List<PairwiseDifferencesWindow> windows = new ArrayList<PairwiseDifferencesWindow>();
        // keep removing the first windows until the queue is empty or the variant is in the window
        while ((!queue.isEmpty()) && !queue.peek().getInterval().overlaps(variant)) {
            windows.add(queue.pop());
        }
        // if the queue is empty, generate the first queue from this variant and add
        if (queue.isEmpty()) {
            // logger.debug("Queue is empty");
            addFirstVariant(variant);
        } else {
            // compute necessary windows for the ws and ss configuration and without the already in the queue
            final int moreWin = maximumWindow - queue.size();
            // if there are still some windows to add
            if (moreWin > 0) {
                logger.debug("We need {} more windows", moreWin);
                // the next window start at the last window start + the step
                final int nextWinStart = queue.peekLast().getStart() + windowStep;
                // and the last window to add should end at
                final int end = nextWinStart + (windowStep * (moreWin));
                // add the rest of the windows
                initPairWiseDifferencesWindows(variant.getContig(), nextWinStart, end);
            }
            // add the variant
            addVariantToQueue(variant);
        }
        return windows;
    }

    /**
     * Get the current queue. Useful for processing the last variants.
     *
     * WARNING: it is not a copy, so it will be modified
     *
     * @return the current queue
     */
    Deque<PairwiseDifferencesWindow> getWindows() {
        return queue;
    }

    /**
     * Adds the first variant, checking what is the queue
     *
     * @param variant the variant to add
     */
    private void addFirstVariant(final VariantContext variant) {
        logger.debug("Generating first queue");
        // calculate the modulo for the step size
        final int mod = variant.getStart() % windowStep;
        switch (mod) {
            // case of the variant being the end of the window
            case 0:
                initPairWiseDifferencesWindow(variant.getContig(),
                        variant.getStart() - windowSize + 1);
                break;
            // case the variant start a window
            case 1:
                initPairWiseDifferencesWindow(variant.getContig(), variant.getStart());
                break;
            // case it is in the middle of a window
            default:
                initPairWiseDifferencesWindow(variant.getContig(), variant.getStart() - mod + 1);
        }
        addVariantToQueue(variant);
        logger.debug("Added first variant ({}:{}) is added, {} windows are kept in RAM",
                variant.getContig(), variant.getStart(), queue.size());
    }

    /**
     * Add the variant to the loaded queue
     */
    private void addVariantToQueue(final VariantContext variant) {
        queue.stream().forEach(win -> function.accept(win, variant));
    }

    /**
     * Create a single window, starting for that position
     *
     * @param contig the contig
     * @param start  the start
     */
    private void initPairWiseDifferencesWindow(final String contig, final int start) {
        initPairWiseDifferencesWindows(contig, start, start + windowSize);
    }

    /**
     * Create all the windows advancing from start to the end.
     *
     * Note: the end is not included
     *
     * @param contig the contig
     * @param start  the start
     * @param end    the end (non inclusive)
     */
    private void initPairWiseDifferencesWindows(final String contig, final int start,
            final int end) {
        for (int i = start; i < end; i += (windowStep)) {
            final PairwiseDifferencesWindow lastWindow = queue.peekLast();
            if (lastWindow != null && lastWindow.isLast()) {
                logger.debug("Could not create more windows after {}", lastWindow.getInterval());
                break;
            }
            queue.add(windowsCreator.apply(contig, i));
        }
    }
}

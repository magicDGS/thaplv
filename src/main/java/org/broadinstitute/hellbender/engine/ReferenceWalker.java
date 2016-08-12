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

package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.Collections;

/**
 * Traverse intervals in the reference.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public abstract class ReferenceWalker extends GATKTool {

    @Override
    public final boolean requiresReference() {
        return true;
    }

    /**
     * Customize initialization of the Feature data source for this traversal type to disable query
     * lookahead.
     */
    @Override
    void initializeFeatures() {
        // Disable query lookahead in our FeatureManager for this traversal type. Query lookahead helps
        // when our query intervals are overlapping and gradually increasing in position (as they are
        // with ReadWalkers, typically), but with ReferenceWalker our query intervals are guaranteed
        // to be non-overlapping, since our interval parsing code always merges overlapping intervals.
        features = new FeatureManager(this, 0);
        if (features.isEmpty()) {  // No available sources of Features for this tool
            features = null;
        }
    }

    /**
     * Initialize intervals with the whole genome if not provided
     */
    void initializeIntervals() {
        final SAMSequenceDictionary dictionary = getReferenceDictionary();
        if (dictionary == null) {
            throw new UserException(
                    "We currently require a sequence dictionary from the reference to process intervals for a "
                            + ReferenceWalker.class.getSimpleName()
                            + ". This restriction may be removed in the future.");
        }
        if (intervalArgumentCollection.intervalsSpecified()) {
            intervalsForTraversal = intervalArgumentCollection.getIntervals(dictionary);
        } else {
            intervalsForTraversal =
                    IntervalUtils.getAllIntervalsForReference(dictionary);
        }
    }

    /**
     * Initialize data sources.
     *
     * Marked final so that tool authors don't override it. Tool authors should override
     * onTraversalSuccess() instead.
     */
    @Override
    protected final void onStartup() {
        // Overridden only to make final so that concrete tool implementations don't override
        super.onStartup();
    }

    private SimpleInterval currentInterval;

    @Override
    public void traverse() {
        for (final SimpleInterval interval : intervalsForTraversal) {
            currentInterval = interval;
            IntervalUtils.cutToShards(Collections.singletonList(interval), 1).stream()
                    .forEach(position -> {
                        apply(position,
                                new ReadsContext(reads, position),
                                new ReferenceContext(reference, position),
                                new FeatureContext(features, position));
                        progressMeter.update(position);
                    });

        }
    }

    public SimpleInterval getCurrentInterval() {
        return currentInterval;
    }

    public abstract void apply(SimpleInterval position, ReadsContext gatkReads,
            ReferenceContext referenceContext, FeatureContext featureContext);

}

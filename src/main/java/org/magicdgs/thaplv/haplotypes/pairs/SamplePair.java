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

package org.magicdgs.thaplv.haplotypes.pairs;

import org.broadinstitute.hellbender.utils.Utils;

/**
 * Simple implementation of a sample pair with common methods for checking if samples belongs to
 * the pair.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class SamplePair {

    // the name of the first sample
    protected final String sample1;

    // the name of the second sample
    protected final String sample2;

    /**
     * Initialize the pair of samples with their names
     *
     * @param sample1 name for the first sample
     * @param sample2 name for the second sample
     */
    public SamplePair(final String sample1, final String sample2) {
        Utils.nonNull(sample1, "null sample1");
        Utils.nonNull(sample2, "null sample2");
        Utils.validateArg(!sample1.equals(sample2), "equal samples");
        this.sample1 = sample1;
        this.sample2 = sample2;
    }

    /**
     * Get the name for the first sample
     *
     * @return the name for the first sample
     */
    public String getSample1() {
        return sample1;
    }

    /**
     * Get the name for the second sample
     *
     * @return the name for the second sample
     */
    public String getSample2() {
        return sample2;
    }

    /**
     * Get the name for the pair
     *
     * @return the formatted name for the pairs
     */
    public String getPairNames() {
        return String.format("(%s,%s)", sample1, sample2);
    }

    /**
     * Check if the pair contain a concrete sample
     *
     * @param sample the sample to test
     *
     * @return {@code true} if it contains the sample; {@code false} otherwise
     */
    public boolean containSample(final String sample) {
        return sample1.equals(sample) || sample2.equals(sample);
    }

    /**
     * Check if the pair have these two samples
     *
     * @param name1   the name for the first sample
     * @param name2   the name for the second sample
     * @param ordered if {@code true}, name1 should be the first sample and name2 the second
     *
     * @return {@code true} if the pair represent these samples (ordered or not, depending on {@code
     * ordered}); {@code false} otherwise
     */
    public boolean containNames(final String name1, final String name2, boolean ordered) {
        if (ordered) {
            return sample1.equals(name1) && sample2.equals(name2);
        }
        return containSample(name1) && containSample(name2);
    }

    /**
     * Check if the pair represent these two samples (independently of the order)
     *
     * @param name1 the name for the first sample
     * @param name2 the the name for the second sample
     *
     * @return {@code true} if the pair represent these sample names; {@code false} otherwise
     */
    public boolean containNames(final String name1, final String name2) {
        return containNames(name1, name2, false);
    }
}

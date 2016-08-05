/*
 * Copyright (c) 2016, Daniel Gomez-Sanchez <daniel.gomez.sanchez@hotmail All rights reserved.
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

package org.magicdgs.thaplv.haplotypes.filters;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.filters.VariantFilter;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * Filter based in number of samples with missing genotypes (containing only {@link
 * Allele#NO_CALL}).
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class NumberOfMissingFilter implements VariantFilter {

    // maximum number of missing individuals
    private final int maxMissing;

    /**
     * Constructor
     *
     * @param maxMissing maximum number of samples missing (inclusive)
     */
    public NumberOfMissingFilter(final int maxMissing) {
        Utils.validateArg(maxMissing >= 0, "maxMisisng < 0");
        this.maxMissing = maxMissing;
    }

    /**
     * Returns {@code false} if there is more than maxMissing no-called genotypes; {@code true}
     * otherwise.
     */
    @Override
    public boolean test(final VariantContext variant) {
        int missing = 0;
        for (final Genotype geno : variant.getGenotypes()) {
            if (geno.isNoCall() && ++missing > maxMissing) {
                // early termination
                return false;
            }
        }
        return true;
    }

    /**
     * Construct a filter from a percentage of missing samples.
     *
     * @param maxPercentage maximum percentage of missing samples [0, 100). When 0, no missing data
     *                      is allowed
     * @param nSamples      the number of samples
     *
     * @return a filter for the number of samples
     */
    public static NumberOfMissingFilter fromMaxPercentage(final double maxPercentage,
            final int nSamples) {
        Utils.validateArg(maxPercentage >= 0 && maxPercentage < 100,
                "percentage is not in range [0, 100)");
        final int maxMissing = (maxPercentage == 0) ? 0 : (int) (maxPercentage / 100 * nSamples);
        return new NumberOfMissingFilter(maxMissing);
    }

    @Override
    public String toString() {
        return String.format("%s for %s max. number of missing genotypes",
                this.getClass().getSimpleName(), maxMissing);
    }

    /** Implemented to assess equal filters using percentage or number or maximum samples */
    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (!(o instanceof NumberOfMissingFilter)) {
            return false;
        }

        NumberOfMissingFilter that = (NumberOfMissingFilter) o;

        return maxMissing == that.maxMissing;

    }

    @Override
    public int hashCode() {
        return maxMissing;
    }
}

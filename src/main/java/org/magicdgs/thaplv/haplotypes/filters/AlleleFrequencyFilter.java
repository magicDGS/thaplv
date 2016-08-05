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

import org.magicdgs.thaplv.utils.AlleleUtils;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.filters.VariantFilter;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.HashMap;

/**
 * Filter based in allele frequency found in the genotypes (the ones with {@link Allele#NO_CALL}
 * will be ignored).
 *
 * WARNING: This filter only works for haploid/homozygous individuals, not taking into account
 * the second allele
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class AlleleFrequencyFilter implements VariantFilter {

    // minimum allele frequency
    private final double minFreq;

    // maximum allele frequency
    private final double maxFreq;

    /**
     * Constructor
     *
     * @param minFreq the minimum frequency in the samples to pass the filter (included). Set to
     *                0 to remove
     * @param maxFreq the maximum frequency in the samples to pass the filter (included). Set to
     *                1 to remove
     */
    public AlleleFrequencyFilter(final double minFreq, final double maxFreq) {
        Utils.validateArg(minFreq >= 0d && minFreq <= 1d, "minFreq not in range [0, 1]");
        Utils.validateArg(maxFreq >= 0d && maxFreq <= 1d, "maxFreq not in range [0, 1]");
        Utils.validateArg(minFreq < maxFreq, "minFreq >= maxFreq");
        this.minFreq = minFreq;
        this.maxFreq = maxFreq;
    }

    /**
     * The allele frequency for all the alleles present in the haplotypes (first allele in the
     * genotypes) is computed.
     *
     * Returns {@code true} if the range of frequencies is within [minFreq, maxFreq], {@code false}
     * otherwise
     */
    @Override
    public boolean test(final VariantContext variant) {
        final HashMap<Allele, Double> alleleCounts = AlleleUtils.getAlleleFrequencies(
                variant.getGenotypes().stream().map(g -> g.getAllele(0)), false);
        return alleleCounts.values().stream().noneMatch(f -> f < minFreq || f > maxFreq);
    }

    @Override
    public String toString() {
        return String.format("%s for allele frequencies in range [%s, %s]",
                this.getClass().getSimpleName(), minFreq, maxFreq);
    }
}

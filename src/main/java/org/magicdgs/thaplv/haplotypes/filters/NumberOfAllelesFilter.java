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

package org.magicdgs.thaplv.haplotypes.filters;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.filters.VariantFilter;

/**
 * Filter based on the number of distinct alleles found in the genotypes (except {@link
 * Allele#NO_CALL}.
 *
 * WARNING: This filter only works for haploid/homozygous individuals, not taking into account
 * the second allele
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class NumberOfAllelesFilter implements VariantFilter {

    // minimum number of distinct alleles
    private final int minAlleles;

    // maximum number of distinct alleles
    private final int maxAlleles;

    /**
     * Constructor
     *
     * @param minAlleles the minimum number of alleles to pass the filter (included). Set to -1 to
     *                   remove
     * @param maxAlleles the maximum number of alleles to pass the filter (included). Set to
     *                   {@link Integer#MAX_VALUE} to remove
     */
    public NumberOfAllelesFilter(final int minAlleles, final int maxAlleles) {
        this.minAlleles = minAlleles;
        this.maxAlleles = maxAlleles;
    }

    /**
     * Returns {@code true} if the number of distinct alleles in the haplotypes (first allele in
     * genotypes) is within the range [minAlleles, maxAlleles]; {@code false} otherwise.
     */
    @Override
    public boolean test(final VariantContext variant) {
        final int distinctAlleles = (int) variant.getGenotypes().stream()
                .map(g -> g.getAllele(0)) // map to the first alleles
                .filter(Allele::isCalled) // filter the no called
                .distinct().count(); // count the distinct ones
        return distinctAlleles >= minAlleles && distinctAlleles <= maxAlleles;
    }

    @Override
    public String toString() {
        return String
                .format("%s for number of different alleles in range [%s, %s]",
                        this.getClass().getSimpleName(), minAlleles, maxAlleles);
    }

}

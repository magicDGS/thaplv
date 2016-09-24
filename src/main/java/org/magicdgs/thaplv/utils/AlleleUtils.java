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

package org.magicdgs.thaplv.utils;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.HashMap;
import java.util.Map;
import java.util.stream.Stream;

/**
 * Static utils for {@link htsjdk.variant.variantcontext.Genotype} and {@link VariantContext} for
 * haploid individuals.
 *
 * WARNING: only works for haploid/homozygous individuals, not taking into account the second
 * allele when parsing genotypes
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class AlleleUtils {

    private AlleleUtils() { }

    /**
     * Get the difference value between two alleles (scaled by a factor of 2)
     *
     * @param allele1 the first allele
     * @param allele2 the second allele
     *
     * @return 1 if they are missing, 0 if they are equals, 2 otherwise
     */
    public static int scaledDifference(final Allele allele1, final Allele allele2) {
        if (allele1.isNoCall() || allele2.isNoCall()) {
            return 1;
        } else if (allele1.equals(allele2)) {
            return 0;
        } else {
            return 2;
        }
    }

    /**
     * Get the difference value between two alleles
     *
     * @param allele1 the first allele
     * @param allele2 the second allele
     *
     * @return 0 if they are equals, 0.5 if one of them is missing, 1 otherwise
     */
    public static double difference(Allele allele1, Allele allele2) {
        final int diff = scaledDifference(allele1, allele2);
        return (diff == 0) ? 0 : diff / 2d;
    }

    /**
     * Get the allele counts from an allele stream
     *
     * @param alleles       alleles to count
     * @param includeNoCall include in the counts the no call allele
     *
     * @return map of alleles to counts
     */
    public static final Map<Allele, Integer> getAlleleCounts(final Stream<Allele> alleles,
            final boolean includeNoCall) {
        final HashMap<Allele, Integer> alleleCounts = new HashMap<>();
        alleles.filter(a -> includeNoCall || a.isCalled()) // filter the no called
                .forEach(al -> alleleCounts.merge(al, 1, Integer::sum));
        return alleleCounts;
    }

    /**
     * Get the allele frequencies from an allele stream
     *
     * @param alleles       alleles to get the allele frequencies from
     * @param includeNoCall include in the counts the no call allele
     *
     * @return map of alleles to counts
     */
    public static final Map<Allele, Double> getAlleleFrequencies(final Stream<Allele> alleles,
            final boolean includeNoCall) {
        final HashMap<Allele, Double> alleleFreqs = new HashMap<>();
        final int total = (int) alleles
                .filter(a -> includeNoCall || a.isCalled()) // filter the no called
                .map(al -> alleleFreqs.merge(al, 1d, Double::sum)) // update the map
                .count(); // count how many alleles we used
        alleleFreqs.entrySet().stream()
                .forEach(e -> alleleFreqs.put(e.getKey(), e.getValue() / total));
        return alleleFreqs;
    }

}

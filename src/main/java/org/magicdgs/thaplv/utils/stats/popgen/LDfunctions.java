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

package org.magicdgs.thaplv.utils.stats.popgen;

import org.magicdgs.thaplv.haplotypes.light.SNPpair;
import org.magicdgs.thaplv.utils.stats.StatUtils;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;

/**
 * Functions for compute linkage disequilibrium (LD) related statistics. The methods works for
 * already polarized data, which means that the frequency of the major allele is the one provided.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
// TODO: almost all this functions could be memoized
// TODO: and the data structure for input could be a container with frequencies of AB, A and B
public class LDfunctions {

    /** Cannot be instantiated. */
    private LDfunctions() { }

    /** Chi-square distribution with 1 degree of freedom. */
    private static final ChiSquaredDistribution chiSqr = new ChiSquaredDistribution(null, 1);

    /**
     * Returned LD statistics with rStatistics are in the following order.
     *
     * @deprecated return a map instead
     */
    @Deprecated
    public static final String[] rStatisticsOrder = new String[] {"r2max", "r2", "r2norm", "rw"};

    /**
     * Computes signed correlation (rw). It needs already polarized frequencies.
     *
     * @param pAB allele frequency of haplotypes with major alleles in locus A and B.
     * @param pA  major allele frequency of locus A.
     * @param pB  the allele frequency of locus B.
     *
     * @return signed correlation.
     */
    public static double rw(final double pAB, final double pA, final double pB) {
        StatUtils.validateFrequency(pAB);
        StatUtils.validateFrequency(pA);
        StatUtils.validateFrequency(pB);
        return (pAB - pA * pB) / Math.sqrt(pA * (1 - pA) * pB * (1 - pB));
    }

    /**
     * Computes rw for a SNPpair.
     *
     * @deprecated use {@link #rw(double, double, double)} instead
     */
    @Deprecated
    public static double rw(final SNPpair pair) {
        return rw(pair.getFreqForHAB(), pair.getMajorAlleleFrequencyForA(),
                pair.getMajorAlleleFrequencyForB());
    }

    /**
     * Computes r2 based on the signed rw.
     *
     * @param rw the signed correlation.
     *
     * @return pearson correlation.
     */
    private static double r2(final double rw) {
        return Math.pow(Math.abs(rw), 2);
    }

    /**
     * Computes pearson correlation (r2) from allele frequencies.
     *
     * @param pAB allele frequency of haplotypes with major alleles in locus A and B.
     * @param pA  major allele frequency of locus A.
     * @param pB  the allele frequency of locus B.
     *
     * @return pearson correlation.
     */
    public static double r2(final double pAB, final double pA, final double pB) {
        return r2(rw(pAB, pA, pB));
    }

    /**
     * Computes r2 for a SNPpair.
     *
     * @deprecated use {@link #r2(double, double, double)} directly.
     */
    @Deprecated
    public static double r2(final SNPpair pair) {
        return r2(pair.getFreqForHAB(), pair.getMajorAlleleFrequencyForA(),
                pair.getMajorAlleleFrequencyForB());
    }

    /**
     * Computes maximum correlation with the frequencies. Facilitates the computation of maximum r2
     * when both are already polarized.
     *
     * @param smaller the smaller frequency.
     * @param higher  the higher frequency.
     *
     * @return the maximum correlation.
     */
    private static double maxR2polarized(final double smaller, final double higher) {
        return (smaller * (1 - higher)) / ((1 - smaller) * higher);
    }

    /**
     * Computes maximum correlation with the frequencies.
     *
     * @param pA major allele frequency of locus A.
     * @param pB the allele frequency of locus B.
     *
     * @return the maximum correlation.
     */
    public static double maxR2(final double pA, final double pB) {
        StatUtils.validateFrequency(pA);
        StatUtils.validateFrequency(pB);
        if (pA == pB) {
            return 1;
        }
        return (pA < pB) ? maxR2polarized(pA, pB) : maxR2polarized(pB, pA);
    }

    // TODO: javadoc
    public static double r2norm(final double pAB, final double pA, final double pB) {
        final double maxR2 = maxR2(pA, pB);
        final double r2 = r2(pAB, pA, pB);
        return r2 / maxR2;
    }

    /**
     * Computes maximum correlation for a SNPpair.
     *
     * @deprecated use {@link #maxR2(double, double)} directly
     */
    @Deprecated
    public static double maxR2(final SNPpair pair) {
        return maxR2(pair.getMajorAlleleFrequencyForA(), pair.getMajorAlleleFrequencyForB());
    }

    /**
     * Computes all the r statistics for a SNPpair with a threshold for minimum counts and
     * chi-squared quantile for the significance of maxR2.
     *
     * @param pair           polarized pair.
     * @param chiSqrQuantile the chi-squared quantile to asses significance of maxR2.
     *
     * @return array with the values if rmax is upper the threshold; null otherwise
     *
     * @deprecated use the single statistics directly.
     */
    @Deprecated
    public static double[] rStatistics(final SNPpair pair, final double chiSqrQuantile) {
        // maf for A and B
        final double pA = pair.getMajorAlleleFrequencyForA();
        final double pB = pair.getMajorAlleleFrequencyForB();
        // compute the maximum r2 that can be reached
        final double maxR2 = maxR2(pA, pB);
        // if the maximum r2 is not significant, it is not worthy to compute the rest of the stuff
        if (!isSignificantR2(maxR2, pair.getTotalCounts(), chiSqrQuantile)) {
            return null;
        }
        // compute rw and r2 (fast way)
        final double pAB = pair.getFreqForHAB();
        final double rw = rw(pAB, pA, pB);
        if (rw == 0) {
            return new double[] {maxR2, 0, 0, 0};
        }
        final double r2 = r2(rw);
        return new double[] {maxR2, r2, r2 / maxR2, rw};
    }

    /**
     * Tests if a value for pearson correlation (r2) is significant based on the chi-square
     * distribution for the provided number of haplotypes.
     *
     * Note: if performed to a maximum pearson correlation it will test for the significant level
     * that could be reached. This may be useful for filter pairs.
     *
     * @param r2                 the r2 value to test.
     * @param numberOfHaplotypes the number of haplotypes.
     * @param chiSqrQuantile     the chi-squared quantile to asses significance of the test.
     *
     * @return {@code true} if r2 is significant; {@code false} otherwise.
     */
    public static boolean isSignificantR2(final double r2, final int numberOfHaplotypes,
            final double chiSqrQuantile) {
        return r2 >= significantThreshold(numberOfHaplotypes, chiSqrQuantile);
    }

    /**
     * Return the minimum pearson correlation that could be significant for the provided number of
     * haplotypes.
     *
     * @param numberOfHaplotypes the number of haplotypes.
     * @param chiSqrQuantile     the chi-squared quantile to asses significance of the test.
     *
     * @return the minimum r2 that can be considered significantly different from 0
     */
    public static double significantThreshold(final int numberOfHaplotypes,
            final double chiSqrQuantile) {
        // TODO: this could be memoized to improve performance because it is expected to
        // TODO: be called for a range from 0 to some number of haplotypes and the same quantile.
        // TODO: nevertheless it requires a double memoizer...
        return chiSqr.inverseCumulativeProbability(chiSqrQuantile) / numberOfHaplotypes;
    }

}

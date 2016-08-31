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

package org.magicdgs.thaplv.utils.ld;

import org.magicdgs.thaplv.haplotypes.light.SNPpair;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;

/**
 * Function for compute LD related statistics
 *
 * @author Daniel Gómez-Sánchez
 */
public class LDfunctions {

    /**
     * Chi-square distribution with 1 degree of fredom
     */
    private static final ChiSquaredDistribution chiSqr = new ChiSquaredDistribution(null, 1);

    /**
     * Returned LD statistics with rStatistics are in the following order
     */
    public static final String[] rStatisticsOrder = new String[] {"r2max", "r2", "r2norm", "rw"};

    /**
     * Compute rw (signed correlation). It needs already polarized frequencies
     *
     * @param pAB the frequency of AB
     * @param pA  the frequency of A
     * @param pB  the frequency of B
     *
     * @return rw
     */
    private static double rw(double pAB, double pA, double pB) {
        return (pAB - pA * pB) / Math.sqrt(pA * (1 - pA) * pB * (1 - pB));
    }

    /**
     * Compute rw for a SNPpair (it should be polarized)
     *
     * @param pair polarized pair
     *
     * @return rw
     */
    public static double rw(SNPpair pair) {
        double pA = pair.getMajorAF_A();
        double pB = pair.getMajorAF_B();
        double pAB = pair.getFreq_AB();
        return rw(pAB, pA, pB);
    }

    /**
     * Compute r2 based on the signed rw
     *
     * @param rw the signed correlation
     *
     * @return r2
     */
    private static double r2(double rw) {
        return Math.pow(Math.abs(rw), 2);
    }

    /**
     * Compute r2 from allele frequencies
     *
     * @param pAB the frequency of AB
     * @param pA  the frequency of A
     * @param pB  the frequency of B
     *
     * @return rw
     */
    private static double r2(double pAB, double pA, double pB) {
        return r2(rw(pAB, pA, pB));
    }

    /**
     * Compute r2 for a SNPpair (it should be polarized)
     *
     * @param pair polarized pair
     *
     * @return rw
     */
    public static double r2(SNPpair pair) {
        double pA = pair.getMajorAF_A();
        double pB = pair.getMajorAF_B();
        double pAB = pair.getFreq_AB();
        return r2(pAB, pA, pB);
    }

    /**
     * Compute maximumR2 with the frequencies
     *
     * @param smaller the smaller frequency
     * @param higher  the higher frequency
     *
     * @return the maximum R2
     */
    private static double maxR2(double smaller, double higher) {
        return (smaller * (1 - higher)) / ((1 - smaller) * higher);
    }

    /**
     * Compute the maximumR2 for a SNPpair (it should be polarized)
     *
     * @param pair polarized pair
     *
     * @return maxR2
     */
    public static double maxR2(SNPpair pair) {
        // Because SNPpair is polarized, it is easier
        double pA = pair.getMajorAF_A();
        double pB = pair.getMajorAF_B();
        if (pA == pB) {
            return 1;
        }
        return (pA < pB) ? maxR2(pA, pB) : maxR2(pB, pA);
    }

    /**
     * Compute the r statistics for a SNPpair with a threshold for minimum counts and chi-squared
     * quantile for the significance of maxR2
     *
     * @param pair           polarized pair
     * @param chiSqrQuantile the chi-squared quantile to asses significance of maxR2
     *
     * @return array with the values if rmax is upper the threshold; null otherwise
     */
    public static double[] rStatistics(SNPpair pair, double chiSqrQuantile) {
        // maf for A and B
        double pA = pair.getMajorAF_A();
        double pB = pair.getMajorAF_B();
        // compute the maximum r2 that can be reached
        double maxR2 = 1;
        if (pA != pB) {
            // logger.debug("Differennt pA=", pA, " and pB=", pB, "; for ", pair);
            maxR2 = (pA < pB) ? maxR2(pA, pB) : maxR2(pB, pA);
        }
        // if the maximum r2 is not significant, it is not worthy to compute the rest of the stuff
        if (!isSignificantR2(maxR2, pair.getTotalCounts(), chiSqrQuantile)) {
            return null;
        }
        // compute rw and r2 (fast way)
        double pAB = pair.getFreq_AB();
        double rw = rw(pAB, pA, pB);
        if (rw == 0) {
            return new double[] {maxR2, 0, 0, 0};
        }
        double r2 = r2(rw);
        return new double[] {maxR2, r2, r2 / maxR2, rw};
    }

    /**
     * Check if some value of r2 is significant based on the distribution of chi-square based on the
     * number of haplotypes
     *
     * @param r2                 the r2 value to test
     * @param numberOfHaplotypes the number of haplotypes
     * @param chiSqrQuantile     the chi-squared quantile to asses significance of maxR2
     *
     * @return    <code>true</code> if r2 is significant; <code>false</code> otherwise
     */

    public static boolean isSignificantR2(double r2, int numberOfHaplotypes,
            double chiSqrQuantile) {
        return r2 >= significantThreshold(numberOfHaplotypes, chiSqrQuantile);
    }

    /**
     * Return the minimum r2 that could be significan for the current number of haplotypes
     *
     * @param numberOfHaplotypes the number of haplotypes
     * @param chiSqrQuantile     the chi-squared quantile to asses significance of maxR2
     *
     * @return the minimum r2 that can be considered significantly different from 0
     */
    public static double significantThreshold(int numberOfHaplotypes, double chiSqrQuantile) {
        return chiSqr.inverseCumulativeProbability(chiSqrQuantile) / numberOfHaplotypes;
    }

}

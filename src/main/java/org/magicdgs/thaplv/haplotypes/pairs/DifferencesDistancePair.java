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

import org.magicdgs.thaplv.utils.AlleleUtils;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Set;

/**
 * Distance for a pair of samples based on the number of differences between them. To account for
 * the uncertainty of missing calls, half of the value of them in the pair (in one or both of the
 * samples) is added to the number of differences.
 *
 * The number of sites considered does not count as distance, but could be retieved using {@link
 * #getNumberOfSites()} to weight.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class DifferencesDistancePair extends DistancePair {

    // logger for the class
    private static Logger logger = LogManager.getLogger(DifferencesDistancePair.class);

    // the number of differences
    private int differences = 0;

    // the number of positions considered
    private int nPositions = 0;

    // the number of missing (in one or both)
    private int nMissing = 0;

    /**
     * Initialize the pair of samples with their names
     *
     * @param sample1 name for the first sample
     * @param sample2 name for the second sample
     */
    public DifferencesDistancePair(final String sample1, final String sample2) {
        super(sample1, sample2);
    }

    /**
     * Get the number of real differences included, without accounting for the missing calls
     */
    public int getNumberOfDifferences() {
        return differences;
    }

    /**
     * Get the number of sites included
     */
    public int getNumberOfSites() {
        return nPositions;
    }

    /**
     * Get the number of missing calls included
     */
    public int getNumberOfMissing() {
        return nMissing;
    }

    /**
     * Get the distance for the pair. No differences counts as 0, differences as 1 and missing data
     * as 0.5
     *
     * @return the distance for the pair; {@link Double#NaN} if not sites (including missing) were
     * included
     */
    @Override
    public double getDistance() {
        return (nPositions + nMissing == 0) ? Double.NaN : differences + (0.5D * nMissing);
    }

    /**
     * Get the number of differences per called site (omitting missing)
     *
     * @return the number of differences per site; {@link Double#NaN} if no called sites were
     * included
     */
    public double getDifferencesPerSite() {
        return (nPositions == 0) ? Double.NaN : (double) differences / nPositions;
    }

    /**
     * Add a pair of haplotypes to the distance computation
     *
     * WARNING: only works for haploid genotypes.
     *
     * @param haplotype1 the first haplotype
     * @param haplotype2 the second haplotype
     * @param checkNames if {@code true}, check the names in the genotypes
     *
     * @throws IllegalArgumentException if the sample names does not match the ones in the pair
     */
    protected void add(final Genotype haplotype1, final Genotype haplotype2, boolean checkNames) {
        Utils.nonNull(haplotype1, "null haplotype1");
        Utils.nonNull(haplotype2, "null haplotype2");
        if (checkNames && !containNames(haplotype1.getSampleName(), haplotype2.getSampleName())) {
            logger.debug("BUG: adding {} vs. {} pairs to {}", haplotype1.getSampleName(),
                    haplotype2.getSampleName(), getPairNames());
            throw new IllegalArgumentException(
                    "Trying to compute a difference between haplotypes that are not in this pair");
        }
        // only compute the difference between called positions
        // TODO: this is a putative broken change with respect to the previous repo
        final int difference = AlleleUtils.scaledDifference(
                haplotype1.getAllele(0), haplotype2.getAllele(0));
        switch (difference) {
            case 1: // missing
                nMissing++;
                break; // no site is count
            case 2: // difference
                differences++;
            case 0: // no difference
                nPositions++; // positions are count for both
                break;
            default:
                // this should never happen
                throw new RuntimeException("Unreacheable code");

        }
    }

    /**
     * Add a pair of haplotypes to the distance computation, checking for the names by default
     *
     * WARNING: only works for haploid genotypes.
     *
     * @param haplotype1 the first haplotype
     * @param haplotype2 the second haplotype
     *
     * @throws IllegalArgumentException if the sample names does not match the ones in the pair
     */
    public void add(final Genotype haplotype1, final Genotype haplotype2) {
        add(haplotype1, haplotype2, true);
    }

    /**
     * Add a variant (which should contains both sample1 and sample2) to the distance computation
     *
     * WARNING: only works for haploid genotypes
     *
     * @param variant the variant which contains the haplotypes for the samples
     *
     * @throws java.lang.IllegalArgumentException if both sample names are not in the variant
     */
    public void add(final VariantContext variant) {
        Utils.nonNull(variant, "null variant");
        // false because it is already checked that this corresponds to this pair
        add(variant.getGenotype(sample1), variant.getGenotype(sample2), false);
    }

    /**
     * Add a difference of haplotypes to the distance computation, checking for the names by
     * default
     *
     * WARNING: only works for haploid genotypes.
     *
     * @param haplotype  the sample haplotype
     * @param checkNames {@code true} for check the names in the genotypes; {@code false }otherwise
     *
     * @throws IllegalArgumentException if the sample name does not match the ones in the pair
     */
    protected void addReference(final Genotype haplotype, final boolean checkNames) {
        Utils.nonNull(haplotype, "null haplotype");
        if (checkNames && !containSample(haplotype.getSampleName())) {
            logger.debug("BUG: Computing reference difference for ", haplotype.getSampleName(), " in ",
                    getPairNames());
            throw new IllegalArgumentException(
                    "Trying to compute a difference between haplotypes that are not in this pair");
        }
        // only compute differences if it is called
        // TODO: this is a putative broken change with respect to the previous repo
        if (haplotype.isCalled()) {
            differences += haplotype.isHomRef() ? 0 : 1;
            nPositions++;
        } else {
            nMissing++;
        }
    }

    /**
     * Compute the differences for a sample haplotype (genotypes should be homozygous) w.r.t.
     * reference, checking for the names.
     *
     * WARNING: only works for haploid genotypes.
     *
     * @param sample the sample haplotype
     *
     * @throws IllegalArgumentException if the sample names does not match the ones in the pair
     */
    public void addReference(final Genotype sample) {
        addReference(sample, true);
    }

    /**
     * Compute the differences for an haplotype (genotypes should be homozygous) w.r.t. reference.
     *
     * WARNING: only works for haploid genotypes
     *
     * @param variant the variant which contains the haplotypes for the samples
     *
     * @throws IllegalArgumentException if the sample names does not match the ones in thepair
     */
    public void addReference(final VariantContext variant) {
        Utils.nonNull(variant);
        final Set<String> names = variant.getSampleNames();
        // TODO: this is a putative broken change with respect to the previous repo
        // only checks once
        final boolean containsSample1 = names.contains(sample1);
        final boolean containsSample2 = names.contains(sample2);
        if (containsSample1 && containsSample2) {
            logger.debug("BUG: Computing reference difference for {}:{} in {}.",
                    variant.getContig(), variant.getStart(), getPairNames());
            throw new IllegalArgumentException(
                    "Only one of the samples in a reference computation could be included in the pair.");
        }
        if (containsSample1) {
            addReference(variant.getGenotype(sample1), false);
        } else if (containsSample2) {
            addReference(variant.getGenotype(sample2), false);
        } else { // the lack of this else was a bug in previous version
            logger.debug("BUG: Computing reference difference for {}:{} in {}",
                    variant.getContig(), variant.getStart(), getPairNames());
            throw new IllegalArgumentException(
                    "None of the samples in the pair is included for reference computation");
        }
    }

    public String toString() {
        return "[" + getPairNames() + " -> differences=" + differences + "; sites=" + nPositions
                + "; missing="
                + nMissing + "]";
    }
}

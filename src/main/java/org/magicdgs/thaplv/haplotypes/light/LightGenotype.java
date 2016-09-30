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

package org.magicdgs.thaplv.haplotypes.light;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;

import java.io.Serializable;
import java.util.HashSet;
import java.util.List;
import java.util.stream.Stream;

/**
 * Class for a light-weight representation of genotypes (when storage in memory is needed).
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 * @deprecated AlleleVector is a better light representation
 */
@Deprecated
public class LightGenotype implements Serializable {
    private static final long serialVersionUID = 1L;

    private final SNP[] orderedSNPlist;
    private final String contig;
    private final int pos;

    /**
     * Public constructor from a variant
     *
     * @param variant the variant to convert
     */
    public LightGenotype(final VariantContext variant) {
        this(variant.getContig(), variant.getStart(), variant.getGenotypes(),
                variant.getReference(), variant.getAlternateAlleles());
    }

    /**
     * Public constructor with all the parameters necessary
     *
     * @param contig       the contig for the genotype
     * @param position     the position for the genotype
     * @param genotypes    the genotypes
     * @param ref          the reference allele
     * @param alternatives the alternative alleles
     */
    public LightGenotype(final String contig, final int position, final GenotypesContext genotypes,
            final Allele ref, final List<Allele> alternatives) {
        this.contig = contig;
        this.pos = position;
        orderedSNPlist = new SNP[genotypes.size()];
        switch (alternatives.size()) {
            case 0:
                throw new IllegalArgumentException(
                        "Trying to obtain a light genotype from a variant with only reference");
            case 1:
                initBiallelic(genotypes);
                break;
            default:
                initMultiAllelic(genotypes, ref, alternatives);
                break;
        }
    }

    /**
     * Constructor for testing purposes.
     *
     * @param contig         the contig for the genotype
     * @param pos            the position for the genotype
     * @param orderedSNPlist the list of SNPs
     */
    @VisibleForTesting
    public LightGenotype(final String contig, final int pos, final SNP[] orderedSNPlist) {
        this.contig = contig;
        this.pos = pos;
        this.orderedSNPlist = orderedSNPlist;
    }

    /**
     * Initialize the SNPlist from a ref/alt biallelic SNP. The reference will be considered the
     * SNP.A
     *
     * @param genotypes the genotypes
     */
    private void initBiallelic(final GenotypesContext genotypes) {
        int i = 0;
        for (final Genotype geno : genotypes) {
            if (geno.isHomRef()) {
                orderedSNPlist[i++] = SNP.A;
            } else if (geno.isHomVar()) {
                orderedSNPlist[i++] = SNP.a;
            } else {
                orderedSNPlist[i++] = SNP.N;
            }
        }
    }

    /**
     * Initialize the SNPlist from various alternatives. The one for the first sample will be
     * considered the SNP.A
     *
     * @param genotypes    the genotypes
     * @param ref          the reference allele
     * @param alternatives the alternative alleles
     */
    private void initMultiAllelic(final GenotypesContext genotypes, final Allele ref,
            final List<Allele> alternatives) {
        // create a set with the alleles
        final HashSet<Allele> alleles = new HashSet<>();
        alleles.add(ref);
        alleles.addAll(alternatives);
        Allele snpA = null;
        Allele snpa = null;
        int i = 0;
        for (final Genotype geno : genotypes) {
            // if it is not homozygous, is a missing SNP
            if (!geno.isHom()) {
                orderedSNPlist[i++] = SNP.N;
            } else {
                final Allele current = geno.getAllele(0);
                if (snpA == null) {
                    snpA = current;
                    alleles.remove(snpA);
                }
                if (current.equals(snpA)) {
                    orderedSNPlist[i++] = SNP.A;
                } else if (snpa != null && alleles.contains(current)) {
                    throw new IllegalArgumentException(
                            "Trying to construct a light genotype from a multi-allelic variant when is not. SNP_A="
                                    + snpA + "; SNP_a=" + snpa + "; rest=" + alleles);
                } else if (snpa == null) {
                    snpa = current;
                    alleles.remove(snpa);
                    orderedSNPlist[i++] = SNP.a;
                } else if (current.equals(snpa)) {
                    orderedSNPlist[i++] = SNP.a;
                } else {
                    throw new RuntimeException("Unreachable code");
                }
            }
        }
    }

    /**
     * Get the position for the genotypes
     *
     * @return the position
     */
    public int getPosition() {
        return this.pos;
    }

    /**
     * Get the contig reference for the genotypes
     *
     * @return the contig
     */
    public String getContig() {
        return this.contig;
    }

    /**
     * Get the number of alleles A
     *
     * @return the count
     */
    public int getNumberOfA() {
        return getNumberOf(SNP.A);
    }

    /**
     * Get the number of alleles a
     *
     * @return the count
     */
    public int getNumberOfa() {
        return getNumberOf(SNP.a);
    }

    /**
     * Get the number of missing
     *
     * @return the count
     */
    public int getNumberOfMissing() {
        return getNumberOf(SNP.N);
    }

    /**
     * Check if the counts for any of the alleles is only one. An invariant site is not singleton
     *
     * @return true if the count for one of the alleles is 1; false otherwise
     */
    public boolean isSingleton() {
        return getNumberOfA() == 1 || getNumberOfa() == 1;
    }

    /**
     * Compute the number of SNPs of certain category
     *
     * @param snp the SNP to count
     *
     * @return the count of "snp"
     */
    private int getNumberOf(final SNP snp) {
        return (int) Stream.of(orderedSNPlist).filter(x -> x == snp).count();
    }

    /**
     * Get the number of genotypes
     *
     * @return the number of genotypes
     */
    public int size() {
        return orderedSNPlist.length;
    }

    /**
     * Get the genotype at index position. The first index is 0
     *
     * @param index the index to retrieve
     *
     * @return the genotype from the index sample
     */
    public SNP getGenotypeAt(final int index) {
        return this.orderedSNPlist[index];
    }

    /**
     * String representation of the LightGenotype
     *
     * @return the string representation
     */
    public String toString() {
        final StringBuilder builder = new StringBuilder();
        builder.append('[');
        for (final SNP snp : orderedSNPlist) {
            builder.append(snp);
        }
        builder.append(']');
        return builder.toString();
    }

    /**
     * Enum to light storage
     */
    public static enum SNP {
        A,
        a,
        N
    }

}

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

package org.magicdgs.thaplv.haplotypes;

import org.magicdgs.thaplv.utils.AlleleUtils;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Light-weight representation of alleles from genotypes as vector.
 *
 * The vector is polarized by allele frequency in decreasing order. In case of ties, the alleles are
 * ordered by {@link Allele#compareTo(Object)}.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class AlleleVector {

    // the list of alleles encoded as the index in the allele set
    private final int[] encodedAlleles;

    // the alleles included in this genotype
    private final List<Allele> alleleSet;

    /**
     * Constructs an allele vector.
     *
     * @param encodedAlleles the alleles encoded as 0, 1, ... for know alleles (from most frequent
     *                       to less frequent) and -1 for unknown ({@link Allele#NO_CALL}).
     * @param alleleSet      list of unique alleles in {@link #encodedAlleles}, where the index
     *                       correspond to the encoding
     */
    @VisibleForTesting
    AlleleVector(final int[] encodedAlleles, final List<Allele> alleleSet) {
        this.encodedAlleles = encodedAlleles;
        this.alleleSet = Collections.unmodifiableList(alleleSet);
    }

    /**
     * Constructs an allele vector from an array of alleles.
     *
     * @param alleles the array of alleles.
     */
    public AlleleVector(final Allele[] alleles) {
        Utils.nonNull(alleles, "null alleles");
        Utils.validateArg(alleles.length != 0, "empty alleles");
        // count the alleles
        final Map<Allele, Integer> alleleCounts =
                AlleleUtils.getAlleleCounts(Stream.of(alleles), false);
        // allele counts polarizing by major allele
        final Comparator<Map.Entry<Allele, Integer>> reversePolirized =
                Map.Entry.comparingByValue(Comparator.reverseOrder());
        this.alleleSet = Collections.unmodifiableList(
                alleleCounts.entrySet().stream()
                        // sort by counts, the major allele first; solve ties with allele comparator
                        .sorted(reversePolirized.thenComparing(Map.Entry.comparingByKey()))
                        .map(Map.Entry::getKey).collect(Collectors.toList()));
        this.encodedAlleles = Stream.of(alleles).mapToInt(alleleSet::indexOf).toArray();
    }

    /**
     * Constructs an allele vector from a list of haploid genotypes.
     *
     * @param genotypes haploid genotypes.
     */
    public AlleleVector(final List<Genotype> genotypes) {
        this(Utils.nonNull(genotypes, "null genotypes")
                .stream().map(g -> g.getAllele(0)).toArray(Allele[]::new));
    }

    /**
     * Constructs an allele vector from a haploid samples contained in a {@link VariantContext}.
     *
     * @param variant the variant context.
     */
    public AlleleVector(final VariantContext variant) {
        this(Utils.nonNull(variant, "null variant").getGenotypes());
    }

    /**
     * Returns the number of samples in the vector.
     */
    public int size() {
        return encodedAlleles.length;
    }

    /**
     * Gets all sample alleles contained in the vector.
     *
     * Note: calling this generates a fresh array.
     */
    public Allele[] getAlleles() {
        return IntStream.of(encodedAlleles).mapToObj(this::alleleAtIndex).toArray(Allele[]::new);
    }

    /**
     * Gets the allele for a certain sample.
     *
     * @param index the index of the sample in the vector.
     *
     * @return the allele at that position.
     *
     * @throws IndexOutOfBoundsException if the index is out of range.
     */
    public Allele getAlleleAt(final int index) {
        if (index >= encodedAlleles.length) {
            throw new IndexOutOfBoundsException("index out of range: " + index);
        }
        return alleleAtIndex(encodedAlleles[index]);
    }

    // helper method for handling no-call alleles at index 1
    private Allele alleleAtIndex(final int i) {
        return (i == -1) ? Allele.NO_CALL : alleleSet.get(i);
    }

    /** Gets the number of unique alleles for this vector. */
    public int numberOfAlleles() {
        return alleleSet.size();
    }

    /**
     * Returns the same vector if it could be represented as a biallelic already (1 or two alleles);
     * if not, it use only the two major alleles.
     */
    public AlleleVector asBiallelic() {
        if (numberOfAlleles() < 3) {
            return this;
        }
        final int[] biallelic = IntStream.of(encodedAlleles).map(i -> (i < 2) ? i : -1).toArray();
        return new AlleleVector(biallelic, alleleSet.subList(0, 2));
    }

    /**
     * Returns {@code true} if the two allele vectors have equals allelic states for all the
     * samples; {@code false} otherwise.
     */
    public boolean equalsState(final AlleleVector other) {
        if (!this.alleleSet.equals(other.alleleSet)) {
            return false;
        }
        return Arrays.equals(this.encodedAlleles, other.encodedAlleles);
    }

    /**
     * WARNING: equals does not include the allele information. Use {@link
     * #equalsState(AlleleVector)} in case of the allele equality is important.
     *
     * This is a convenience for use equals/hashCode to cache results.
     */
    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (!(o instanceof AlleleVector)) {
            return false;
        }
        return Arrays.equals(encodedAlleles, ((AlleleVector) o).encodedAlleles);
    }

    /**
     * Only the encoded alleles are included in the hashCode.
     *
     * @see {@link #equals(Object)} for more information.
     */
    @Override
    public int hashCode() {
        // TODO: putative performance improvement -> cache the hashCode
        return Arrays.hashCode(encodedAlleles);
    }

    // for testing purposes
    public String toString() {
        return String.format("Vector of %s: %s", alleleSet, Arrays.toString(encodedAlleles));
    }
}

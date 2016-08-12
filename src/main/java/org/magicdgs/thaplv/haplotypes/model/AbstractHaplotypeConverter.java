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

package org.magicdgs.thaplv.haplotypes.model;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Abstract {@link HaplotypeConverter}, with default methods suitable for most of the
 * implementations of {@link HaplotypeModel}
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public abstract class AbstractHaplotypeConverter implements HaplotypeConverter {

    /** Unmodifiable representation of missing alleles for the output ploidy provided */
    protected final List<Allele> missingGenotypeAllele;

    /**
     * Construct an abstract converter. The ouptut ploidy indicates the ploidy for the returned
     * {@link Genotype}.
     *
     * WARNING: {@link AbstractHaplotypeConverter} only handle haploid and diploids
     */
    protected AbstractHaplotypeConverter(final int outputPloidy) {
        Utils.validateArg(!(outputPloidy != 1 && outputPloidy != 2), "outputPloidy != 1 or 2");
        missingGenotypeAllele = Collections.unmodifiableList(
                Collections.nCopies(outputPloidy, Allele.NO_CALL));
    }

    /**
     * Get the haplotype for a single genotype. If the genotype is already an haploid, return it as
     * it is (and duplicate alleles if the output ploidy is not the same). If not, first check if
     * it is considered missing by the model ({@link #isConsideredMissing(Genotype)}) and if it is
     * not, try to get the alleles using {@link #getHaplotypeAllele(List)}.
     *
     * The returned genotype will have {@link #getOutputPlotidy()} alleles, and phased.
     */
    public Genotype apply(final Genotype genotype) {
        // early termination for already haploid individuals
        if (genotype.getAlleles().size() == 1) {
            return (getOutputPlotidy() == 1) ? genotype :
                    new GenotypeBuilder(genotype).phased(true)
                            .alleles(Collections.nCopies(getOutputPlotidy(), genotype.getAllele(0)))
                            .make();
        }
        final GenotypeBuilder builder = new GenotypeBuilder(genotype).phased(true);
        if (isConsideredMissing(genotype)) {
            return builder.alleles(missingGenotypeAllele).make();
        } else {
            return builder.alleles(Collections.nCopies(missingGenotypeAllele.size(),
                    getHaplotypeAllele(genotype.getAlleles())))
                    .make();
        }
    }

    /** {@inheritDoc} */
    // The output ploidy is the length of missingGenotypeAlleles
    public int getOutputPlotidy() {
        return missingGenotypeAllele.size();
    }

    /**
     * Convert a list of alleles from a single individual into the alleles for the haplotype.
     *
     * WARNING: {@link AbstractHaplotypeConverter} only consider genotypes with 2 alleles or less.
     *
     * @param alleles the alleles from the individual
     *
     * @return the allele in the haplotype from the individual
     */
    protected Allele getHaplotypeAllele(final List<Allele> alleles) {
        // converting the allele to a set to quick check if they are not equals
        final Set<Allele> alleleSet = new HashSet<>(alleles);
        switch (alleleSet.size()) {
            // if the set is 1, return the allele
            case 1:
                return alleles.get(0);
            case 2:
                return getHaplotypeAllele(alleles.get(0), alleles.get(1));
            default:
                throw new IllegalArgumentException(
                        "This software cannot handle genotypes with more than 2 alleles");
        }
    }

    /**
     * Is this genotype considered missing by this model? It should return always true if the
     * genotype is completely missing.
     *
     * Behaviour of this method may be unexpected for non-diploid individuals in implementations.
     *
     * @param genotype the genotype to test
     *
     * @return {@code true} if it is considered missing; {@code false} otherwise
     */
    protected abstract boolean isConsideredMissing(final Genotype genotype);

    /**
     * Get the haplotype for two alleles that are different and are not considered missing by this
     * HaplotypeConverter.
     *
     * Implementations don't need to handle situations where alleles do not pass the  the
     * {@link #isConsideredMissing(Genotype)} filter.
     *
     * @param allele1 the first allele
     * @param allele2 the second allele
     *
     * @return the allele that configures the haplotype
     */
    protected abstract Allele getHaplotypeAllele(final Allele allele1, final Allele allele2);
}

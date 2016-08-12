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
import org.apache.logging.log4j.Logger;

/**
 * {@link HaplotypeConverter} for individuals that come from a back-cross. We assume that the
 * reference allele is the one defined as reference in the {@link htsjdk.variant.variantcontext.VariantContext}.
 *
 * Two models could be used:
 *
 * 1) Conservative: homozygous alternative are considered missing. The idea behind this behaviour
 * is that the crossed reference genome is completely knonw (pure homozygous/haploid genome) so
 * this kind of variants represent sequencing errors.
 *
 * 2) Non-conservative: homozygous alternative are considered alternative. The idea behind this
 * behaviour is that some variability could be found in the crossed reference, not tracked in the
 * FASTTA file, and thus one of the alleles comes from the reference and the other from the
 * haplotype.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public final class BackCross extends AbstractHaplotypeConverter {

    private final boolean conservative;

    /**
     * Constructor assuming a conservative or non-conservative model
     *
     * @param outputPloidy final ploidy for the returned haplotype
     * @param conservative should homozygous alternative being considered missing?
     */
    public BackCross(final int outputPloidy, final boolean conservative) {
        super(outputPloidy);
        this.conservative = conservative;
    }

    /**
     * Constructor for a conservative back-cross model
     *
     * @param finalPloidy final ploidy for the haplotype
     */
    public BackCross(final int finalPloidy) {
        this(finalPloidy, true);
    }

    @Override
    public void log(final Logger logger) {
        super.log(logger);
        logger.warn("Reference homozygous are considered reference");
        logger.warn("Heterozygous positions are considered alternative");
        logger.warn("Homozygous alternative are considered {}",
                (conservative) ? "missing" : "alternative");
        logger.warn("Heterozygous non-reference are considered missing");
    }

    /**
     * In the back-cross model a genotype is considered missing if it is heterozygous reference or
     * pseudo-called (one of them is called and the other is not).
     *
     * In the conservative model, also homozygous alternative are missing
     */
    @Override
    public boolean isConsideredMissing(final Genotype genotype) {
        return genotype.isNoCall() // first check if it is no call
                // heterozygous reference are never considered
                || (genotype.isHetNonRef())
                // if it is the conservative model, non homozygous alternative
                || (conservative && genotype.isHomVar());
    }

    /**
     * In the back-cross model, the non-reference allele is the one returned, and if one of them is
     * missing it can't be called.
     *
     * In the conservative mode, one of the alleles should be reference to call the allele.
     */
    @Override
    protected Allele getHaplotypeAllele(final Allele allele1, final Allele allele2) {
        // do this checking only if it is the conservative mode
        if (conservative // both alleles should be called and non-reference
                && allele1.isNonReference() && allele2.isNonReference()) {
            return Allele.NO_CALL;
        }
        // if one of them is no-called
        if (allele1.isNoCall() || allele2.isNoCall()) {
            return Allele.NO_CALL;
        }
        // the non-reference is the call
        return (allele1.isNonReference()) ? allele1 : allele2;
    }

}

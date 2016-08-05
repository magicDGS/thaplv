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

package org.magicdgs.thaplv.haplotypes.model;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import org.apache.logging.log4j.Logger;

/**
 * {@link HaplotypeConverter} for individuals that come from an haploid sequencing or inbred lines.
 *
 * Heterozygous genotypes will be considered missing, because either is a sequencing error or a
 * segregating site in the inbred line.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public final class InbredLine extends AbstractHaplotypeConverter {

    /**
     * Constructor for an haploid/inbred line converter. The ouptut ploidy indicates the ploidy for
     * the returned {@link Genotype} for the converter.
     */
    protected InbredLine(final int outputPloidy) {
        super(outputPloidy);
    }

    @Override
    public void log(final Logger logger) {
        super.log(logger);
        logger.warn("Only homozygous positions are considered as called");
        logger.warn("Heterozygous positions are considered missing");
    }

    /**
     * Return {@code true} if the genotype is not homozygous. If one of the alleles is missing,
     * return {@code false}
     */
    @Override
    public boolean isConsideredMissing(final Genotype genotype) {
        return !genotype.isHom() && !genotype.isMixed();
    }

    /**
     * Both alleles should be either equals or one of them no-called, because otherwise {@link
     * #isConsideredMissing(Genotype)} will return false for this converter.
     *
     * Return the allele that is called (if both are called, they should be the same)
     */
    @Override
    protected Allele getHaplotypeAllele(final Allele allele1, final Allele allele2) {
        return (allele1.isCalled()) ? allele1 : allele2;
    }
}

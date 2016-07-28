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

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * Haplotype model that stores information for converting {@link htsjdk.variant.variantcontext.Genotype}
 * to haplotypes.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public enum HaplotypeModel {

    INBRED_LINE(new InbredLine(1), new InbredLine(2)),
    HAPLOID(INBRED_LINE.haploid, INBRED_LINE.diploid),
    BACK_CROSS(new BackCross(1), new BackCross(2)),
    DONT_CHECK(DontCheckHaplotypeConverter.SINGLETON, DontCheckHaplotypeConverter.SINGLETON);

    // return haploid haplotypes
    private final VariantHaplotypeConverter haploid;
    // return diploid haplotypes
    private final VariantHaplotypeConverter diploid;

    HaplotypeModel(final HaplotypeConverter haploid, final HaplotypeConverter diploid) {
        this.haploid = HaplotypeConverter.getVariantConverter(haploid);
        this.diploid = HaplotypeConverter.getVariantConverter(diploid);
    }

    HaplotypeModel(final VariantHaplotypeConverter haploid,
            final VariantHaplotypeConverter diploid) {
        this.haploid = haploid;
        this.diploid = diploid;
    }

    /**
     * Get the converter for a variant for this model. Only haploid/diploids are allowed
     */
    public static VariantHaplotypeConverter getVariantHaplotypeConverter(
            final HaplotypeModel model, final int ploidy) {
        Utils.nonNull(model, "null model");
        switch (ploidy) {
            case 1:
                return model.haploid;
            case 2:
                return model.diploid;
            default:
                throw new IllegalArgumentException("Ploidy of " + ploidy + " is not allowed");
        }
    }

    /**
     * {@link VariantHaplotypeConverter} for don't check genotypes
     */
    private static final class DontCheckHaplotypeConverter implements VariantHaplotypeConverter {

        /** Use always the singleton for avoid constructing new objects */
        private final static DontCheckHaplotypeConverter SINGLETON =
                new DontCheckHaplotypeConverter();

        @Override
        public void log(final Logger logger) {
            logger.info("Genotypes won't be checked for haploid state.");
            logger.warn("Not checking genotypes is discouraged if the haplotypes were obtained with diploid callers or not converted before.");
            logger.warn("Most of the tools assume haploid individuals (only the first allele in the Genotypes is considered) and results could be wrong if they are not.");
        }

        @Override
        public VariantContext apply(final VariantContext variant) {
            return variant;
        }
    }
}

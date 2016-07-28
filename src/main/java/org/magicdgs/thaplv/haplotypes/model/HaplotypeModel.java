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

import java.util.function.Function;

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
    DONT_CHECK(v -> v, v -> v);

    // return haploid haplotypes
    private final Function<VariantContext, VariantContext> haploid;
    // return diploid haplotypes
    private final Function<VariantContext, VariantContext> diploid;

    HaplotypeModel(final HaplotypeConverter haploid, final HaplotypeConverter diploid) {
        this.haploid = HaplotypeConverter.getVariantConverter(haploid);
        this.diploid = HaplotypeConverter.getVariantConverter(diploid);
    }

    HaplotypeModel(final Function<VariantContext, VariantContext> haploid,
            final Function<VariantContext, VariantContext> diploid) {
        this.haploid = haploid;
        this.diploid = diploid;
    }

    /**
     * Get the converter for a variant for this model. Only haploid/diploids are allowed
     */
    public static Function<VariantContext, VariantContext> getVariantConverter(
            final HaplotypeModel model, final int ploidy) {
        switch (ploidy) {
            case 1:
                return model.haploid;
            case 2:
                return model.diploid;
            default:
                throw new IllegalArgumentException("Ploidy of " + ploidy + " is not allowed");
        }
    }

}

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

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * Function that converts a {@link Genotype} with a diploid individual into haploid. When a {@link
 * Genotype} with an already haploid individual is feed into the function, a new copy of the
 * variant
 * will be returned.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public interface HaplotypeConverter extends Function<Genotype, Genotype> {

    /** Get the output ploidy for the converter. */
    int getOutputPlotidy();

    /**
     * Log the characteristics of this haplotype converter in the provided logger.
     * Default implementation only outputs the name of the HaplotypeConverter.
     */
    default void log(final Logger logger) {
        logger.debug("Using {} as haplotype model", this.getClass().getSimpleName());
    }

    /**
     * Get a function to map a {@link VariantContext} with haploid/diploids to a {@link
     * VariantContext} with haplotypes, using the provided {@link HaplotypeConverter}
     */
    public static VariantHaplotypeConverter getVariantConverter(
            final HaplotypeConverter converter) {
        return new VariantHaplotypeConverter() {

            @Override
            public void log(final Logger logger) {
                converter.log(logger);

            }

            @Override
            public VariantContext apply(final VariantContext variant) {
                Utils.nonNull(variant);
                final VariantContextBuilder builder = new VariantContextBuilder(variant);
                final List<Genotype> genotypes = variant.getGenotypes().stream()
                        .map(converter).collect(Collectors.toList());
                builder.genotypes(genotypes);
                return builder.make();
            }
        };
    }
}

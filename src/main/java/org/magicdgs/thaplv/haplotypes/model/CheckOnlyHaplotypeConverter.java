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

import org.magicdgs.thaplv.cmd.ThaplvArgumentDefinitions;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.util.concurrent.atomic.AtomicInteger;

/**
 * Singleton {@link VariantHaplotypeConverter} for check every {@link #samplingFrequency} variants
 * if the Genotypes are real haploids.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
class CheckOnlyHaplotypeConverter implements VariantHaplotypeConverter {

    // sample 1 variant for every 100 encountered
    @VisibleForTesting
    static final int samplingFrequency = 100;

    @VisibleForTesting
    final AtomicInteger counter = new AtomicInteger(0);

    /** Use always the singleton for avoid constructing new objects */
    private final static CheckOnlyHaplotypeConverter SINGLETON = new CheckOnlyHaplotypeConverter();

    /** Get this variant haplotype converter. */
    public static CheckOnlyHaplotypeConverter getSingleton() {
        return SINGLETON;
    }

    /** Can't be instantiated. */
    private CheckOnlyHaplotypeConverter() {}

    @Override
    public void log(final Logger logger) {
        logger.info("Haploid genotypes are expected an no transformation will be done.");
        logger.info("Genotypes will be checked every {} variants for haploid state.",
                samplingFrequency);
        logger.warn(
                "No transformation is discouraged if the haplotypes were obtained with diploid callers or not converted before.");
        logger.warn(
                "Most of the tools assume haploid individuals (only the first allele in the Genotypes is considered) and results could be wrong if they are not.");
    }

    @Override
    public VariantContext apply(final VariantContext variant) {
        if (counter.getAndIncrement() % samplingFrequency == 0) {
            for (final Genotype g : variant.getGenotypes()) {
                if (g.getAlleles().stream().distinct().count() != 1) {
                    throw new UserException.BadArgumentValue(
                            ThaplvArgumentDefinitions.HAPLOTYPE_MODEL_LONG,
                            HaplotypeModel.CHECK_ONLY.toString(),
                            String.format(
                                    "found no haploid individual (%s %s) at %s:%s-%s. Please, convert to haplotypes before using this option.",
                                    g.getSampleName(), g.getGenotypeString(),
                                    variant.getContig(), variant.getStart(), variant.getEnd()));
                }
            }
        }
        return variant;
    }
}

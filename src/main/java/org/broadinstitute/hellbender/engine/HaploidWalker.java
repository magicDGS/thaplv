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

package org.broadinstitute.hellbender.engine;

import org.magicdgs.thaplv.cmd.argumentcollections.HaplotypeModelArgumentCollection;
import org.magicdgs.thaplv.cmd.argumentcollections.HaplotypeModelNoPloidyArgumentCollection;
import org.magicdgs.thaplv.cmd.argumentcollections.HaplotypeModelWithPloidyArgumentCollection;
import org.magicdgs.thaplv.haplotypes.model.VariantHaplotypeConverter;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.engine.filters.VariantFilter;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.stream.StreamSupport;

/**
 * {@link VariantWalker} that transform the {@link VariantContext} into haploid.
 *
 * The variants passed to {@link #apply(VariantContext, ReadsContext, ReferenceContext,
 * FeatureContext)} are guarantee to have haploid individuals of undetermined ploidy except if
 * {@link #allowsDontCheck()} is {@code true}. Thus, using the first allele from {@link
 * htsjdk.variant.variantcontext.Genotype#getAllele(int)} will be enough for compute statistics.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public abstract class HaploidWalker extends VariantWalker {

    // the haplotypes models
    @ArgumentCollection(doc = "Haplotype models")
    public HaplotypeModelArgumentCollection haplotypeModelArgumentCollection =
            (requiresOutputPloidy())
                    ? new HaplotypeModelWithPloidyArgumentCollection(allowsDontCheck())
                    : new HaplotypeModelNoPloidyArgumentCollection(allowsDontCheck());

    /**
     * Does this tool requires an output ploidy for output VCF?
     *
     * @return {@code true} if it will be output; {@link false} otherwise
     */
    protected abstract boolean requiresOutputPloidy();

    /**
     * Returns {@code true} if {@link org.magicdgs.thaplv.haplotypes.model.HaplotypeModel#DONT_CHECK}
     * model is allowd; {@link false} otherwise.
     */
    protected abstract boolean allowsDontCheck();

    // the variant converter for get the haplotypes
    private VariantHaplotypeConverter converter;

    /** If override, this method should be called anyway for initialize the coverter. */
    @Override
    protected String[] customCommandLineValidation() {
        // get the haplotype model argument collection
        converter = haplotypeModelArgumentCollection.getHaplotypeConverter();
        return super.customCommandLineValidation();
    }

    /**
     * Implementation of variant-based traversal. Iterates over the variants, converting to
     * haplotypes using the parameters from {@link #haplotypeModelArgumentCollection} and filtering
     * afterwards by {@link #makeVariantFilter()}.
     *
     * Variants that pass the filters are passed to {@link #apply(VariantContext, ReadsContext,
     * ReferenceContext, FeatureContext)}
     */
    @Override
    public final void traverse() {
        converter.log(logger);
        final VariantFilter filter = makeVariantFilter();
        // Process each variant in the input stream.
        // ecause drivingVariants are private, we need to have the iterator here
        // TODO: contribute to GATK4 to get access to the drivingVariants
        final Iterable<VariantContext> iterator = ()
                -> features.getFeatureIterator(getDrivingVariantsFeatureInput());
        StreamSupport.stream(iterator.spliterator(), false)
                .map(converter) // converting to haplotypes
                .filter(filter)
                .forEach(variant -> {
                    final SimpleInterval variantInterval = new SimpleInterval(variant);
                    apply(variant,
                            new ReadsContext(reads, variantInterval),
                            new ReferenceContext(reference, variantInterval),
                            new FeatureContext(features, variantInterval));
                    progressMeter.update(variantInterval);
                });
    }
}

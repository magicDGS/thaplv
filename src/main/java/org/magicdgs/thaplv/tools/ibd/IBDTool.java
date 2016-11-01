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

package org.magicdgs.thaplv.tools.ibd;

import org.magicdgs.thaplv.cmd.argumentcollections.SlidingWindowArgumentCollection;
import org.magicdgs.thaplv.haplotypes.filters.NumberOfMissingFilter;
import org.magicdgs.thaplv.tools.ibd.engine.IBDOutput;
import org.magicdgs.thaplv.tools.ibd.engine.IBDcollector;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.HaploidWalker;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.VariantFilter;

/**
 * Abstract tool for IBD computation, both pair-wise or against the reference.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public abstract class IBDTool extends HaploidWalker {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output prefix for the divergence (.diff) and IBD tracks (.ibd)", optional = true)
    public String outputPrefix;

    @ArgumentCollection
    public SlidingWindowArgumentCollection slidingWindowArgumentCollection =
            new SlidingWindowArgumentCollection(500_000, 100_000);

    @Argument(fullName = "minimum-differences", shortName = "md", doc = "Threshold for number of differences per site in the window to output as an IBD", optional = true)
    public double minimumDifferences = 0.0005;

    @Argument(fullName = "output-differences", shortName = "out-diff", doc = "Output the divergence per window (.diff file) in addition to the IBD tracks")
    public boolean outputDiff = false;

    protected IBDOutput output;

    protected VCFHeader header;

    // TODO: remove this requirement to don't count Ns in the window
    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    protected boolean requiresOutputPloidy() {
        return false;
    }

    @Override
    protected boolean allowsCheckOnly() {
        return true;
    }

    /**
     * Make the variant filter to filter number of samples missing
     */
    protected VariantFilter makeVariantFilter() {
        logger.warn("Variants with only one sample genotyped will be skipped");
        return new NumberOfMissingFilter(header.getNGenotypeSamples() - 1);
    }

    /** Get the IBD collector for this tool. */
    protected abstract IBDcollector getIBDcollector();

    @Override
    protected String[] customCommandLineValidation() {
        // validate teh sliding window arguments
        slidingWindowArgumentCollection.validateArguments();
        return super.customCommandLineValidation();
    }

    @Override
    public void onTraversalStart() {
        header = getHeaderForVariants();
        output = new IBDOutput(outputPrefix, getIBDcollector(), minimumDifferences, !outputDiff);
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext,
            final ReferenceContext referenceContext, final FeatureContext featureContext) {
        output.addToCollector(variant);
    }

    /**
     * Close the output
     */
    @Override
    public void closeTool() {
        if (output != null) {
            output.close();
        }
    }

}

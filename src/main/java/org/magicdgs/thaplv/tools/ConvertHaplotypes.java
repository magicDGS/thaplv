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

package org.magicdgs.thaplv.tools;

import org.magicdgs.thaplv.cmd.ThaplvArgumentDefinitions;
import org.magicdgs.thaplv.cmd.programgroups.ConversionProgramGroup;
import org.magicdgs.thaplv.haplotypes.filters.HaplotypeFilterLibrary;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.HaploidWalker;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.VariantFilter;
import org.broadinstitute.hellbender.engine.filters.VariantFilterLibrary;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.io.File;

/**
 * Simple tool to convert haplotypes (either Inbred Lines or Back-crosses) to the haploid format
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
@CommandLineProgramProperties(
        summary = "Check diploid calls obtained with other callers and get the haploid/haplotypes "
                + "using the provided model. It allows to don't use any model afterwards (DONT_CHECK) "
                + "to speed up the analysis without conversion.",
        oneLineSummary = "Get the haplotypes from a backcross or inbred line experiment.",
        programGroup = ConversionProgramGroup.class)
public final class ConvertHaplotypes extends HaploidWalker {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "File to which variants should be written", common = false, optional = false)
    public File outFile;

    @Argument(fullName = ThaplvArgumentDefinitions.ONLY_POLYMORPHIC_LONG, shortName = ThaplvArgumentDefinitions.ONLY_POLYMORPHIC_SHORT, doc = "Only polymorphic sites will be output.", optional = true)
    public boolean onlyPolymorphic = false;

    @Override
    protected boolean requiresOutputPloidy() {
        return true;
    }

    @Override
    protected boolean allowsDontCheck() {
        return false;
    }

    private VariantContextWriter vcfWriter = null;

    /**
     * Make the variant filter to filter invariant sites
     */
    protected VariantFilter makeVariantFilter() {
        if (onlyPolymorphic) {
            return HaplotypeFilterLibrary.POLYMORPHIC_FILTER;
        }
        return VariantFilterLibrary.ALLOW_ALL_VARIANTS;
    }

    /**
     * Set up the VCF writer
     */
    @Override
    public void onTraversalStart() {
        // Initialize VCF header lines
        final VCFHeader header = getHeaderForVariants();
        header.addMetaDataLine(new VCFHeaderLine("source", this.getClass().getSimpleName()));
        vcfWriter = GATKVariantContextUtils
                .createVCFWriter(outFile, getBestAvailableSequenceDictionary(), false);
        vcfWriter.writeHeader(header);
    }

    @Override
    public void apply(VariantContext variant, ReadsContext readsContext,
            ReferenceContext referenceContext, FeatureContext featureContext) {
        // TODO: convert some of the annotations for the VCF for something
        vcfWriter.add(variant);
    }

    /**
     * Close out the new variants file.
     */
    @Override
    public void closeTool() {
        if (vcfWriter != null) {
            vcfWriter.close();
        }
    }
}

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

package org.magicdgs.thaplv.tools.imputation;

import org.magicdgs.thaplv.cmd.programgroups.AlphaProgramGroup;
import org.magicdgs.thaplv.tools.imputation.engine.ImputationCollector;
import org.magicdgs.thaplv.tools.imputation.engine.ImputationOutput;

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

import java.io.File;
import java.io.IOException;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
// TODO: add more details to the usage
@CommandLineProgramProperties(summary = "TODO: KNN method", oneLineSummary = "Impute a VCF file with the KNN algorithm.",
        programGroup = AlphaProgramGroup.class)
public class ImputeKNN extends HaploidWalker {

    @Argument(shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output imputed VCF file.")
    public File outputFile;

    @Argument(shortName = "WS",
            doc = "Number of base pair upstream/downstream a missing position to use in the imputation.", optional = true)
    public int windowSize = 1000;

    @Argument(shortName = "K", doc = "Number of neighbour to use in the KNN algorithm.", optional = true)
    public int kNeighbours = 5;
    //	@Argument(shortName = "stats", doc="Output imputation statistics", optional = true)
    //	public boolean IMPUTATION_STATISTICS = false;
    //	@Argument(shortName = "neigh", doc="Output the neighbours statistics", optional = true)
    //	public boolean NEIGHBOURS_STATISTICS = false;

    @Argument(shortName = "IMP", doc = "Add IMP format tag with the information for the imputed alleles", optional = true)
    public boolean imputeFormat = false;

    private ImputationOutput output;

    @Override
    protected boolean requiresOutputPloidy() {
        return true;
    }

    @Override
    protected boolean allowsDontCheck() {
        return true;
    }

    @Override
    public void onTraversalStart() {
        // initialize the writter
        final VCFHeader header = getHeaderForVariants();
        header.addMetaDataLine(new VCFHeaderLine("source", this.getClass().getSimpleName()));
        final VariantContextWriter writer =
                createVCFWriter(outputFile);
        // TODO: add command line to the header
        // get the imputation collector
        final ImputationCollector collector =
                new ImputationCollector(windowSize, kNeighbours, imputeFormat);
        output = new ImputationOutput(writer, collector);
        output.writeHeader(header, imputeFormat);
    }

    @Override
    public void apply(VariantContext variant, ReadsContext readsContext,
            ReferenceContext referenceContext, FeatureContext featureContext) {
        output.addVariant(variant);
    }

    @Override
    public void closeTool() {
        if (output != null) {
            try {
                output.close();
            } catch (IOException e) {
                // TODO: this should be changed
            }
        }
    }
}

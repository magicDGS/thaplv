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

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import org.vetmeduni.thaplv.cmdline.helpers.AbstractTool;
import org.vetmeduni.thaplv.tools.imputation.imputer.ImputationCollector;
import org.vetmeduni.thaplv.tools.imputation.imputer.ImputationOutput;
import org.vetmeduni.thaplv.utils.VariantLogger;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.Alpha;

import java.io.File;
import java.io.IOException;

/**
 * @author Daniel Gómez-Sánchez
 */
// TODO: add more details to the usage
@CommandLineProgramProperties(usage = "TODO: KNN method", usageShort = "Impute a VCF file with the KNN algorithm.",
        programGroup = Alpha.class)
public class ImputeKNN extends AbstractTool {

    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Output imputed VCF file.", maxElements = 1)
    public File OUTPUT;

    @Option(shortName = "WS",
            doc = "Number of base pair upstream/downstream a missing position to use in the imputation.", optional = true)
    public int WINDOW_SIZE = 1000;

    @Option(shortName = "K", doc = "Number of neighbour to use in the KNN algorithm.", optional = true)
    public int K_NEIGHBOURS = 5;
    //	@Option(shortName = "stats", doc="Output imputation statistics", optional = true)
    //	public boolean IMPUTATION_STATISTICS = false;
    //	@Option(shortName = "neigh", doc="Output the neighbours statistics", optional = true)
    //	public boolean NEIGHBOURS_STATISTICS = false;

    @Option(shortName = "IMP", doc = "Add IMP format tag with the information for the imputed alleles", optional = true)
    public boolean IMPUTE_FORMAT = false;

    @Override
    protected int doWork() {
        try {
            // initialize the reader
            final VCFFileReader reader = getHaplotypeReader();
            final VCFHeader header = new VCFHeader(reader.getFileHeader());
            final VariantContextWriter writer =
                    getVariantContextWriter(OUTPUT, header.getSequenceDictionary());
            // TODO: add command line to the header
            // get the imputation collector
            final ImputationCollector collector =
                    new ImputationCollector(WINDOW_SIZE, K_NEIGHBOURS, IMPUTE_FORMAT);
            final ImputationOutput output = new ImputationOutput(writer, collector);
            output.writeHeader(header, IMPUTE_FORMAT);
            final CloseableIterator<VariantContext> it = getIteratorForReader(reader);
            // initialize logger
            final VariantLogger progress = new VariantLogger(logger, 1000);
            // iterate
            while (it.hasNext()) {
                final VariantContext variant = it.next();
                output.addVariant(variant);
                progress.variant(variant);
            }
            progress.logNumberOfVariantsProcessed();
            CloserUtil.close(it);
            CloserUtil.close(reader);
            output.close();
            return 0;
        } catch (IOException e) {
            logger.error(e, e.getMessage());
            return 1;
        }
    }
}

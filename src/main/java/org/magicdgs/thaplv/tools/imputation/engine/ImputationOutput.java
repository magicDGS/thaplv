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

package org.magicdgs.thaplv.tools.imputation.engine;

import htsjdk.samtools.util.Log;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import java.io.Closeable;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Output for the imputation process
 *
 * @author Daniel Gómez-Sánchez
 */
public class ImputationOutput implements Closeable {

    private final static Log logger = Log.getInstance(ImputationOutput.class);

    private final VariantContextWriter output;

    private final ImputationCollector collector;

    // TODO: add info
    protected static VCFFormatHeaderLine knnFormatLine = new VCFFormatHeaderLine(
            ImputationCollector.KNN_FORMAT_ID, ImputationCollector.KNN_FORMAT_ORDER.length,
            VCFHeaderLineType.Integer,
            "Used input for imputation for this haplotype: " + Arrays
                    .toString(ImputationCollector.KNN_FORMAT_ORDER)
                    + ". TODO");

    protected static VCFInfoHeaderLine
            imputedInfoLine = new VCFInfoHeaderLine(ImputationCollector.IMPUTED_INFO_ID, 1,
            VCFHeaderLineType.Integer, "Number of imputed haplotypes");

    protected static VCFFormatHeaderLine distanceFormatLine =
            new VCFFormatHeaderLine(ImputationCollector.DISTANCE_FORMAT_ID,
                    VCFHeaderLineCount.R, VCFHeaderLineType.Float,
                    "K-Nearest Neighbours average distance for each allele");

    protected static VCFFormatHeaderLine frequencyFormatLine =
            new VCFFormatHeaderLine(ImputationCollector.FREQUENCY_FORMAT_ID,
                    VCFHeaderLineCount.R, VCFHeaderLineType.Integer,
                    "K-Nearest Neighbours counts for each allele");

    /**
     * Constructor
     *
     * @param writer    the write for the results (should be already initialized)
     * @param collector the collector to add the imputation result
     */
    public ImputationOutput(VariantContextWriter writer, final ImputationCollector collector) {
        this.output = writer;
        this.collector = collector;
    }

    /**
     * Write the header for the output, adding the information for the IMP format
     *
     * @param header the header to print
     */
    public void writeHeader(final VCFHeader header) {
        writeHeader(header, true);
    }

    /**
     * Write the header for the output
     *
     * @param header       the header to print
     * @param imputeFormat if <code>true</code>  imputed samples will have information
     */
    public void writeHeader(final VCFHeader header, boolean imputeFormat) {
        if (imputeFormat) {
            header.addMetaDataLine(distanceFormatLine);
            header.addMetaDataLine(knnFormatLine);
            header.addMetaDataLine(frequencyFormatLine);
        }
        header.addMetaDataLine(imputedInfoLine);
        this.output.writeHeader(header);
    }

    /**
     * Add a variant and output what is needed
     *
     * @param variant the variant to add
     */
    public void addVariant(VariantContext variant) {
        ArrayList<VariantContext> imputed = collector.addVariant(variant);
        for (VariantContext v : imputed) {
            output.add(v);
        }
    }

    @Override
    public void close() throws IOException {
        final ArrayList<VariantContext> last = new ArrayList<>();
        collector.finalizeImputation(last);
        for (VariantContext v : last) {
            output.add(v);
        }
        output.close();
    }
}

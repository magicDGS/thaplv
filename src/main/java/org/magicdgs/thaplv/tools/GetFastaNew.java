
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

package org.magicdgs.thaplv.tools;

import org.magicdgs.thaplv.cmd.argumentcollections.HaplotypeModelArgumentCollection;
import org.magicdgs.thaplv.cmd.argumentcollections.HaplotypeModelNoPloidyArgumentCollection;
import org.magicdgs.thaplv.cmd.programgroups.ConversionProgramGroup;
import org.magicdgs.thaplv.io.FastaWriter;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceWalker;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Simple tool for extract a FASTA file for each sample in a VCF substituting the variants in a
 * reference.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
@CommandLineProgramProperties(
        oneLineSummary = "Get a FASTA file for each haplotype in a VCF file.",
        summary =
                "Substitute in a provided FASTA file the variants in the VCF for each sample. Any missing haplotype will be output as an N. "
                        + "Provided intervals are regions of the VCF that will be used for update the reference, but the whole chromosome will be output.",
        programGroup = ConversionProgramGroup.class)
public final class GetFastaNew extends ReferenceWalker {

    @Argument(fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME, shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME, doc = "A VCF file containing variants", common = false, optional = false)
    public FeatureInput<VariantContext> drivingVariantFile;

    @ArgumentCollection(doc = "Haplotype models")
    public HaplotypeModelArgumentCollection haplotypeModelArgumentCollection =
            new HaplotypeModelNoPloidyArgumentCollection(true);

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output prefix before ${sampleName}.fasta", optional = true)
    public String outputPrefix = "";

    @Argument(fullName = "width", shortName = "w", doc = "Output width for sequence lines", optional = true)
    int sequenceWidth = 80;

    // maps sample and it's writer
    private HashMap<String, FastaWriter> writerMap;

    // maps sample and sequence for mutation)
    private HashMap<String, ReferenceSequence> sequenceMap;

    // stort the samples instead of passing them to the methods
    private Set<String> samples;

    private SimpleInterval currentSequence;

    @Override
    public void onTraversalStart() {
        final VCFHeader header = (VCFHeader) getHeaderForFeatures(drivingVariantFile);
        samples = new HashSet<>(header.getSampleNamesInOrder());
        writerMap = new HashMap<>(samples.size());
        for (final String s : samples) {
            final File file = new File(outputPrefix + s + ".fasta");
            IOUtil.assertFileIsWritable(file);
            try {
                writerMap.put(s, new FastaWriter(file, sequenceWidth));
            } catch (FileNotFoundException e) {
                throw new UserException.CouldNotCreateOutputFile(file, e.getMessage(), e);
            }
        }
        sequenceMap = new HashMap<>(samples.size());
    }

    @Override
    public void apply(SimpleInterval position, ReadsContext gatkReads,
            ReferenceContext referenceContext, FeatureContext featureContext) {
        // check if the current sequence is not initialize
        if (currentSequence == null) {
            logger.debug("Initializing current sequence");
            setNewCurrentSequence();
        }
        // check the chromosome
        if (!currentSequence.equals(getCurrentInterval())) {
            logger.debug("Change to the next current sequence");
            endSequences(currentSequence);
            setNewCurrentSequence();
        }
        // TODO: probably is better to filter for indels here
        // get the values for the variants
        final List<VariantContext> variants = featureContext.getValues(drivingVariantFile);
        if (variants.isEmpty()) {
            writeRefBaseForAll(referenceContext.getBase());
            return;
        } else if (variants.size() != 1) {
            logger.warn("Ommiting multiple variants located at {}", position);
            writeRefBaseForAll(referenceContext.getBase());
            return;
        }
        final VariantContext variant = variants.get(0);
        if (variant.isIndel()) {
            // log only if it is the first base that spans this location
            if (position.getStart() == variant.getStart()) {
                logger.warn("Indel variant ignored at position {}", position);
            }
            // TODO: is this correct? probably indels should write the reference position
            // TODO: if so, this will be in only inside the previous if...
            return;
        }
        for (final Genotype genotype : variant.getGenotypes()) {
            final FastaWriter sampleWriter = writerMap.get(genotype.getSampleName());
            // TODO: putative difference with previous regions, assuming haploids directly
            final Allele al = genotype.getAllele(0);
            if (al.isNoCall()) {
                sampleWriter.writeBase((byte) 'N');
            } else {
                sampleWriter.writeBase(al.getBases()[0]);
            }
        }
    }

    private void setNewCurrentSequence() {
        currentSequence = new SimpleInterval(getCurrentInterval());
        final String sequenceName =
                (isWholeChromosome(currentSequence)) ? currentSequence.getContig()
                        : currentSequence.toString();
        logger.info("Starting {}.", sequenceName);
        for (final FastaWriter writer : writerMap.values()) {
            writer.writeHeader(sequenceName);
        }
    }

    private boolean isWholeChromosome(final SimpleInterval interval) {
        final SAMSequenceRecord record = getReferenceDictionary().getSequence(interval.getContig());
        final SimpleInterval wholeRefInterval =
                new SimpleInterval(record.getSequenceName(), 1, record.getSequenceLength());
        return interval.equals(wholeRefInterval);
    }

    private void endSequences(final SimpleInterval interval) {
        logger.info("Ending {}", interval);
        for (final Map.Entry<String, FastaWriter> entry : writerMap.entrySet()) {
            logger.debug("Ending for {}", entry.getKey());
            entry.getValue().endSequence();
        }
    }

    private void writeRefBaseForAll(final byte ref) {
        for (final FastaWriter writer : writerMap.values()) {
            writer.writeBase(ref);
        }
    }

    @Override
    public void closeTool() {
        for (final FastaWriter w : writerMap.values()) {
            CloserUtil.close(w);
        }
    }
}
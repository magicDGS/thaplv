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

import org.magicdgs.thaplv.cmd.programgroups.ConversionProgramGroup;
import org.magicdgs.thaplv.io.FastaWriter;

import htsjdk.samtools.SAMException;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.HaploidWalker;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * Simple tool for extract a FASTA file for each sample in a VCF substituting the variants in a
 * reference
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
@CommandLineProgramProperties(
        oneLineSummary = "Get a FASTA file for each haplotype in a VCF file.",
        summary =
                "Substitute in a provided FASTA file the variants in the VCF for each sample. Any missing haplotype will be output as an N. "
                        + "Provided intervals are regions of the VCF that will be used for update the reference, but the whole chromosome will be output "
                        + "(regions outside the provided intervals will be the same as the reference).",
        programGroup = ConversionProgramGroup.class)
public final class GetFasta extends HaploidWalker {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output prefix before sample specification. If a folder, it should end with '/'", optional = true)
    public String outputPrefix = "";

    @Argument(fullName = "width", shortName = "w", doc = "Output width for sequence lines", optional = true)
    public int sequenceWidth = 80;

    // maps sample and it's writer
    private Map<String, FastaWriter> writerMap;

    // maps sample and sequence for mutation)
    private Map<String, ReferenceSequence> sequenceMap;

    // stort the samples instead of passing them to the methods
    private Set<String> samples;

    // the reference file
    private IndexedFastaSequenceFile reference;

    // last visited contig
    private String lastContig = null;

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

    @Override
    public void onTraversalStart() {
        try {
            reference = new IndexedFastaSequenceFile(referenceArguments.getReferenceFile());
        } catch (FileNotFoundException e) {
            throw new UserException.MissingReferenceFaiFile(
                    new File(referenceArguments.getReferenceFileName() + ".fai"),
                    referenceArguments.getReferenceFile());
        }
        samples = new HashSet<>(getHeaderForVariants().getSampleNamesInOrder());
        writerMap = new HashMap<>(samples.size());
        for (final String s : samples) {
            final File file = new File(outputPrefix + s + ".fasta");
            try {
                IOUtil.assertFileIsWritable(file);
                writerMap.put(s, new FastaWriter(file, sequenceWidth));
            } catch (FileNotFoundException | SAMException e) {
                throw new UserException.CouldNotCreateOutputFile(file, e.getMessage(), e);
            }
        }
        sequenceMap = new HashMap<>(samples.size());
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext,
            final ReferenceContext referenceContext, final FeatureContext featureContext) {
        final int startPos = variant.getStart();
        // ignoring indels
        if (variant.isIndel()) {
            logger.warn("Indel variant ignored at position {}:{}", variant.getContig(), startPos);
            return;
        }
        lastContig = updateSequenceMap(variant.getContig(), lastContig, reference);
        for (final Genotype genotype : variant.getGenotypes(samples)) {
            // TODO: putative difference with previous regions, assuming haploids directly
            final Allele al = genotype.getAllele(0);
            if (al.isNoCall()) {
                mutateSequenceForSample(genotype.getSampleName(), startPos, (byte) 'N');
            } else {
                mutateSequenceForSample(genotype.getSampleName(), startPos, al.getBases()[0]);
            }
        }
    }

    @Override
    public Object onTraversalSuccess() {
        outputSequences(lastContig);
        return null;
    }

    /**
     * Mutate the sequence for this sample witht he base provided
     *
     * @param sampleName the sample to mutate
     * @param position   1-based position
     * @param base       the new base in the position
     */
    private void mutateSequenceForSample(final String sampleName, final int position,
            final byte base) {
        sequenceMap.get(sampleName).getBases()[position - 1] = base;
    }

    /**
     * Update the sequence map and output sequences if possible
     *
     * @param currentContig the current contig
     * @param lastContig    the last contig
     * @param sequence      the sequence file
     *
     * @return the contig in current use
     */
    private String updateSequenceMap(final String currentContig, final String lastContig,
            final IndexedFastaSequenceFile sequence) {
        if (currentContig.equals(lastContig)) {
            return lastContig;
        }
        if (lastContig != null) {
            outputSequences(lastContig);
        } else {
            sequenceMap = new HashMap<>(samples.size());
        }
        // speed-up analysis making a sequence template (only needs to read the file once)
        final ReferenceSequence template = sequence.getSequence(currentContig);
        for (final String s : samples) {
            sequenceMap.put(s, new ReferenceSequence(currentContig, template.getContigIndex(),
                    Arrays.copyOf(template.getBases(), template.length())));
        }
        return currentContig;
    }

    /**
     * Write down all the sequences stored in the map and clear it
     */
    private void outputSequences(final String contig) {
        // TODO: we should only output the requested intervals sequence
        // intervalArgumentCollection.getIntervals(getReferenceDictionary());
        logger.info("Writting sequence {} for each sample.", contig);
        for (final Map.Entry<String, ReferenceSequence> entry : sequenceMap.entrySet()) {
            writerMap.get(entry.getKey()).addWholeSequence(entry.getValue());
            logger.debug("Finished " + entry.getKey());
        }
        sequenceMap.clear();
    }

    @Override
    public void closeTool() {
        CloserUtil.close(reference);
        for (final FastaWriter w : writerMap.values()) {
            CloserUtil.close(w);
        }
    }
}
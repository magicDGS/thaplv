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
package org.magicdgs.thaplv.io;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SequenceUtil;

import java.io.Closeable;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

/**
 * Class for count the number of Ns in the reference sequence provided
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class FastaNsCounter implements Closeable {

    private final IndexedFastaSequenceFile fastaReader;

    /**
     * Default constructor
     *
     * @param fasta the fasta file
     */
    public FastaNsCounter(final File fasta) throws FileNotFoundException {
        fastaReader = new IndexedFastaSequenceFile(fasta);
    }

    /**
     * Count the Ns in this position
     *
     * @param region the coordinates to count the Ns
     *
     * @return the number of N nucleotides in the region
     */
    public int countNsRegion(final Locatable region) {
        int count = 0;
        for (final byte base : fastaReader
                .getSubsequenceAt(region.getContig(), region.getStart(), region.getEnd())
                .getBases()) {
            if (SequenceUtil.isNoCall(base)) {
                count++;
            }
        }
        return count;
    }

    /**
     * Get the chromosome length for this dictionary
     *
     * @param chromosome the chromosome to check
     *
     * @return the length of the chromosome
     */
    public int getChromosomeLength(final String chromosome) {
        if (getDictionary() == null) {
            throw new IllegalStateException("Reference FASTA file needs a sequence dictionary.");
        }
        return getDictionary().getSequence(chromosome).getSequenceLength();
    }

    /**
     * Get the sequence dictionary for this fasta file
     *
     * @return the sequence dictionary
     */
    public SAMSequenceDictionary getDictionary() {
        return fastaReader.getSequenceDictionary();
    }

    @Override
    public void close() throws IOException {
        fastaReader.close();
    }
}

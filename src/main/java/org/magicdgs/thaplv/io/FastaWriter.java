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

import htsjdk.samtools.reference.ReferenceSequence;

import java.io.Closeable;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;

/**
 * Simple class to write a FASTA file with whole sequences
 *
 * @author Daniel Gómez-Sánchez
 */
public class FastaWriter implements Closeable {

    private final PrintStream writer;

    private final int lineWidth;

    private boolean headerAlreadyWritten;

    private int counter;

    /**
     * Constructor with a file
     *
     * @param file      the file to output
     * @param lineWidth the width for the output
     */
    public FastaWriter(final File file, final int lineWidth) throws FileNotFoundException {
        this(new PrintStream(file), lineWidth);
    }

    /**
     * Constructor with a PrintStream
     *
     * @param writer    the writer to use
     * @param lineWidth the width for the output
     */
    public FastaWriter(final PrintStream writer, final int lineWidth) {
        this.writer = writer;
        this.lineWidth = lineWidth;
        this.headerAlreadyWritten = false;
        this.counter = 0;
    }

    /**
     * Add a sequence to print out
     *
     * @param sequence the sequence to print out
     */
    public void addSequence(final ReferenceSequence sequence) {
        writeHeader(sequence.getName());
        int counter = 0;
        for (byte base : sequence.getBases()) {
            writeBase(base);
        }
        endSequence();
    }

    public void writeHeader(final String header) {
        if(headerAlreadyWritten) {
            throw new IllegalStateException("Header already written");
        }
        writer.print(">");
        writer.println(header);
        headerAlreadyWritten = true;
    }

    public void writeBase(final byte base) {
        writer.print((char) base);
        if(++counter == lineWidth) {
            writer.println();
            // TODO: flush for low memory?
            // writer.flush();
            counter = 0;
        }
    }

    public void endSequence() {
        if(counter != 0) {
            writer.println();
        }
        writer.flush();
        counter = 0;
        headerAlreadyWritten = false;
    }

    @Override
    public void close() throws IOException {
        endSequence();
        writer.close();
    }
}

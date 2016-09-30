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

package org.magicdgs.thaplv.tools.ld.engine;

import org.magicdgs.thaplv.utils.stats.LengthBinning;

import java.io.BufferedWriter;
import java.io.Closeable;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class BinLDwriter implements Closeable {

    private BufferedWriter[] writers;
    // format for the values
    private static DecimalFormat roundedFormat = new DecimalFormat("#.#######");

    public BinLDwriter(String outputPrefix, LDbinning emptyBins) throws IOException {
        writers = new BufferedWriter[LDbinning.statsToBin.size()];
        for (int i = 0; i < LDbinning.statsToBin.size(); i++) {
            writers[i] = new BufferedWriter(new FileWriter(
                    String.format("%s.%s", outputPrefix, LDbinning.statsToBin.get(i))));
            writers[i].write("chr\t");
            writers[i].write(emptyBins.getBin(i).formatHeaderBins());
            writers[i].newLine();
        }
    }

    public void write(String contig, LDbinning bins) throws IOException {
        for (int i = 0; i < LDbinning.statsToBin.size(); i++) {
            // get the length binning
            LengthBinning lb = bins.getBin(i);
            writers[i].write(contig);
            writers[i].write('\t');
            writers[i].write(lb.formatBins(roundedFormat));
        }

    }

    @Override
    public void close() throws IOException {
        for (BufferedWriter w : writers) {
            w.close();
        }
    }
}

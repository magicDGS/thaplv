/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2015 Daniel Gómez-Sánchez
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

package org.magicdgs.thaplv.tools.ld.engine;

import org.magicdgs.thaplv.utils.stats.LengthBinning;

import java.io.BufferedWriter;
import java.io.Closeable;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;

/**
 * @author Daniel Gómez-Sánchez
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

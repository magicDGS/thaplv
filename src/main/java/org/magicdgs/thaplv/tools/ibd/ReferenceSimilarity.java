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

package org.magicdgs.thaplv.tools.ibd;

import org.magicdgs.thaplv.cmd.programgroups.AnalysisProgramGroup;
import org.magicdgs.thaplv.io.FastaNsCounter;
import org.magicdgs.thaplv.tools.ibd.engine.IBDcollector;

import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.File;
import java.io.FileNotFoundException;

/**
 * Tool for compute IBD regions against a reference.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
@CommandLineProgramProperties(
        summary =
                "Identify identity-by-descent (IBD) regions for the reference using the approach described in Langley et. al 2012, Genetics 192(2). "
                        + "Differences between individuals and the reference (from the VCF reference allele) are computed in windows, and called as IBD tracks if the "
                        + "differences per site are smaller than a threshold.",
        oneLineSummary = "Compute reference similarity (IBD regions w.r.t. reference) using the number of differences.",
        programGroup = AnalysisProgramGroup.class)
public final class ReferenceSimilarity extends IBDTool {

    @Override
    protected IBDcollector getIBDcollector() {
        final FastaNsCounter nCounter;
        try {
            nCounter = new FastaNsCounter(referenceArguments.getReferenceFile());
        } catch (FileNotFoundException e) {
            throw new UserException.MissingReferenceFaiFile(
                    new File(referenceArguments.getReferenceFileName() + ".fai"),
                    referenceArguments.getReferenceFile());
        }
        return new IBDcollector(header.getSampleNamesInOrder(), nCounter,
                slidingWindowArgumentCollection.windowSize,
                slidingWindowArgumentCollection.stepSize);
    }
}

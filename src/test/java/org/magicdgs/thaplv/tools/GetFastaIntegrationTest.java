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


import org.magicdgs.thaplv.cmd.ThaplvArgumentDefinitions;
import org.magicdgs.thaplv.utils.test.BaseTest;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class GetFastaIntegrationTest extends CommandLineProgramTest {

    @Test
    public void testGetFasta() throws IOException {
        final String[] samples = {"sample1", "sample2"};
        final File tmpDir = BaseTest.createTempDir("getFasta");
        final ArgumentsBuilder arguments = new ArgumentsBuilder();
        // add the input
        arguments.addVCF(new File(TEST_ROOT_FILE_DIRECTORY
                + this.getClass().getPackage().getName().replace(".", "/")
                + "/ConvertHaplotypes/example.vcf"));
        // add the reference
        arguments.addReference(DROSOPHILA_SIMULANS_2L_REFERENCE);
        // add the model
        arguments.addArgument(ThaplvArgumentDefinitions.HAPLOTYPE_MODEL_LONG, "HAPLOID");
        // the output
        arguments.addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, tmpDir.getAbsolutePath() + "/");
        // run command line
        // TODO: this test takes a long time, because the current implementation does not allow
        // TODO: output only a region. This implementation should be changed
        runCommandLine(arguments);
        for (final String sampl : samples) {
            final File sampleFile = new File(tmpDir, sampl + ".fasta");
            Assert.assertTrue(sampleFile.exists(), "output file does not exists: " + sampleFile.toString());
            IntegrationTestSpec.assertEqualTextFiles(
                    sampleFile, getTestFile("expected_haploid." + sampl + ".fasta"));
        }
    }

}
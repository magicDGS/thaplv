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

import org.magicdgs.thaplv.cmd.ThaplvArgumentDefinitions;
import org.magicdgs.thaplv.utils.test.BaseTest;
import org.magicdgs.thaplv.utils.test.CommandLineProgramTest;

import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class IBDcomputeIntegrationTest extends CommandLineProgramTest {

    private static final String expectedPrefix = "expected";

    /** Base arguments for the reference and the haploid model */
    static ArgumentsBuilder getReferenceAndHaplotypeModelForIBDCompute() {
        return new ArgumentsBuilder().addReference(DROSOPHILA_SIMULANS_2L_REFERENCE)
                .addArgument(ThaplvArgumentDefinitions.HAPLOTYPE_MODEL_LONG, "HAPLOID");
    }

    /** Base arguments includes the example.vcf, the reference and the haploid model */
    private ArgumentsBuilder getBaseArgumentsForIBDCompute() {
        return getReferenceAndHaplotypeModelForIBDCompute().addVCF(getTestFile("example.vcf.gz"));
    }

    @Test
    public void testIBDcompute() throws Exception {
        final File tmpDir = BaseTest.createTempDir("ibdCompute");
        final ArgumentsBuilder arguments = getBaseArgumentsForIBDCompute();
        arguments.addBooleanArgument("output-differences", true);
        // the output
        final String outputPrefix = "testIBDcompute";
        arguments.addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME,
                tmpDir.getAbsolutePath() + "/" + outputPrefix);
        // run command line
        runCommandLine(arguments);
        for (final String extension : new String[] {".ibd", ".diff"}) {
            final File outputFile = new File(tmpDir, outputPrefix + extension);
            Assert.assertTrue(outputFile.exists(),
                    "output file does not exists: " + outputFile.toString());
            IntegrationTestSpec
                    .assertEqualTextFiles(outputFile, getTestFile(expectedPrefix + extension));
        }
    }

    @Test
    public void testIBDcomputeCombiningWindows() throws Exception {
        final File tmpDir = BaseTest.createTempDir("ibdCompute");
        final ArgumentsBuilder arguments = getBaseArgumentsForIBDCompute();
        // the window-size
        arguments.addArgument(ThaplvArgumentDefinitions.WINDOW_SIZE_LONG, "200")
                .addArgument(ThaplvArgumentDefinitions.WINDOW_STEP_LONG, "100");
        // the output
        final String outputPrefix = "testIBDcompute_200_100";
        arguments.addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME,
                tmpDir.getAbsolutePath() + "/" + outputPrefix);
        // run command line
        runCommandLine(arguments);
        File outputFile = new File(tmpDir, outputPrefix + ".diff");
        // .diff does not exists
        Assert.assertFalse(outputFile.exists(), "diff file does not exists");
        // .ibd exists and the file is the same
        outputFile = new File(tmpDir, outputPrefix + ".ibd");
        Assert.assertTrue(outputFile.exists(), "ibd file does not exists");
        IntegrationTestSpec
                .assertEqualTextFiles(outputFile,
                        getTestFile(expectedPrefix + "_200_100.ibd"));
    }

    @Test
    public void testNotOutputDiff() throws Exception {
        final File tmpDir = BaseTest.createTempDir("ibdCompute");
        final ArgumentsBuilder arguments = getBaseArgumentsForIBDCompute();
        // the output
        final String outputPrefix = "testNotOutputDiff";
        arguments.addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME,
                tmpDir.getAbsolutePath() + "/" + outputPrefix);
        // run command line
        runCommandLine(arguments);
        File outputFile = new File(tmpDir, outputPrefix + ".diff");
        // .diff does not exists
        Assert.assertFalse(outputFile.exists(), "diff file exists");
        // .ibd exists and the file is the same
        outputFile = new File(tmpDir, outputPrefix + ".ibd");
        Assert.assertTrue(outputFile.exists(), "ibd file does not exists");
        IntegrationTestSpec
                .assertEqualTextFiles(outputFile, getTestFile(expectedPrefix + ".ibd"));
    }

    @Test(expectedExceptions = UserException.BadArgumentValue.class)
    public void testBadWindowSize() throws Exception {
        final ArgumentsBuilder arguments = getBaseArgumentsForIBDCompute()
                .addArgument(ThaplvArgumentDefinitions.WINDOW_SIZE_LONG, "0")
                .addOutput(createTempFile("testBadWindowSize", ""));
        runCommandLine(arguments);
    }

    @Test(expectedExceptions = UserException.BadArgumentValue.class)
    public void testBadStepSize() throws Exception {
        final ArgumentsBuilder arguments = getBaseArgumentsForIBDCompute()
                .addArgument(ThaplvArgumentDefinitions.WINDOW_STEP_LONG, "0")
                .addOutput(createTempFile("testBadStepSize", ""));
        runCommandLine(arguments);
    }

    @Test(expectedExceptions = UserException.BadArgumentValue.class)
    public void testBadStepSizeWithWindowSize() throws Exception {
        final ArgumentsBuilder arguments = getBaseArgumentsForIBDCompute()
                .addArgument(ThaplvArgumentDefinitions.WINDOW_SIZE_LONG, "100")
                .addArgument(ThaplvArgumentDefinitions.WINDOW_STEP_LONG, "200")
                .addOutput(createTempFile("testBadStepSizeWithWindowSize", ""));
        runCommandLine(arguments);
    }

}
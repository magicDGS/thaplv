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

package org.magicdgs.thaplv.tools.ld;


import org.magicdgs.thaplv.cmd.ThaplvArgumentDefinitions;
import org.magicdgs.thaplv.tools.ld.engine.LDdecayOutput;
import org.magicdgs.thaplv.utils.test.CommandLineProgramTest;

import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.text.TextFormattingUtils;
import org.broadinstitute.hellbender.utils.text.XReadLines;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.List;
import java.util.stream.Collectors;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class LDdecayIntegrationTest extends CommandLineProgramTest {

    /** This is the delta for the heuristic values. */
    private final static double APPROXIMATION_DELTA = 0.001;

    /** This is the maximum percentage of heuristic values that are allowed to fail. */
    private final static double MAXIMUM_PERCENTAGE_HEURISTICS = 10;

    private final static File vcfInput = getCommonTestFile("10samples.vcf.gz");

    private final static File tmpDir = createTempDir("ldDefault");

    private static final ArgumentsBuilder getBaseArgumentsBuilder() {
        return new ArgumentsBuilder()
                .addArgument(ThaplvArgumentDefinitions.HAPLOTYPE_MODEL_LONG, "HAPLOID")
                .addVCF(vcfInput);
    }


    @DataProvider
    public Object[][] badArgs() {
        return new Object[][] {
                {getBaseArgumentsBuilder().addArgument(LDdecay.CHI_SQR_QUANTILE_ARGNAME, "1")},
                {getBaseArgumentsBuilder().addArgument(LDdecay.CHI_SQR_QUANTILE_ARGNAME, "-1")},
                {getBaseArgumentsBuilder().addArgument(LDdecay.CHI_SQR_QUANTILE_ARGNAME, "100")},
                {getBaseArgumentsBuilder().addArgument(LDdecay.MINIMUM_SAMPLES_ARGNAME, "-1")}
        };
    }

    @Test(dataProvider = "badArgs", expectedExceptions = UserException.BadArgumentValue.class)
    public void testBadArguments(final ArgumentsBuilder args) throws Exception {
        args.addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, "badArguments");
        runCommandLine(args);
    }

    /** The expected files where generated with single thread versions. */
    @DataProvider(name = "singleThreadTests")
    public Object[][] getSingleThreadDataForTests() {
        return new Object[][] {
                // test including singletons
                {"testIncludeSingletons", "expected_IncludeSingletons", 384824,
                        getBaseArgumentsBuilder()
                                .addBooleanArgument(LDdecay.INCLUDE_SINGLETONS_ARGNAME, true)},
                // default arguments
                {"testDefaultArguments", "expected_default", 49169, getBaseArgumentsBuilder()},
                // test the number of missing argument
                {"testMinimumMissing5", "expected_missing5", 107461, getBaseArgumentsBuilder()
                        .addArgument(LDdecay.MINIMUM_SAMPLES_ARGNAME, "5")},
                // test the maximum-minimum distance ranges
                {"testMaxMinRange", "expected_MaxMinRange", 5307, getBaseArgumentsBuilder()
                        .addArgument("minimum-distance", "100")
                        .addArgument("maximum-distance", "1000")}

        };
    }

    /** Note: The arguments should not contain the output. */
    @Test(dataProvider = "singleThreadTests")
    public void testSingleThread(final String testName, final String expectedPrefix,
            final int expectedPairs, final ArgumentsBuilder args) throws Exception {
        log("Running test: " + testName);
        args.addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME,
                tmpDir.getAbsolutePath() + "/" + testName); // output
        // run command line
        final int computedPairs = (int) runCommandLine(args);
        Assert.assertEquals(computedPairs, expectedPairs, "wrong number of computed pairs");
        // testing equal files
        for (final String suffix : LDdecayOutput.statsToBin) {
            final File expected = getTestFile(expectedPrefix + "." + suffix);
            final File result = new File(tmpDir, testName + "." + suffix);
            Assert.assertTrue(result.exists());
            IntegrationTestSpec.assertEqualTextFiles(result, expected);
        }
    }

    @DataProvider(name = "multiThreadTests")
    public Object[][] getMultiThreadDataForTests() {
        return new Object[][] {
                // test including singletons
                {"testIncludeSingletons", "expected_IncludeSingletons", 384824,
                        getBaseArgumentsBuilder()
                                .addBooleanArgument(LDdecay.INCLUDE_SINGLETONS_ARGNAME, true)},
                // default arguments
                {"testDefaultArguments", "expected_default", 49169, getBaseArgumentsBuilder()},
                // test the number of missing argument
                {"testMinimumMissing5", "expected_missing5", 107461, getBaseArgumentsBuilder()
                        .addArgument(LDdecay.MINIMUM_SAMPLES_ARGNAME, "5")},
                // TODO: this test only have 5 heuristic values
                // test the maximum-minimum distance ranges
                // {"testMaxMinRange", "expected_MaxMinRange", 5307, getBaseArgumentsBuilder()
                //        .addArgument("minimum-distance", "100")
                //        .addArgument("maximum-distance", "1000")}
        };
    }

    /** Note: The arguments should not contain the output. */
    @Test(dataProvider = "multiThreadTests")
    public void testMultiThread(final String testName, final String expectedPrefix,
            final int expectedPairs, final ArgumentsBuilder args) throws Exception {
        final String testPrefix = testName + ".multi-thread";
        log("Running multi-thread test: " + testName);
        args.addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME,
                tmpDir.getAbsolutePath() + "/" + testPrefix) // output
                .addArgument("threads", "2");
        // run command line
        final int computedPairs = (int) runCommandLine(args);
        Assert.assertEquals(computedPairs, expectedPairs, "wrong number of computed pairs");
        // testing all files are created
        for (final String suffix : LDdecayOutput.statsToBin) {
            Assert.assertTrue(new File(tmpDir, testPrefix + "." + suffix).exists());
        }
        // testing only the file with r2, because it have more data points
        // this means that the multi-thread version is working properly
        final File result = new File(tmpDir, testPrefix + ".r2");
        final File expected = getTestFile(expectedPrefix + ".r2");
        assertEqualOutputFromApproximation(result, expected);
    }

    /**
     * This is used for test the output with approximate heuristic values derived from multi-thread
     * processing. First, we compare doubles with a less stringent delta
     * ({@link #APPROXIMATION_DELTA}) and count how many times this fail; if it fails more than
     * {@link #MAXIMUM_PERCENTAGE_HEURISTICS} of the tests in the file, it fails.
     */
    private static void assertEqualOutputFromApproximation(final File resultFile,
            final File expectedFile) throws Exception {

        // read all the lines
        final List<String[]> actualLines = new XReadLines(resultFile).readLines().stream()
                .map(TextFormattingUtils::splitWhiteSpace).collect(Collectors.toList());
        final List<String[]> expectedLines = new XReadLines(expectedFile).readLines().stream()
                .map(TextFormattingUtils::splitWhiteSpace).collect(Collectors.toList());

        final String[] actualHeader = actualLines.remove(0);
        final String[] expectedHeader = expectedLines.remove(0);
        Assert.assertEquals(actualHeader, expectedHeader, "not equal headers");
        //For ease of debugging, we look at the lines first and only then check their counts
        final int minLen = Math.min(actualLines.size(), expectedLines.size());
        // heuristic values counts
        double failedHeuristic = 0;
        double totalHeuristic = 0;
        for (int i = 0; i < minLen; i++) {
            final String[] actual = actualLines.get(i);
            final String[] expected = expectedLines.get(i);
            // non-heuristic values
            Assert.assertEquals(actual[0], expected[0],
                    "not equal" + actualHeader[0] + "at line " + i + 1);
            Assert.assertEquals(Integer.valueOf(actual[1]), Integer.valueOf(expected[1]),
                    "not equal " + actualHeader[1] + " at line " + i + 1);
            Assert.assertEquals(Integer.valueOf(actual[2]), Integer.valueOf(expected[2]),
                    "not equal " + actualHeader[2] + " at line " + i + 1);
            Assert.assertEquals(Double.valueOf(actual[3]), Double.valueOf(expected[3]),
                    DEFAULT_TOLERANCE, "not equal " + actualHeader[3] + " value at line " + i + 1);
            Assert.assertEquals(Double.valueOf(actual[4]), Double.valueOf(expected[4]),
                    DEFAULT_TOLERANCE, "not equal " + actualHeader[4] + " value at line " + i + 1);
            for (int j = 5; j < actual.length; j++) {
                try {
                    totalHeuristic++;
                    Assert.assertEquals(Double.valueOf(actual[j]), Double.valueOf(expected[j]),
                            APPROXIMATION_DELTA,
                            "not equal heuristic " + actualHeader[j] + " at line " + i + 1);
                } catch (AssertionError e) {
                    failedHeuristic++;
                }
            }
        }
        Assert.assertEquals(actualLines.size(), expectedLines.size(), "line counts");
        // heuristics
        final double heuristicPercentage = 100 * failedHeuristic / totalHeuristic;
        log(String.format("%s%% failed heuristics: %s/%s",
                heuristicPercentage, failedHeuristic, totalHeuristic));
        Assert.assertTrue(heuristicPercentage < MAXIMUM_PERCENTAGE_HEURISTICS);
    }

}
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

package org.magicdgs.thaplv.utils.test;

import htsjdk.samtools.util.Log;
import htsjdk.variant.variantcontext.Genotype;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.LoggingUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.Reporter;
import org.testng.annotations.BeforeSuite;

import java.io.File;

/**
 * Base test class for thaplv tests. All test cases should extend from this class;
 * it sets up the logger, and resolves the location of directories that we rely on.
 *
 * This class have a lot of code from {@link org.broadinstitute.hellbender.utils.test.BaseTest} and
 * also contains several static methods useful for testing.
 *
 * Test should be named as following:
 *
 * 1) ClassNameUnitTest for unit tests of single classes
 * 2) ClassNameIntegrationTest for tool tests
 * 3) ClassNameTest for other tests
 *
 * In addition, any resource file specific for a class should be sited under
 * {@link #TEST_ROOT_FILE_DIRECTORY}/package_path/ClassName/testFile
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public abstract class BaseTest {

    @BeforeSuite
    public void setTestVerbosity() {
        LoggingUtils.setLoggingLevel(Log.LogLevel.DEBUG);
    }

    public static final Logger logger = LogManager.getLogger("org.magicdgs.thaplv");

    private static final String CURRENT_DIRECTORY = System.getProperty("user.dir");

    // directory for all test files
    public static final String TEST_ROOT_FILE_DIRECTORY =
            new File(CURRENT_DIRECTORY, "src/test/resources/").getAbsolutePath() + "/";

    /** Log this message so that it shows up inline during output as well as in html reports */
    public static void log(final String message) {
        Reporter.log(message, true);
    }

    /** Default tolerance for float and double tests */
    public static final double DEFAULT_TOLERANCE = 1e-6; // quite stringent

    /** Returns the location of the resource directory for the tested class */
    public String getTestDataDir() {
        return TEST_ROOT_FILE_DIRECTORY + "/"
                + this.getClass().getPackage().getName().replace(".", "/")
                + "/" + getTestedClassName() + "/";
    }

    /** Reference genome from the D. simulans chromosome arm 2L */
    public static final File DROSOPHILA_SIMULANS_2L_REFERENCE =
            new File(TEST_ROOT_FILE_DIRECTORY + "org/magicdgs/thaplv/drosophila.2L.fa");

    /**
     * Get a test file resolved using {@link #getTestDataDir()} as the parent
     *
     * @param fileName the name of a file
     *
     * @return the test file
     */
    public File getTestFile(String fileName) {
        return new File(getTestDataDir(), fileName);
    }

    /**
     * Creates a temp file that will be deleted on exit after tests are complete.
     *
     * This will also mark the corresponding Tribble/Tabix/BAM indices matching the temp file for
     * deletion.
     *
     * @param name      Prefix of the file.
     * @param extension Extension to concat to the end of the file.
     *
     * @return A file in the temporary directory starting with name, ending with extension, which
     * will be deleted after the program exits.
     */
    public static File createTempFile(final String name, final String extension) {
        return IOUtils.createTempFile(name, extension);
    }

    /**
     * Creates an empty temp directory which will be deleted on exit after tests are complete
     *
     * @param prefix prefix for the directory name
     *
     * @return an empty directory starting with prefix which will be deleted after the program exits
     */
    public static File createTempDir(final String prefix) {
        final File dir = IOUtils.tempDir(prefix, "");
        IOUtils.deleteRecursivelyOnExit(dir);
        return dir;
    }

    /** Returns the name of the class being tested */
    public final String getTestedClassName() {
        if (getClass().getSimpleName().contains("IntegrationTest")) {
            return getClass().getSimpleName().replaceAll("IntegrationTest$", "");
        } else if (getClass().getSimpleName().contains("UnitTest")) {
            return getClass().getSimpleName().replaceAll("UnitTest$", "");
        } else {
            return getClass().getSimpleName().replaceAll("Test$", "");
        }
    }

    /**
     * Assert that two genotypes are equals
     *
     * @param actual      the actual genotype
     * @param expected    the expected genotype
     * @param ignorePhase if {@code true} phasing will be ignored when comparing haplotypes
     * @param msg         message to pass if the test fails
     */
    public static void assertGenotypesEquals(final Genotype actual, final Genotype expected,
            final boolean ignorePhase, final String msg) {
        if (!expected.sameGenotype(actual, ignorePhase)) {
            String formatted = "";
            if (null != msg) {
                formatted = msg + " ";
            }
            formatted = formatted + "expected [" + expected + "] but found [" + actual + "]";
            Assert.fail(formatted);
        }
    }

    /**
     * Assert that two genotypes are equals
     *
     * @param actual      the actual genotype
     * @param expected    the expected genotype
     * @param ignorePhase if {@code true} phasing will be ignored when comparing haplotypes
     */
    public static void assertGenotypesEquals(final Genotype actual, final Genotype expected,
            final boolean ignorePhase) {
        assertGenotypesEquals(actual, expected, ignorePhase, null);
    }

    /**
     * Assert that two genotypes are equals, ignoring the phase
     *
     * @param actual   the actual genotype
     * @param expected the expected genotype
     */
    public static void assertGenotypesEquals(final Genotype actual, final Genotype expected) {
        assertGenotypesEquals(actual, expected, false);
    }
}

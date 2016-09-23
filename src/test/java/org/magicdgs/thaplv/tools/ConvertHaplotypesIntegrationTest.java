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

import org.magicdgs.thaplv.utils.test.CommandLineProgramTest;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Arrays;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class ConvertHaplotypesIntegrationTest extends CommandLineProgramTest {

    @Test
    public void testConvertHaplotypesHaploid() throws IOException {
        IntegrationTestSpec testSpec = new IntegrationTestSpec(
                " -V " + getTestDataDir() + "example.vcf"
                        + " -H HAPLOID"
                        + " -O %s",
                Arrays.asList(getTestDataDir() + "example_haploid.vcf")
        );
        testSpec.executeTest("testConvertHaplotypesHaploid", this);
    }

    @Test
    public void testConvertHaplotypesOutputPloidy1() throws IOException {
        IntegrationTestSpec testSpec = new IntegrationTestSpec(
                " -V " + getTestDataDir() + "example.vcf"
                        + " -H HAPLOID"
                        + " -ploidy 1"
                        + " -O %s",
                Arrays.asList(getTestDataDir() + "example_haploid_1.vcf")
        );
        testSpec.executeTest("testConvertHaplotypesOutputPloidy1", this);
    }

    @Test
    public void testConvertHaplotypesBackCross() throws IOException {
        IntegrationTestSpec testSpec = new IntegrationTestSpec(
                " -V " + getTestDataDir() + "example.vcf"
                        + " -H BACK_CROSS"
                        + " -O %s",
                Arrays.asList(getTestDataDir()
                        + "example_backcross.vcf")
        );
        testSpec.executeTest("testConvertHaplotypesBackCross", this);
    }

    @Test
    public void testConvertHaplotypesDontCheckError() throws IOException {
        IntegrationTestSpec testSpec = new IntegrationTestSpec(
                " -H DONT_CHECK"
                        + " -V " + getTestDataDir() + "example.vcf"
                        + " -O %s",
                1,
                UserException.class
        );
        testSpec.executeTest("testConvertHaplotypesDontCheckError", this);
    }

    @Test
    public void testConvertHaplotypesWrongPloidy() throws IOException {
        IntegrationTestSpec testSpec = new IntegrationTestSpec(
                " -H HAPLOID"
                        + " -ploidy 10"
                        + " -V " + getTestDataDir() + "example.vcf"
                        + " -O %s",
                1,
                UserException.class
        );
        testSpec.executeTest("testConvertHaplotypesDontCheck", this);
    }

}
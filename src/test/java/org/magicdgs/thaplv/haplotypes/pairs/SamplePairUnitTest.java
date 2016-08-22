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

package org.magicdgs.thaplv.haplotypes.pairs;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class SamplePairUnitTest {

    @DataProvider(name = "badArguments")
    public Object[][] getBadArguments() {
        return new Object[][]{
                {"sample1", null},
                {null, "sample2"},
                {null, null},
                {"sample1", "sample1"}
        };
    }

    @Test(dataProvider = "badArguments", expectedExceptions = IllegalArgumentException.class)
    public void testBadArguments(final String sample1, final String sample2) throws Exception {
        new SamplePair(sample1, sample2);
    }

    @Test
    public void testGetPairNames() throws Exception {
        Assert.assertEquals(new SamplePair("sample1", "sample2").getPairNames(),
                "(sample1,sample2)");
    }

    @DataProvider(name = "samplePairs")
    public Object[][] getSamplesPairs() {
        return new Object[][]{
                {"sample1", "sample2"},
                {"you", "me"},
                {"europe", "africa"}
        };
    }

    @Test(dataProvider = "samplePairs")
    public void testContainSample(final String sample1, final String sample2) throws Exception {
        final SamplePair pair = new SamplePair(sample1, sample2);
        Assert.assertTrue(pair.containSample(sample1));
        Assert.assertTrue(pair.containSample(sample2));
        Assert.assertFalse(pair.containSample("notInSamples"));
    }

    @Test(dataProvider = "samplePairs")
    public void testContainNames(final String sample1, final String sample2) throws Exception {
        final SamplePair pair = new SamplePair(sample1, sample2);
        // contained
        Assert.assertTrue(pair.containNames(sample1, sample2, false));
        Assert.assertTrue(pair.containNames(sample2, sample1, false));
        Assert.assertTrue(pair.containNames(sample1, sample2, true));
        // not contained because of order
        Assert.assertFalse(pair.containNames(sample2, sample1, true));
        // one sample is not in this
        Assert.assertFalse(pair.containNames("notInSamples", sample2, false));
        Assert.assertFalse(pair.containNames(sample1, "notInSamples", false));
        Assert.assertFalse(pair.containNames("notInSamples", sample2, true));
        Assert.assertFalse(pair.containNames(sample1, "notInSamples", true));
    }

}
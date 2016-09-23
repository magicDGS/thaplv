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

package org.magicdgs.thaplv.utils.stats.popgen;

import static org.magicdgs.thaplv.haplotypes.light.LightGenotype.SNP;

import org.magicdgs.thaplv.haplotypes.light.LightGenotype;
import org.magicdgs.thaplv.haplotypes.light.SNPpair;
import org.magicdgs.thaplv.utils.test.BaseTest;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class LDfunctionsUnitTest extends BaseTest {

    @DataProvider(name = "SNPpairs")
    public Object[][] getSNPpairs() throws Exception {
        // SNP_1 A=7; a=4
        final LightGenotype SNP_1 = new LightGenotype("2L", 1,
                new SNP[] {SNP.A, SNP.A, SNP.A, SNP.A, SNP.A, SNP.A, SNP.A, SNP.a, SNP.a, SNP.a,
                        SNP.a});
        // SNP_2 A=7; a=4
        final LightGenotype SNP_2 = new LightGenotype("2L", 2,
                new SNP[] {SNP.A, SNP.a, SNP.a, SNP.A, SNP.A, SNP.A, SNP.A, SNP.A, SNP.a, SNP.a,
                        SNP.A});
        // SNP_3 A=5; a=6
        final LightGenotype SNP_3 = new LightGenotype("2L", 2,
                new SNP[] {SNP.a, SNP.A, SNP.A, SNP.a, SNP.A, SNP.A, SNP.a, SNP.a, SNP.A, SNP.a,
                        SNP.a});
        // singleton SNP A=10; a=1
        final LightGenotype singleton = new LightGenotype("2L", 2,
                new SNP[] {SNP.A, SNP.A, SNP.A, SNP.A, SNP.A, SNP.A, SNP.A, SNP.A, SNP.A, SNP.A,
                        SNP.a});


        return new Object[][] {
                // constructing a fully linked pair of SNPs with snp1
                {new SNPpair(SNP_1, SNP_1), 1d, 1d, 1d},
                // SNP_1 vs. SNP_2: manually computed; AB=5; Ab=2; aB=2; ab=2
                {new SNPpair(SNP_1, SNP_2), 1d, 0.2142857, 0.04591837},
                // SNP_1 vs. SNP_3: manually computed; AB=4; Ab=3; aB=1; ab=3 => AB=3; Ab=4; aB=3; ab=1
                {new SNPpair(SNP_1, SNP_3), 0.6857143, -0.3105295, 0.09642857},
                // SNP_2 vs SNP_3: manually computed; AB=2; Ab=5; aB=3; ab=1 =>  AB=5; Ab=2; aB=1; ab=3
                {new SNPpair(SNP_2, SNP_3), 0.6857143, 0.4485426, 0.2011905},
                // maxR2=same as previous because the frequency is the same
                // SNP_1 vs. singleton: manually computed; AB=7; Ab=0; aB=3; ab=1
                {new SNPpair(SNP_1, singleton), 0.175, 0.41833, 0.175},
                // SNP_2 vs. singleton: manually computed; AB=6; Ab=1; aB=4; ab=0
                {new SNPpair(SNP_2, singleton), 0.175, -0.2390457, 0.05714286},
                // maxR2=same as previous because the frequency is the same
                // SNP_3 vs. singleton: manually computed; AB=5; Ab=0; aB=5; ab=1 => AB=5; aB=5; Ab=1; ab=0
                {new SNPpair(SNP_3, singleton), 0.12, -0.2886751, 0.08333333}

        };

    }

    @Test(dataProvider = "SNPpairs")
    public void testWithSNPPair(final SNPpair pair, final double expectedR2max,
            final double expectedRW, final double expectedR2) {
        Assert.assertEquals(LDfunctions.maxR2(pair), expectedR2max, DEFAULT_TOLERANCE,
                "wrong max. r2");
        Assert.assertEquals(LDfunctions.rw(pair), expectedRW, DEFAULT_TOLERANCE, "wrong rw");
        Assert.assertEquals(LDfunctions.r2(pair), expectedR2, DEFAULT_TOLERANCE, "wrong r2");
    }


    @Test
    public void testIsSignificantR2() throws Exception {
        Assert.assertTrue(LDfunctions.isSignificantR2(0.25, 36, 0.95),
                "r2=0.25 is not significant for 95%-quantile of 36 haplotypes");
        Assert.assertFalse(LDfunctions.isSignificantR2(0.10, 36, 0.95),
                "r2=0.10 is significant for 95%-quantile of 36 haplotypes");
    }

    @DataProvider(name = "significantTest")
    public Object[][] getSignificantThresholdData() {
        // computed with R
        return new Object[][] {
                {28, 0.95, 0.137195},
                {36, 0.95, 0.1067072},
                {205, 0.95, 0.01873882},
                {28, 0.99, 0.2369606},
                {36, 0.99, 0.1843027},
                {205, 0.99, 0.03236535}
        };
    }

    @Test(dataProvider = "significantTest")
    public void testSignificantThreshold(final int nHaplotypes, final double quantile,
            final double expected) throws Exception {
        Assert.assertEquals(LDfunctions.significantThreshold(nHaplotypes, quantile), expected,
                DEFAULT_TOLERANCE);
    }

    @DataProvider(name = "rwData")
    public Object[][] getRwData() {
        return new Object[][] {
                // completely linked
                {7d / 11, 7d / 11, 7d / 11, 1d},
                // positive rw
                {5d / 11, 7d / 11, 7d / 11, 0.2142857},
                {5d / 11, 7d / 11, 6d / 11, 0.4485426},
                // negative rw
                {3d / 11, 7d / 11, 6d / 11, -0.3105295},
                // with singletons
                {7d / 11, 7d / 11, 10d / 11, 0.41833},
                {6d / 11, 7d / 11, 10d / 11, -0.2390457},
                {5d / 11, 6d / 11, 10d / 11, -0.2886751}

        };
    }

    @Test(dataProvider = "rwData")
    public void testRw(final double pAB, final double pA, final double pB, final double expected)
            throws Exception {
        Assert.assertEquals(LDfunctions.rw(pAB, pA, pB), expected, DEFAULT_TOLERANCE);
    }

    @DataProvider(name = "r2Data")
    public Object[][] getR2Data() {
        return new Object[][] {
                // completely linked
                {7d / 11, 7d / 11, 7d / 11, 1d},
                // positive rw
                {5d / 11, 7d / 11, 7d / 11, 0.04591837},
                {5d / 11, 7d / 11, 6d / 11, 0.2011905},
                // negative rw
                {3d / 11, 7d / 11, 6d / 11, 0.09642857},
                // with singletons
                {7d / 11, 7d / 11, 10d / 11, 0.175},
                {6d / 11, 7d / 11, 10d / 11, 0.05714286},
                {5d / 11, 6d / 11, 10d / 11, 0.08333333}

        };
    }


    @Test(dataProvider = "r2Data")
    public void testR2(final double pAB, final double pA, final double pB, final double expected)
            throws Exception {
        Assert.assertEquals(LDfunctions.r2(pAB, pA, pB), expected, DEFAULT_TOLERANCE);
    }

    @DataProvider(name = "maxR2data")
    public Object[][] getMaxR2Data() {
        return new Object[][] {
                // equal frequencies
                {7d / 11, 7d / 11, 1d},
                // sligthly different frequencies
                {7d / 11, 6d / 11, 0.6857143},
                {6d / 11, 7d / 11, 0.6857143},
                // with singletons
                {7d / 11, 10d / 11, 0.175},
                {6d / 11, 10d / 11, 0.12}

        };
    }

    @Test(dataProvider = "maxR2data")
    public void testMaxR2(final double pA, final double pB, final double expected)
            throws Exception {
        Assert.assertEquals(LDfunctions.maxR2(pA, pB), expected, DEFAULT_TOLERANCE);
    }

    @DataProvider(name = "norR2data")
    public Object[][] getNormR2Data() {
        return new Object[][] {
                // completely linked
                {7d / 11, 7d / 11, 7d / 11, 1d},
                // positive rw
                {5d / 11, 7d / 11, 7d / 11, 0.04591837},
                {5d / 11, 7d / 11, 6d / 11, 0.2011905 / 0.6857143},
                // negative rw
                {3d / 11, 7d / 11, 6d / 11, 0.09642857 / 0.6857143},
                // with singletons
                {7d / 11, 7d / 11, 10d / 11, 0.175 / 0.175},
                {6d / 11, 7d / 11, 10d / 11, 0.05714286 / 0.175},
                {5d / 11, 6d / 11, 10d / 11, 0.08333333 / 0.12}

        };
    }

    @Test(dataProvider = "norR2data")
    public void testNormR2(final double pAB, final double pA, final double pB,
            final double expected) throws Exception {
        Assert.assertEquals(LDfunctions.r2norm(pAB, pA, pB), expected, DEFAULT_TOLERANCE);
    }

}
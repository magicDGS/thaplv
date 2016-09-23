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

package org.magicdgs.thaplv.haplotypes.light;

import static org.magicdgs.thaplv.haplotypes.light.LightGenotype.SNP;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 * @deprecated tested class is deprecated
 */
@Deprecated
public class SNPpairTest {

    @DataProvider(name = "pairsToTest")
    public Object[][] pairsToTest() throws Exception {
        final LightGenotype variant1 = new LightGenotype("2L", 1,
                new SNP[] {SNP.A, SNP.A, SNP.A, SNP.a, SNP.a, SNP.a, SNP.N});
        final LightGenotype variant2 = new LightGenotype("2L", 2,
                new SNP[] {SNP.a, SNP.a, SNP.A, SNP.A, SNP.a, SNP.N, SNP.N});
        final LightGenotype variant3 = new LightGenotype("2L", 100,
                new SNP[] {SNP.A, SNP.N, SNP.N, SNP.A, SNP.N, SNP.a, SNP.a});
        final SNPpair pair1_2 = new SNPpair(variant1, variant2);
        final SNPpair pair1_3 = new SNPpair(variant1, variant3);
        final SNPpair pair2_3 = new SNPpair(variant2, variant3);
        final SNPpair pair3_2 = new SNPpair(variant3, variant2);
        return new Object[][] {
                {pair1_2, pair1_3, pair2_3, pair3_2}
        };
    }

    @Test(dataProvider = "pairsToTest")
    public void testGetDistance(final SNPpair pair1_2, final SNPpair pair1_3, final SNPpair pair2_3,
            final SNPpair pair3_2) throws Exception {
        Assert.assertEquals(pair1_2.getDistance(), 1);
        Assert.assertEquals(pair1_3.getDistance(), 99);
        Assert.assertEquals(pair2_3.getDistance(), 98);
        Assert.assertEquals(pair2_3.getDistance(), pair3_2.getDistance());
    }

    @Test(dataProvider = "pairsToTest")
    public void testIsInvariantA(final SNPpair pair1_2, final SNPpair pair1_3,
            final SNPpair pair2_3, final SNPpair pair3_2) throws Exception {
        Assert.assertFalse(pair1_2.isInvariantA());
        Assert.assertFalse(pair1_3.isInvariantA());
        Assert.assertFalse(pair2_3.isInvariantA());
        Assert.assertTrue(pair3_2.isInvariantA());
    }

    @Test(dataProvider = "pairsToTest")
    public void testIsInvariantB(final SNPpair pair1_2, final SNPpair pair1_3,
            final SNPpair pair2_3, final SNPpair pair3_2) throws Exception {
        Assert.assertFalse(pair1_2.isInvariantB());
        Assert.assertFalse(pair1_3.isInvariantB());
        Assert.assertTrue(pair2_3.isInvariantB());
        Assert.assertFalse(pair3_2.isInvariantB());
    }

    @Test(dataProvider = "pairsToTest")
    public void testGetTotalCounts(final SNPpair pair1_2, final SNPpair pair1_3,
            final SNPpair pair2_3, final SNPpair pair3_2) throws Exception {
        Assert.assertEquals(pair1_2.getTotalCounts(), 5);
        Assert.assertEquals(pair1_3.getTotalCounts(), 3);
        Assert.assertEquals(pair2_3.getTotalCounts(), 2);
        Assert.assertEquals(pair2_3.getTotalCounts(), pair2_3.getTotalCounts());
    }

    @Test(dataProvider = "pairsToTest")
    public void testGetMajorAlleleFrequencyForA(final SNPpair pair1_2, final SNPpair pair1_3,
            final SNPpair pair2_3, final SNPpair pair3_2) throws Exception {
        Assert.assertEquals(pair1_2.getMajorAlleleFrequencyForA(), 3d / 5d, Double.MIN_VALUE);
        Assert.assertEquals(pair1_3.getMajorAlleleFrequencyForA(), 2d / 3d, Double.MIN_VALUE);
        Assert.assertEquals(pair2_3.getMajorAlleleFrequencyForA(), 1d / 2d, Double.MIN_VALUE);
        Assert.assertEquals(pair3_2.getMajorAlleleFrequencyForA(), 1d, Double.MIN_VALUE);
    }

    @Test(dataProvider = "pairsToTest")
    public void testGetMajorAlleleFrequencyForB(final SNPpair pair1_2, final SNPpair pair1_3,
            final SNPpair pair2_3, final SNPpair pair3_2) throws Exception {
        Assert.assertEquals(pair1_2.getMajorAlleleFrequencyForB(), 3d / 5d, Double.MIN_VALUE);
        Assert.assertEquals(pair1_3.getMajorAlleleFrequencyForB(), 2d / 3d, Double.MIN_VALUE);
        Assert.assertEquals(pair2_3.getMajorAlleleFrequencyForB(), 1d, Double.MIN_VALUE);
        Assert.assertEquals(pair3_2.getMajorAlleleFrequencyForB(), 1d / 2d, Double.MIN_VALUE);
    }

    @Test(dataProvider = "pairsToTest")
    public void testGetMinorAlleleFrequencyForA(final SNPpair pair1_2, final SNPpair pair1_3,
            final SNPpair pair2_3, final SNPpair pair3_2) throws Exception {
        Assert.assertEquals(pair1_2.getMinorAlleleFrequencyForA(), 2d / 5d, Double.MIN_VALUE);
        Assert.assertEquals(pair1_3.getMinorAlleleFrequencyForA(), 1d / 3d, Double.MIN_VALUE);
        Assert.assertEquals(pair2_3.getMinorAlleleFrequencyForA(), 1d / 2d, Double.MIN_VALUE);
        Assert.assertEquals(pair3_2.getMinorAlleleFrequencyForA(), 0d, Double.MIN_VALUE);
    }

    @Test(dataProvider = "pairsToTest")
    public void testGetMinorAlleleFrequencyForB(final SNPpair pair1_2, final SNPpair pair1_3,
            final SNPpair pair2_3, final SNPpair pair3_2) throws Exception {
        Assert.assertEquals(pair1_2.getMinorAlleleFrequencyForB(), 2d / 5d, Double.MIN_VALUE);
        Assert.assertEquals(pair1_3.getMinorAlleleFrequencyForB(), 1d / 3d, Double.MIN_VALUE);
        Assert.assertEquals(pair2_3.getMinorAlleleFrequencyForB(), 0d, Double.MIN_VALUE);
        Assert.assertEquals(pair3_2.getMinorAlleleFrequencyForB(), 1d / 2d, Double.MIN_VALUE);
    }

    @Test(dataProvider = "pairsToTest")
    public void testOneIsSingleton(final SNPpair pair1_2, final SNPpair pair1_3,
            final SNPpair pair2_3, final SNPpair pair3_2) throws Exception {
        Assert.assertFalse(pair1_2.oneIsSingleton());
        Assert.assertTrue(pair1_3.oneIsSingleton());
        Assert.assertTrue(pair2_3.oneIsSingleton());
        Assert.assertTrue(pair3_2.oneIsSingleton());
    }

}
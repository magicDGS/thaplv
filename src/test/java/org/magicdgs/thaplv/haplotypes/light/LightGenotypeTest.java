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

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.LinkedList;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 * @deprecated tested class is deprecated
 */
@Deprecated
public class LightGenotypeTest {

    private final static int testSize = 100;
    private final static SNP[] result = new SNP[testSize];

    private final static LightGenotype simpleTest;
    private final static LightGenotype fromBiallelic;
    private final static LightGenotype fromMultiAllelic;

    static {
        // create a biallelic and multiallelic alleles
        final GenotypesContext biallelic = GenotypesContext.create();
        final GenotypesContext multiallelic = GenotypesContext.create();
        final Allele ref = Allele.create("A", true);
        final Allele alt1 = Allele.create("T", false);
        final Allele alt2 = Allele.create("C", false);
        for (int i = 0; i < testSize; i++) {
            if (i % 2 == 0) {
                biallelic.add(GenotypeBuilder.create("sample" + i, new LinkedList<Allele>() {{
                    add(ref);
                }}));
                multiallelic.add(GenotypeBuilder.create("sample" + i, new LinkedList<Allele>() {{
                    add(alt1);
                }}));
                result[i] = SNP.A;
            } else if (i % 3 == 0) {
                biallelic.add(
                        GenotypeBuilder.create("sample" + i, new LinkedList<Allele>() {{
                            add(Allele.NO_CALL);
                        }}));
                multiallelic.add(
                        GenotypeBuilder.create("sample" + i, new LinkedList<Allele>() {{
                            add(Allele.NO_CALL);
                        }}));
                result[i] = SNP.N;
            } else {
                biallelic.add(GenotypeBuilder.create("sample" + i, new LinkedList<Allele>() {{
                    add(alt1);
                }}));
                multiallelic.add(GenotypeBuilder.create("sample" + i, new LinkedList<Allele>() {{
                    add(alt2);
                }}));
                result[i] = SNP.a;
            }
        }
        simpleTest = new LightGenotype("2L", 1, new SNP[] {SNP.A, SNP.a, SNP.N});
        fromBiallelic = new LightGenotype("2L", 1, biallelic, ref, new LinkedList<Allele>() {{
            add(alt1);
        }});
        fromMultiAllelic = new LightGenotype("2L", 1, multiallelic, ref, new LinkedList<Allele>() {{
            add(alt1);
            add(alt2);
        }});
    }


    @Test
    public void testGetNumberOfA() throws Exception {
        Assert.assertEquals(simpleTest.getNumberOfA(), 1);
    }

    @Test
    public void testGetNumberOfa() throws Exception {
        Assert.assertEquals(simpleTest.getNumberOfa(), 1);
    }

    @Test
    public void testGetNumberOfMissing() throws Exception {
        Assert.assertEquals(simpleTest.getNumberOfMissing(), 1);
    }

    @Test
    public void testSize() throws Exception {
        Assert.assertEquals(simpleTest.size(), 3);
        Assert.assertEquals(fromBiallelic.size(), testSize);
        Assert.assertEquals(fromMultiAllelic.size(), testSize);
    }

    @Test
    public void testGetGenotypeAt() throws Exception {
        Assert.assertEquals(simpleTest.getGenotypeAt(0), SNP.A);
        Assert.assertEquals(simpleTest.getGenotypeAt(1), SNP.a);
        Assert.assertEquals(simpleTest.getGenotypeAt(2), SNP.N);
    }

    @Test
    public void testConstructors() throws Exception {
        int countA = 0;
        int countN = 0;
        // testing if the constructor is making the correct result
        for (int i = 0; i < testSize; i++) {
            Assert.assertEquals(fromBiallelic.getGenotypeAt(i), result[i]);
            Assert.assertEquals(fromMultiAllelic.getGenotypeAt(i), result[i]);
            if (i % 2 == 0) {
                countA++;
            } else if (i % 3 == 0) {
                countN++;
            }
        }
        // checking the counts
        Assert.assertEquals(fromBiallelic.getNumberOfA(), countA);
        Assert.assertEquals(fromBiallelic.getNumberOfa(), testSize - countA - countN);
        Assert.assertEquals(fromBiallelic.getNumberOfMissing(), countN);
        Assert.assertEquals(fromMultiAllelic.getNumberOfA(), countA);
        Assert.assertEquals(fromMultiAllelic.getNumberOfa(), testSize - countA - countN);
        Assert.assertEquals(fromMultiAllelic.getNumberOfMissing(), countN);
    }

    @Test(expectedExceptions = ArrayIndexOutOfBoundsException.class)
    public void testGetGenotypeAtIndexException() throws Exception {
        simpleTest.getGenotypeAt(3);
    }
}
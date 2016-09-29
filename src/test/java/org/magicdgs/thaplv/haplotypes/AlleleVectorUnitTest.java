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

package org.magicdgs.thaplv.haplotypes;

import org.magicdgs.thaplv.utils.test.BaseTest;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class AlleleVectorUnitTest extends BaseTest {

    private final Allele refAllele = Allele.create((byte) 'A', true);
    private final Allele altT = Allele.create((byte) 'T');
    private final Allele altG = Allele.create((byte) 'G');

    @DataProvider
    public Object[][] badAlleles() {
        return new Object[][] {
                {null},
                {new Allele[0]}
        };
    }

    @Test(dataProvider = "badAlleles", expectedExceptions = IllegalArgumentException.class)
    public void testIllegalAlleles(final Allele[] alleles) throws Exception {
        new AlleleVector(alleles);
    }

    @DataProvider
    public Object[][] badGenotypes() {
        return new Object[][] {
                {GenotypesContext.create()},
                {null}
        };
    }

    @Test(dataProvider = "badGenotypes", expectedExceptions = IllegalArgumentException.class)
    public void testIllegalGenotypes(final List<Genotype> genotypes) throws Exception {
        new AlleleVector(genotypes);
    }

    @DataProvider
    public Object[][] badVariants() {
        return new Object[][] {
                {new VariantContextBuilder("source", "1", 1, 1, Collections.singleton(refAllele))
                        .make()},
                {null}
        };
    }

    @Test(dataProvider = "badVariants", expectedExceptions = IllegalArgumentException.class)
    public void testIllegalVariants(final VariantContext variant) throws Exception {
        new AlleleVector(variant);
    }

    @Test
    public void testSimpleAlleleVector() throws Exception {
        final Allele[] alleleArray =
                new Allele[] {altT, altT, refAllele, refAllele, altG, Allele.NO_CALL};
        final AlleleVector vector = new AlleleVector(alleleArray);
        // check if retrieving the alleles are the same
        Assert.assertEquals(vector.getAlleles(), alleleArray);
        Assert.assertNotSame(vector.getAlleles(), alleleArray);
        // test simple methods
        Assert.assertEquals(vector.numberOfAlleles(), 3);
        Assert.assertEquals(vector.size(), alleleArray.length);
        for (int i = 0; i < alleleArray.length; i++) {
            Assert.assertEquals(vector.getAlleleAt(i), alleleArray[i]);
        }
        // create an expected allele and check that everything is equals
        final AlleleVector expected =
                new AlleleVector(new int[] {1, 1, 0, 0, 2, -1},
                        Arrays.asList(refAllele, altT, altG));
        Assert.assertEquals(vector, expected);
        Assert.assertTrue(vector.equalsState(expected));
    }

    @DataProvider(name = "biallelic")
    public Object[][] getBiallelicData() {
        return new Object[][] {
                {new AlleleVector(new Allele[] {refAllele, refAllele}),
                        new AlleleVector(new int[] {0, 0}, Collections.singletonList(refAllele))},
                {new AlleleVector(new Allele[] {altG, altG}),
                        new AlleleVector(new int[] {0, 00}, Collections.singletonList(altG))},
                {new AlleleVector(new Allele[] {Allele.NO_CALL, Allele.NO_CALL}),
                        new AlleleVector(new int[] {-1, -1}, Collections.emptyList())},
                {new AlleleVector(new Allele[] {refAllele, altG}),
                        new AlleleVector(new int[] {0, 1}, Arrays.asList(refAllele, altG))},
                {new AlleleVector(new Allele[] {refAllele, altG, altG}),
                        new AlleleVector(new int[] {1, 0, 0}, Arrays.asList(altG, refAllele))},
        };
    }

    @Test(dataProvider = "biallelic")
    public void testToBiallelicNoTransform(final AlleleVector vector, final AlleleVector expected)
            throws Exception {
        Assert.assertEquals(vector, expected);
        Assert.assertEquals(vector.asBiallelic(), expected);
        Assert.assertSame(vector.asBiallelic(), vector);
    }

    @DataProvider(name = "nonBiallelic")
    public Object[][] getNoBiallelicData() {
        return new Object[][] {
                {new AlleleVector(new Allele[] {refAllele, altG, altT}),
                        new AlleleVector(new int[] {0, 1, -1},
                                Arrays.asList(refAllele, altG))},
                {new AlleleVector(new Allele[] {refAllele, altG, altT, altG}),
                        new AlleleVector(new int[] {1, 0, -1, 0},
                                Arrays.asList(refAllele, altG))},
                {new AlleleVector(new Allele[] {refAllele, altG, altT, altG, Allele.NO_CALL}),
                        new AlleleVector(new int[] {1, 0, -1, 0, -1},
                                Arrays.asList(refAllele, altG))}
        };
    }

    @Test(dataProvider = "nonBiallelic")
    public void testToBiallelicTransform(final AlleleVector vector, final AlleleVector expected)
            throws Exception {
        Assert.assertNotEquals(vector, expected);
        Assert.assertEquals(vector.asBiallelic(), expected);
        Assert.assertNotSame(vector.asBiallelic(), vector);
        Assert.assertNotSame(vector.asBiallelic(), expected);
    }

    @Test
    public void testEqualsButNotInState() throws Exception {
        final AlleleVector vector1 =
                new AlleleVector(new Allele[] {refAllele, Allele.NO_CALL, refAllele, altG, altT});
        final AlleleVector vector2 =
                new AlleleVector(new Allele[] {altG, Allele.NO_CALL, altG, refAllele, altT});
        Assert.assertEquals(vector1, vector2);
        Assert.assertFalse(vector1.equalsState(vector2));
        Assert.assertFalse(vector2.equalsState(vector1));
    }

}
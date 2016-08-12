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

package org.magicdgs.thaplv.haplotypes.model;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class InbredLineUnitTest extends HaplotypeConverterUnitTest {

    private static AbstractHaplotypeConverter haploid = new InbredLine(1);
    private static AbstractHaplotypeConverter diploid = new InbredLine(2);

    @DataProvider(name = "genotypesData")
    public Object[][] getTestGenotypes() {
        return new Object[][] {
                // missing
                {GenotypeBuilder.createMissing("test", 2), true, Allele.NO_CALL},
                // homozygous ref
                {GenotypeBuilder.create("test", Arrays.asList(refA, refA)), false, refA},
                // homozygous alt
                {GenotypeBuilder.create("test", Arrays.asList(altT, altT)), false, altT},
                {GenotypeBuilder.create("test", Arrays.asList(altC, altC)), false, altC},
                // heterozygous ref
                {GenotypeBuilder.create("test", Arrays.asList(refA, altT)), true, Allele.NO_CALL},
                {GenotypeBuilder.create("test", Arrays.asList(refA, altC)), true, Allele.NO_CALL},
                // heterozygous alt
                {GenotypeBuilder.create("test", Arrays.asList(altC, altT)), true, Allele.NO_CALL},
                // heterozygous no-call
                {GenotypeBuilder.create("test", Arrays.asList(refA, Allele.NO_CALL)), false, refA},
                {GenotypeBuilder.create("test", Arrays.asList(altT, Allele.NO_CALL)), false, altT},
                {GenotypeBuilder.create("test", Arrays.asList(altC, Allele.NO_CALL)), false, altC},
                {GenotypeBuilder.create("test", Arrays.asList(Allele.NO_CALL, refA)), false, refA},
                {GenotypeBuilder.create("test", Arrays.asList(Allele.NO_CALL, altT)), false, altT},
                {GenotypeBuilder.create("test", Arrays.asList(Allele.NO_CALL, altC)), false, altC}
        };
    }


    @Test(dataProvider = "genotypesData")
    public void testModelMethods(final Genotype genotype, final boolean isMissing,
            final Allele expectedAllele) throws Exception {
        testIsConsideredMissing(haploid, genotype, isMissing, "wrong missing");
        if (isMissing) {
            Assert.assertEquals(expectedAllele, Allele.NO_CALL, "wrong data provider");
        } else {
            testGetHaplotypeAllele(haploid, genotype.getAllele(0), genotype.getAllele(1),
                    expectedAllele, "wrong allele");
        }
    }

    @Test(dataProvider = "genotypesData")
    public void testApply(final Genotype genotype, final boolean isMissing,
            final Allele expectedAllele) throws Exception {
        // haploid expected
        Genotype expected = new GenotypeBuilder(genotype).phased(true)
                .alleles(Arrays.asList(expectedAllele)).make();
        assertGenotypesEquals(haploid.apply(genotype), expected, false);
        // expected diploid
        expected = new GenotypeBuilder(genotype).phased(true)
                .alleles(Arrays.asList(expectedAllele, expectedAllele)).make();
        assertGenotypesEquals(diploid.apply(genotype), expected, true);
    }

    @Test(dataProvider = "haploidData", dataProviderClass = HaplotypeConverterUnitTest.class)
    public void testApplyHaploids(final Genotype genotype, final Allele allele) throws Exception {
        Genotype expected = new GenotypeBuilder(genotype).phased(true)
                .alleles(Arrays.asList(allele)).make();
        assertGenotypesEquals(haploid.apply(genotype), expected, false);
        // expected diploid
        expected = new GenotypeBuilder(genotype).phased(true)
                .alleles(Arrays.asList(allele, allele)).make();
        assertGenotypesEquals(diploid.apply(genotype), expected, true);
    }
}
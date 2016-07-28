/*
 * Copyright (c) 2016, Daniel Gomez-Sanchez <daniel.gomez.sanchez@hotmail All rights reserved.
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
import java.util.Collections;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class BackCrossUnitTest extends HaplotypeConverterUnitTest {

    private static AbstractHaplotypeConverter conservative = new BackCross(1);
    private static AbstractHaplotypeConverter conservativeDiploid = new BackCross(1);
    private static AbstractHaplotypeConverter relaxed = new BackCross(1, false);
    private static AbstractHaplotypeConverter relaxedDiploid = new BackCross(1, false);

    @DataProvider(name = "genotypesData")
    public Object[][] getTestGenotypes() {
        return new Object[][] {
                // missing
                {GenotypeBuilder.createMissing("test", 2), true, true,
                        Allele.NO_CALL, Allele.NO_CALL},
                // homozygous ref
                {GenotypeBuilder.create("test", Arrays.asList(refA, refA)), false, false,
                        refA, refA},
                // homozygous alt
                {GenotypeBuilder.create("test", Arrays.asList(altT, altT)), true, false,
                        Allele.NO_CALL, altT},
                {GenotypeBuilder.create("test", Arrays.asList(altC, altC)), true, false,
                        Allele.NO_CALL, altC},
                // heterozygous ref
                {GenotypeBuilder.create("test", Arrays.asList(refA, altT)), false, false,
                        altT, altT},
                {GenotypeBuilder.create("test", Arrays.asList(refA, altC)), false, false,
                        altC, altC},
                // heterozygous alt
                {GenotypeBuilder.create("test", Arrays.asList(altC, altT)), true, true,
                        Allele.NO_CALL, Allele.NO_CALL},
                // heterozygous no-call are always no call
                {GenotypeBuilder.create("test", Arrays.asList(refA, Allele.NO_CALL)), false, false,
                        Allele.NO_CALL, Allele.NO_CALL},
                {GenotypeBuilder.create("test", Arrays.asList(altT, Allele.NO_CALL)), false, false,
                        Allele.NO_CALL, Allele.NO_CALL},
                {GenotypeBuilder.create("test", Arrays.asList(altC, Allele.NO_CALL)), false, false,
                        Allele.NO_CALL, Allele.NO_CALL},
                {GenotypeBuilder.create("test", Arrays.asList(Allele.NO_CALL, refA)), false, false,
                        Allele.NO_CALL, Allele.NO_CALL},
                {GenotypeBuilder.create("test", Arrays.asList(Allele.NO_CALL, altT)), false, false,
                        Allele.NO_CALL, Allele.NO_CALL},
                {GenotypeBuilder.create("test", Arrays.asList(Allele.NO_CALL, altC)), false, false,
                        Allele.NO_CALL, Allele.NO_CALL}
        };
    }

    @Test(dataProvider = "genotypesData")
    public void testConservative(final Genotype genotype, final boolean isMissingConservative,
            final boolean isMissingRelaxed, final Allele expectedAlleleConservative,
            final Allele expectedAlleleRelaxed) throws Exception {
        testModelMethods(conservative, genotype, isMissingConservative, expectedAlleleConservative);
    }

    @Test(dataProvider = "genotypesData")
    public void testRelaxed(final Genotype genotype, final boolean isMissingConservative,
            final boolean isMissingRelaxed, final Allele expectedAlleleConservative,
            final Allele expectedAlleleRelaxed) throws Exception {
        testModelMethods(relaxed, genotype, isMissingRelaxed, expectedAlleleRelaxed);
    }

    private void testModelMethods(final AbstractHaplotypeConverter converter,
            final Genotype genotype, final boolean isMissing, final Allele expectedAllele)
            throws Exception {
        testIsConsideredMissing(converter, genotype, isMissing, "wrong missing");
        if (isMissing) {
            Assert.assertEquals(expectedAllele, Allele.NO_CALL, "wrong data provider");
        } else {
            testGetHaplotypeAllele(converter, genotype.getAllele(0), genotype.getAllele(1),
                    expectedAllele, "wrong allele");
        }
    }

    @Test(dataProvider = "genotypesData")
    public void testConservativeApply(final Genotype genotype, final boolean isMissingConservative,
            final boolean isMissingRelaxed, final Allele expectedAlleleConservative,
            final Allele expectedAlleleRelaxed) throws Exception {
        // haploid
        testApply(conservative, conservative.getOutputPlotidy(), genotype,
                expectedAlleleConservative);
        // diploid
        testApply(conservativeDiploid, conservativeDiploid.getOutputPlotidy(), genotype,
                expectedAlleleConservative);
    }

    @Test(dataProvider = "genotypesData")
    public void testRelaxedApply(final Genotype genotype, final boolean isMissingConservative,
            final boolean isMissingRelaxed, final Allele expectedAlleleConservative,
            final Allele expectedAlleleRelaxed) throws Exception {
        // haploid
        testApply(relaxed, relaxed.getOutputPlotidy(), genotype, expectedAlleleRelaxed);
        // diploid
        testApply(relaxedDiploid, relaxedDiploid.getOutputPlotidy(), genotype,
                expectedAlleleRelaxed);
    }

    private void testApply(final AbstractHaplotypeConverter converter, final int ploidy,
            final Genotype genotype, final Allele expectedAllele) throws Exception {
        final Genotype expected = new GenotypeBuilder(genotype).phased(true)
                .alleles(Collections.nCopies(ploidy, expectedAllele)).make();
        assertGenotypesEquals(converter.apply(genotype), expected, true);
    }

    @Test(dataProvider = "haploidData", dataProviderClass = HaplotypeConverterUnitTest.class)
    public void testApplyHaploids(final Genotype genotype, final Allele allele) throws Exception {
        // conservative
        testApply(conservative, conservative.getOutputPlotidy(), genotype, allele);
        testApply(conservativeDiploid, conservativeDiploid.getOutputPlotidy(), genotype, allele);
        testApply(relaxed, relaxed.getOutputPlotidy(), genotype, allele);
        testApply(relaxed, relaxedDiploid.getOutputPlotidy(), genotype, allele);
    }

}
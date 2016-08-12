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

import org.magicdgs.thaplv.utils.test.BaseTest;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;

import java.util.Collections;

/**
 * Base test for all haplotype converters, with helper methods and data.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class HaplotypeConverterUnitTest extends BaseTest {

    public static final Allele refA = Allele.create("A", true);
    public static final Allele altT = Allele.create("T", false);
    public static final Allele altC = Allele.create("C", false);

    /**
     * Already haploids should return the same allele
     */
    @DataProvider(name = "haploidData")
    public static Object[][] getTestHaploidGenotypes() {
        return new Object[][] {
                // test haploids
                {GenotypeBuilder.create("test", Collections.singletonList(refA)), refA},
                {GenotypeBuilder.create("test", Collections.singletonList(altC)), altC},
                {GenotypeBuilder.create("test", Collections.singletonList(altT)), altT},
                {GenotypeBuilder.create("test", Collections.singletonList(Allele.NO_CALL)),
                        Allele.NO_CALL}
        };
    }

    public static void testGetHaplotypeAllele(final AbstractHaplotypeConverter converter,
            final Allele allele1, final Allele allele2, final Allele expected, final String msg)
            throws Exception {
        Assert.assertEquals(converter.getHaplotypeAllele(allele1, allele2), expected, msg);
    }

    public static void testIsConsideredMissing(final AbstractHaplotypeConverter converter,
            final Genotype genotype, final boolean isMissing, final String msg) throws Exception {
        Assert.assertEquals(converter.isConsideredMissing(genotype), isMissing, msg);
    }
}
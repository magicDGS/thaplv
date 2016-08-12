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

package org.magicdgs.thaplv.haplotypes.filters;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.engine.filters.VariantFilter;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class NumberOfMissingFilterUnitTest extends VariantFilterBaseTest {

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testBadArgument() {
        new NumberOfMissingFilter(-1);
    }

    @DataProvider(name = "missingData")
    public Object[][] getMissingData() {
        final VariantContextBuilder builder = new VariantContextBuilder().chr("1").start(1).stop(1)
                .alleles(Arrays.asList(refA, altT, altC));

        final VariantFilter noMissing = new NumberOfMissingFilter(0);
        final VariantFilter oneMissing = new NumberOfMissingFilter(2);
        final VariantFilter neverFilter = new NumberOfMissingFilter(Integer.MAX_VALUE);

        return new Object[][] {
                // passing if there is no missing
                {noMissing, builder.genotypes(getHmGenotypes(refA, 5)).make(), true},
                {oneMissing, builder.genotypes(getHmGenotypes(refA, 5)).make(), true},
                {neverFilter, builder.genotypes(getHmGenotypes(refA, 5)).make(), true},
                // if there is missing
                {noMissing, builder.genotypes(Stream.concat(
                        getHmGenotypes(refA, 10).stream(),
                        getHmGenotypes(Allele.NO_CALL, 2).stream())
                        .collect(Collectors.toList())).make(), false},
                {oneMissing, builder.genotypes(Stream.concat(
                        getHmGenotypes(refA, 10).stream(),
                        getHmGenotypes(Allele.NO_CALL, 2).stream())
                        .collect(Collectors.toList())).make(), true},
                {neverFilter, builder.genotypes(Stream.concat(
                        getHmGenotypes(refA, 10).stream(),
                        getHmGenotypes(Allele.NO_CALL, 2).stream())
                        .collect(Collectors.toList())).make(), true},
                {noMissing, builder.genotypes(Stream.concat(
                        getHmGenotypes(refA, 10).stream(),
                        getHmGenotypes(Allele.NO_CALL, 3).stream())
                        .collect(Collectors.toList())).make(), false},
                {oneMissing, builder.genotypes(Stream.concat(
                        getHmGenotypes(refA, 10).stream(),
                        getHmGenotypes(Allele.NO_CALL, 3).stream())
                        .collect(Collectors.toList())).make(), false},
                {neverFilter, builder.genotypes(Stream.concat(
                        getHmGenotypes(refA, 10).stream(),
                        getHmGenotypes(Allele.NO_CALL, 3).stream())
                        .collect(Collectors.toList())).make(), true}
        };
    }

    @Test(dataProvider = "missingData")
    public void testTest(final NumberOfMissingFilter filter, final VariantContext context,
            final boolean pass) throws Exception {
        Assert.assertEquals(filter.test(context), pass);
    }

    @DataProvider(name = "badArguments")
    public Object[][] getBadArgumentsFromMaxPercentage() {
        return new Object[][] {{-1}, {102}};
    }

    @Test(dataProvider = "badArguments", expectedExceptions = IllegalArgumentException.class)
    public void testFromMaxPercentageBadArguments(final int maxPercentage) throws Exception {
        NumberOfMissingFilter.fromMaxPercentage(maxPercentage, 10);
    }

    @DataProvider(name = "fromMaxPercentage")
    public Object[][] getFromMaxPercentageData() {
        return new Object[][] {
                {0, 100, new NumberOfMissingFilter(0)},
                {10, 100, new NumberOfMissingFilter(10)},
                {1, 100, new NumberOfMissingFilter(1)}
        };
    }

    @Test(dataProvider = "fromMaxPercentage")
    public void testFromMaxPercentage(final int maxPercentage, final int nSamples,
            final NumberOfMissingFilter expected) throws Exception {
        Assert.assertEquals(NumberOfMissingFilter.fromMaxPercentage(maxPercentage, nSamples),
                expected);
    }

}
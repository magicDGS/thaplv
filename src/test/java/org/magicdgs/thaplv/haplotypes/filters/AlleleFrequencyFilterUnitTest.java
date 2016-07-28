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

package org.magicdgs.thaplv.haplotypes.filters;

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
public class AlleleFrequencyFilterUnitTest extends VariantFilterBaseTest {

    @DataProvider(name = "badArguments")
    public Object[][] getBadArguments() {
        return new Object[][] {
                {-1, 1},
                {0, 10},
                {10, 0.9},
                {0.5, -1},
                {0.5, 0.3},
                {0.9, 0.9}
        };
    }

    @Test(dataProvider = "badArguments", expectedExceptions = IllegalArgumentException.class)
    public void testBadConstruction(final double minFreq, final double maxFreq) throws Exception {
        new AlleleFrequencyFilter(minFreq, maxFreq);
    }

    @DataProvider(name = "frequencyFilter")
    public Object[][] getFrequencyFilterData() throws Exception {
        final VariantContextBuilder builder = new VariantContextBuilder().chr("1").start(1).stop(1)
                .alleles(Arrays.asList(refA, altT, altC));

        final VariantFilter lowFrequencyFilter = new AlleleFrequencyFilter(0.4, 1);
        final VariantFilter highFrequencyFilter = new AlleleFrequencyFilter(0, 0.7);
        final VariantFilter bothFilters = new AlleleFrequencyFilter(0.4, 0.7);

        return new Object[][] {
                // fixed alleles
                {lowFrequencyFilter, builder.genotypes(getHmGenotypes(refA, 5)).make(), true},
                {highFrequencyFilter, builder.genotypes(getHmGenotypes(refA, 5)).make(), false},
                {bothFilters, builder.genotypes(getHmGenotypes(refA, 5)).make(), false},
                {lowFrequencyFilter, builder.genotypes(getHmGenotypes(altC, 5)).make(), true},
                {highFrequencyFilter, builder.genotypes(getHmGenotypes(altC, 5)).make(), false},
                {bothFilters, builder.genotypes(getHmGenotypes(altC, 5)).make(), false},
                // minor frequency
                {lowFrequencyFilter, builder.genotypes(Stream.concat(
                        getHmGenotypes(refA, 2).stream(),
                        getHmGenotypes(altC, 1).stream())
                        .collect(Collectors.toList())).make(), false},
                {highFrequencyFilter, builder.genotypes(Stream.concat(
                        getHmGenotypes(refA, 2).stream(),
                        getHmGenotypes(altC, 1).stream())
                        .collect(Collectors.toList())).make(), true},
                {lowFrequencyFilter, builder.genotypes(Stream.concat(
                        getHmGenotypes(refA, 2).stream(),
                        getHmGenotypes(altC, 1).stream())
                        .collect(Collectors.toList())).make(), false}
        };

    }

    @Test(dataProvider = "frequencyFilter")
    public void testTest(final VariantFilter filter, final VariantContext context,
            final boolean pass) throws Exception {
        Assert.assertEquals(filter.test(context), pass);
    }

}
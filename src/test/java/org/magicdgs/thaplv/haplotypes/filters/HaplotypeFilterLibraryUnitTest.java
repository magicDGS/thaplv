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
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class HaplotypeFilterLibraryUnitTest extends VariantFilterBaseTest {

    @DataProvider(name = "biallelicFilter")
    public Object[][] getBiallelicFilterData() throws Exception {
        final VariantContextBuilder builder = new VariantContextBuilder().chr("1").start(1).stop(1)
                .alleles(Arrays.asList(refA, altT, altC));
        return new Object[][] {
                // homozygous don't pass
                {builder.genotypes(getHmGenotypes(refA, 10)).make(), false},
                {builder.genotypes(getHmGenotypes(altT, 10)).make(), false},
                {builder.genotypes(getHmGenotypes(altC, 10)).make(), false},
                // heterozygous don't not pass (because only the first allele)
                {builder.genotypes(getHzGenotypes(refA, altT, 10)).make(), false},
                {builder.genotypes(getHzGenotypes(altC, altT, 10)).make(), false},
                // homozygous of two classes pass
                {builder.genotypes(Stream.concat(getHmGenotypes(refA, 10).stream(),
                        getHmGenotypes(altC, 10).stream())
                        .collect(Collectors.toList())).make(), true},
                {builder.genotypes(Stream.concat(getHmGenotypes(altT, 10).stream(),
                        getHmGenotypes(altC, 10).stream())
                        .collect(Collectors.toList())).make(), true},
                // including some missing samples do not pass
                {builder.genotypes(Stream.concat(getHmGenotypes(refA, 5).stream(),
                        getHmGenotypes(Allele.NO_CALL, 5).stream())
                        .collect(Collectors.toList())).make(), false},
                {builder.genotypes(Stream.concat(getHmGenotypes(altT, 5).stream(),
                        getHmGenotypes(Allele.NO_CALL, 5).stream())
                        .collect(Collectors.toList())).make(), false},
                // triallelic don't pass
                {builder.genotypes(Stream.concat(
                        Stream.concat(getHmGenotypes(refA, 3).stream(),
                                getHmGenotypes(altT, 3).stream()),
                        getHmGenotypes(altC, 3).stream())
                        .collect(Collectors.toList())).make(), false}

        };
    }

    @Test(dataProvider = "biallelicFilter")
    public void testBiallelicFilter(final VariantContext context, final boolean pass)
            throws Exception {
        Assert.assertEquals(HaplotypeFilterLibrary.BIALLELIC_FILTER.test(context), pass);
    }

    @DataProvider(name = "polymorphicFilter")
    public Object[][] getPolymorphicFilterData() throws Exception {
        final VariantContextBuilder builder = new VariantContextBuilder().chr("1").start(1).stop(1)
                .alleles(Arrays.asList(refA, altT, altC));
        return new Object[][] {
                // homozygous don't pass
                {builder.genotypes(getHmGenotypes(refA, 10)).make(), false},
                {builder.genotypes(getHmGenotypes(altT, 10)).make(), false},
                {builder.genotypes(getHmGenotypes(altC, 10)).make(), false},
                // heterozygous don't not pass (because only the first allele)
                {builder.genotypes(getHzGenotypes(refA, altT, 10)).make(), false},
                {builder.genotypes(getHzGenotypes(altC, altT, 10)).make(), false},
                // homozygous of two classes pass
                {builder.genotypes(Stream.concat(getHmGenotypes(refA, 10).stream(),
                        getHmGenotypes(altC, 10).stream())
                        .collect(Collectors.toList())).make(), true},
                {builder.genotypes(Stream.concat(getHmGenotypes(altT, 10).stream(),
                        getHmGenotypes(altC, 10).stream())
                        .collect(Collectors.toList())).make(), true},
                // including some missing samples do not pass
                {builder.genotypes(Stream.concat(getHmGenotypes(refA, 10).stream(),
                        getHmGenotypes(Allele.NO_CALL, 10).stream())
                        .collect(Collectors.toList())).make(), false},
                {builder.genotypes(Stream.concat(getHmGenotypes(altT, 10).stream(),
                        getHmGenotypes(Allele.NO_CALL, 10).stream())
                        .collect(Collectors.toList())).make(), false},
                // triallelic pass
                {builder.genotypes(Stream.concat(
                        Stream.concat(getHmGenotypes(refA, 3).stream(),
                                getHmGenotypes(altT, 3).stream()),
                        getHmGenotypes(altC, 3).stream())
                        .collect(Collectors.toList())).make(), true}

        };
    }

    @Test(dataProvider = "polymorphicFilter")
    public void testPolymorphicFilter(final VariantContext context, final boolean pass)
            throws Exception {
        Assert.assertEquals(HaplotypeFilterLibrary.POLYMORPHIC_FILTER.test(context), pass);
    }

    @DataProvider(name = "noSingletonFilter")
    public Object[][] getNoSingletonFilterData() throws Exception {
        final VariantContextBuilder builder = new VariantContextBuilder().chr("1").start(1).stop(1)
                .alleles(Arrays.asList(refA, altT, altC));
        return new Object[][] {
                // fixed or not, they pass if there is more than 1
                // two homozygous (or heterozygous because of the first allele)
                {builder.genotypes(getHmGenotypes(refA, 2)).make(), true},
                {builder.genotypes(getHmGenotypes(altT, 2)).make(), true},
                {builder.genotypes(getHmGenotypes(altC, 2)).make(), true},
                {builder.genotypes(getHzGenotypes(refA, altT, 2)).make(), true},
                {builder.genotypes(getHzGenotypes(altC, altT, 2)).make(), true},
                // two individuals of each
                {builder.genotypes(Stream.concat(getHmGenotypes(refA, 2).stream(),
                        getHmGenotypes(altC, 2).stream())
                        .collect(Collectors.toList())).make(), true},
                {builder.genotypes(Stream.concat(getHmGenotypes(altT, 2).stream(),
                        getHmGenotypes(altC, 2).stream())
                        .collect(Collectors.toList())).make(), true},
                // including one genotype with no call
                {builder.genotypes(Stream.concat(getHmGenotypes(refA, 2).stream(),
                        getHmGenotypes(Allele.NO_CALL, 2).stream())
                        .collect(Collectors.toList())).make(), true},
                {builder.genotypes(Stream.concat(getHmGenotypes(altT, 2).stream(),
                        getHmGenotypes(Allele.NO_CALL, 2).stream())
                        .collect(Collectors.toList())).make(), true},
                // even if the NO_CALL only appears one
                {builder.genotypes(Stream.concat(getHmGenotypes(refA, 2).stream(),
                        getHmGenotypes(Allele.NO_CALL, 1).stream())
                        .collect(Collectors.toList())).make(), true},
                {builder.genotypes(Stream.concat(getHmGenotypes(altT, 2).stream(),
                        getHmGenotypes(Allele.NO_CALL, 1).stream())
                        .collect(Collectors.toList())).make(), true},
                // three alleles in haplotypes
                {builder.genotypes(Stream.concat(
                        Stream.concat(getHmGenotypes(refA, 2).stream(),
                                getHmGenotypes(altT, 2).stream()),
                        getHmGenotypes(altC, 2).stream())
                        .collect(Collectors.toList())).make(), true},
                // singletons, they do not pass
                // only one haplotype is a singleton by definition
                {builder.genotypes(getHmGenotypes(refA, 1)).make(), false},
                {builder.genotypes(getHmGenotypes(altT, 1)).make(), false},
                {builder.genotypes(getHmGenotypes(altC, 1)).make(), false},
                {builder.genotypes(getHzGenotypes(refA, altT, 1)).make(), false},
                {builder.genotypes(getHzGenotypes(altC, altT, 1)).make(), false},
                // two haplotypes with different alleles, are two singletons
                {builder.genotypes(Stream.concat(getHmGenotypes(refA, 1).stream(),
                        getHmGenotypes(altC, 1).stream())
                        .collect(Collectors.toList())).make(), false},
                {builder.genotypes(Stream.concat(getHmGenotypes(altT, 1).stream(),
                        getHmGenotypes(altC, 1).stream())
                        .collect(Collectors.toList())).make(), false},
                {builder.genotypes(Stream.concat(getHmGenotypes(refA, 1).stream(),
                        getHmGenotypes(Allele.NO_CALL, 1).stream())
                        .collect(Collectors.toList())).make(), false},
                {builder.genotypes(Stream.concat(getHmGenotypes(altT, 1).stream(),
                        getHmGenotypes(Allele.NO_CALL, 1).stream())
                        .collect(Collectors.toList())).make(), false},
                // and also three individuals, one of them with a different allele
                {builder.genotypes(Stream.concat(getHmGenotypes(refA, 2).stream(),
                        getHmGenotypes(altC, 1).stream())
                        .collect(Collectors.toList())).make(), false},
                {builder.genotypes(Stream.concat(getHmGenotypes(altT, 2).stream(),
                        getHmGenotypes(altC, 1).stream())
                        .collect(Collectors.toList())).make(), false},
                {builder.genotypes(Stream.concat(getHmGenotypes(refA, 1).stream(),
                        getHmGenotypes(Allele.NO_CALL, 2).stream())
                        .collect(Collectors.toList())).make(), false},
                {builder.genotypes(Stream.concat(getHmGenotypes(altT, 1).stream(),
                        getHmGenotypes(Allele.NO_CALL, 2).stream())
                        .collect(Collectors.toList())).make(), false},
                // or in the case of tri-allelic with one of them with only one haplotype
                {builder.genotypes(Stream.concat(
                        Stream.concat(getHmGenotypes(refA, 1).stream(),
                                getHmGenotypes(altT, 2).stream()),
                        getHmGenotypes(altC, 2).stream())
                        .collect(Collectors.toList())).make(), false},
                {builder.genotypes(Stream.concat(
                        Stream.concat(getHmGenotypes(refA, 2).stream(),
                                getHmGenotypes(altT, 1).stream()),
                        getHmGenotypes(altC, 2).stream())
                        .collect(Collectors.toList())).make(), false},
                {builder.genotypes(Stream.concat(
                        Stream.concat(getHmGenotypes(refA, 2).stream(),
                                getHmGenotypes(altT, 2).stream()),
                        getHmGenotypes(altC, 1).stream())
                        .collect(Collectors.toList())).make(), false},
                // and only mising does not pass
                {builder.genotypes(getHmGenotypes(Allele.NO_CALL, 5)).make(), false}
        };
    }

    @Test(dataProvider = "noSingletonFilter")
    public void testNoSingletonFilter(final VariantContext context, final boolean pass)
            throws Exception {
        Assert.assertEquals(HaplotypeFilterLibrary.NO_SINGLETON_FILTER.test(context), pass);
    }

}
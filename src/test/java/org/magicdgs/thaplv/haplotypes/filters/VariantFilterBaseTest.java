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

import org.magicdgs.thaplv.utils.test.BaseTest;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Base test for VariantFilters
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class VariantFilterBaseTest extends BaseTest {

    protected static final Allele refA = Allele.create("A", true);
    protected static final Allele altT = Allele.create("T", false);
    protected static final Allele altC = Allele.create("C", false);

    protected static List<Genotype> getHmGenotypes(final Allele al, final int number) {
        final List<Genotype> homozygous = new ArrayList<>();
        for (int i = 0; i < number; i++) {
            homozygous.add(GenotypeBuilder
                    .create("sample_" + al + "_hm_" + +i, Arrays.asList(al, al)));
        }
        return homozygous;
    }

    protected static List<Genotype> getHzGenotypes(final Allele al1, final Allele al2,
            final int number) {
        final List<Genotype> homozygous = new ArrayList<>();
        for (int i = 0; i < number; i++) {
            homozygous.add(GenotypeBuilder
                    .create("sample_" + al1 + al2 + "_hz_" + +i, Arrays.asList(al1, al2)));
        }
        return homozygous;
    }
}

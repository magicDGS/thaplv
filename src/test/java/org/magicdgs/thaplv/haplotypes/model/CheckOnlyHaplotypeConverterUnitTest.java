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

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class CheckOnlyHaplotypeConverterUnitTest extends HaplotypeConverterUnitTest {

    @Test
    public void testFirstExceptionApply() throws Exception {
        final Genotype haploid = GenotypeBuilder.create("haploid", Arrays.asList(refA, refA));
        final Genotype diploid = GenotypeBuilder.create("diploid", Arrays.asList(refA, altC));
        final CheckOnlyHaplotypeConverter singleton = CheckOnlyHaplotypeConverter.getSingleton();
        final VariantContextBuilder builder =
                new VariantContextBuilder("test", "1", 1, 1, Arrays.asList(refA, altC));
        singleton.counter.set(0);
        // the first variant should be checked and throw an error
        try {
            singleton.apply(builder.genotypes(diploid).make());
            Assert.fail("expected UserException.BadArgumentValue");
        } catch (UserException.BadArgumentValue e) {

        }
        // the second should not be checked
        singleton.apply(builder.genotypes(diploid).make());
        // if we set to 0, the haploid should not complain
        singleton.counter.set(0);
        singleton.apply(builder.genotypes(haploid).make());
        // if we set to the frequency, the diploid should throw the error
        // the first variant should be checked and throw an error
        try {
            singleton.counter.set(CheckOnlyHaplotypeConverter.samplingFrequency);
            singleton.apply(builder.genotypes(diploid).make());
            Assert.fail("expected UserException.BadArgumentValue");
        } catch (UserException.BadArgumentValue e) {

        }
        // but not the haploid
        singleton.counter.set(CheckOnlyHaplotypeConverter.samplingFrequency);
        singleton.apply(builder.genotypes(haploid).make());
    }

}
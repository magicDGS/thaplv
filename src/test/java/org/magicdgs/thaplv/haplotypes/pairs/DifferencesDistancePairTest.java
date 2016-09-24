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

package org.magicdgs.thaplv.haplotypes.pairs;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Collections;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class DifferencesDistancePairTest {

    @Test
    public void testEmptyDistance() throws Exception {
        final DifferencesDistancePair pair = new DifferencesDistancePair("sample1", "sample2");
        Assert.assertEquals(pair.getNumberOfDifferences(), 0);
        Assert.assertEquals(pair.getNumberOfSites(), 0);
        Assert.assertEquals(pair.getNumberOfMissing(), 0);
        Assert.assertEquals(pair.getDistance(), Double.NaN);
        Assert.assertEquals(pair.getDifferencesPerSite(), Double.NaN);
    }

    @Test
    public void testComputedDistance() throws Exception {
        // create alleles
        final Allele refAllele = Allele.create((byte) 'A', true);
        final Allele altAllele = Allele.create((byte) 'C', false);
        // create three different genotypes for sample1 (sample2 is always reference)
        final Genotype sample1missing = GenotypeBuilder.createMissing("sample1", 2);
        final Genotype sample1ref = GenotypeBuilder
                .create("sample1", Collections.nCopies(2, refAllele));
        final Genotype sample1alt = GenotypeBuilder
                .create("sample1", Collections.nCopies(2, altAllele));
        final Genotype sample2ref = GenotypeBuilder
                .create("sample2", Collections.nCopies(2, refAllele));
        // create the pair
        final DifferencesDistancePair pair = new DifferencesDistancePair("sample1", "sample2");
        // create differences and missing
        int differences = 0, equals = 0, missing = 0;
        // add differences (1/2 are differences, 1/6 missing and the rest similarities)
        for (int i = 0; i < 567; i++) {
            if (i % 2 == 0) {
                // create a difference
                pair.add(sample1alt, sample2ref);
                differences++;
            } else if (i % 3 == 0) {
                // create a missing
                pair.add(sample1missing, sample2ref);
                missing++;
            } else {
                // create a similarity
                pair.add(sample1ref, sample2ref);
                equals++;
            }
            // check in every iteration
            Assert.assertEquals(pair.getNumberOfDifferences(), differences);
            Assert.assertEquals(pair.getNumberOfSites(), differences + equals);
            Assert.assertEquals(pair.getNumberOfMissing(), missing);
            Assert.assertEquals(pair.getDistance(), differences + 0.5 * missing);
            Assert.assertEquals(pair.getDifferencesPerSite(),
                    (double) differences / (differences + equals));
        }
        Assert.assertEquals(pair.getNumberOfDifferences(), differences);
        Assert.assertEquals(pair.getNumberOfSites(), differences + equals);
        Assert.assertEquals(pair.getNumberOfMissing(), missing);
        Assert.assertEquals(pair.getDistance(), differences + 0.5 * missing);
        Assert.assertEquals(pair.getDifferencesPerSite(),
                (double) differences / (differences + equals));
    }

    @Test
    public void testComputedDistanceAgainstRef() throws Exception {
        // create alleles
        final Allele refAllele = Allele.create((byte) 'A', true);
        final Allele altAllele = Allele.create((byte) 'C', false);
        // create three different genotypes for sample1 (sample2 is always reference)
        final Genotype sample1missing = GenotypeBuilder.createMissing("sample1", 2);
        final Genotype sample1ref = GenotypeBuilder
                .create("sample1", Collections.nCopies(2, refAllele));
        final Genotype sample1alt = GenotypeBuilder
                .create("sample1", Collections.nCopies(2, altAllele));
        // create the pair
        final DifferencesDistancePair pair = new DifferencesDistancePair("Reference", "sample1");
        // create differences and missing
        int differences = 0;
        int equals = 0;
        int missing = 0;
        // add differences (1/2 are differences, 1/6 missing and the rest similarities)
        for (int i = 0; i < 567; i++) {
            if (i % 2 == 0) {
                // create a difference
                pair.addReference(sample1alt);
                differences++;
            } else if (i % 3 == 0) {
                // create a missing
                pair.addReference(sample1missing);
                missing++;
            } else {
                // create a similarity
                pair.addReference(sample1ref);
                equals++;
            }
            // check in every iteration
            Assert.assertEquals(pair.getNumberOfDifferences(), differences);
            Assert.assertEquals(pair.getNumberOfSites(), differences + equals);
            Assert.assertEquals(pair.getNumberOfMissing(), missing);
            Assert.assertEquals(pair.getDistance(), differences + 0.5 * missing);
            Assert.assertEquals(pair.getDifferencesPerSite(),
                    (double) differences / (differences + equals));
        }
        Assert.assertEquals(pair.getNumberOfDifferences(), differences);
        Assert.assertEquals(pair.getNumberOfSites(), differences + equals);
        Assert.assertEquals(pair.getNumberOfMissing(), missing);
        Assert.assertEquals(pair.getDistance(), differences + 0.5 * missing);
        Assert.assertEquals(pair.getDifferencesPerSite(),
                (double) differences / (differences + equals));
    }

     @DataProvider(name = "wrongGenotypes")
    public Object[][] wrongAddGenotype() throws Exception {
         final Genotype sample1 =  GenotypeBuilder.createMissing("sample1", 1);
         final Genotype sample2 =  GenotypeBuilder.createMissing("sample2", 1);
         final Genotype wrongSample = GenotypeBuilder.createMissing("sample10", 1);
         final DifferencesDistancePair pair = new DifferencesDistancePair(sample1.getSampleName(), sample2.getSampleName());
         return new Object[][] {
                 {pair, wrongSample, sample1},
                 {pair, wrongSample, sample2}
         };
     }

    @Test(dataProvider = "wrongGenotypes", expectedExceptions = IllegalArgumentException.class)
    public void testWrongAddFirstSample(final DifferencesDistancePair pair, final Genotype wrong, final Genotype good) throws Exception {
        pair.add(wrong, good);
    }

    @Test(dataProvider = "wrongGenotypes", expectedExceptions = IllegalArgumentException.class)
    public void testWrongAddSecondSample(final DifferencesDistancePair pair, final Genotype wrong, final Genotype good) throws Exception {
        pair.add(good, wrong);
    }

    @Test(dataProvider = "wrongGenotypes", expectedExceptions = IllegalArgumentException.class)
    public void testWrongAddReference(final DifferencesDistancePair pair, final Genotype wrong, final Genotype good) throws Exception {
        pair.addReference(wrong);
    }

}
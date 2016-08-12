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

package org.magicdgs.thaplv.utils;

import org.magicdgs.thaplv.utils.test.BaseTest;

import htsjdk.variant.variantcontext.Allele;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.stream.Stream;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class AlleleUtilsUnitTest extends BaseTest {

    private static final Allele refA = Allele.create("A", true);
    private static final Allele altT = Allele.create("T", false);
    private static final Allele altC = Allele.create("C", false);

    @DataProvider(name = "allelesToCount")
    public Object[][] getAlleleCountData() {
        return new Object[][] {
                // only one allele
                {Stream.of(refA, refA, refA), Collections.singletonMap(refA, 3), true},
                {Stream.of(altT, altT), Collections.singletonMap(altT, 2), true},
                {Stream.of(altC, altC, altC, altC), Collections.singletonMap(altC, 4), true},
                {Stream.of(refA, refA, refA, refA), Collections.singletonMap(refA, 4), false},
                {Stream.of(altT, altT, altT), Collections.singletonMap(altT, 3), false},
                {Stream.of(altC, altC), Collections.singletonMap(altC, 2), false},
                // with NO_CALL
                {Stream.of(refA, refA, Allele.NO_CALL), new HashMap<Allele, Integer>(2) {{
                    put(refA, 2);
                    put(Allele.NO_CALL, 1);
                }}, true},
                {Stream.of(altT, Allele.NO_CALL), new HashMap<Allele, Integer>(2) {{
                    put(altT, 1);
                    put(Allele.NO_CALL, 1);
                }}, true},
                {Stream.of(refA, Allele.NO_CALL), Collections.singletonMap(refA, 1), false},
                {Stream.of(altC, altC, Allele.NO_CALL), Collections.singletonMap(altC, 2), false},
                // several alleles
                {Stream.of(refA, refA, refA, altT, altT, altC, Allele.NO_CALL),
                        new HashMap<Allele, Integer>(4) {{
                            put(refA, 3);
                            put(altT, 2);
                            put(altC, 1);
                            put(Allele.NO_CALL, 1);
                        }}, true},
                {Stream.of(refA, refA, refA, altT, altT, altC, Allele.NO_CALL),
                        new HashMap<Allele, Integer>(3) {{
                            put(refA, 3);
                            put(altT, 2);
                            put(altC, 1);
                        }}, false}
        };
    }

    @Test(dataProvider = "allelesToCount")
    public void testGetAlleleCounts(final Stream<Allele> alleles,
            final Map<Allele, Integer> expected, final boolean countNoCall) throws Exception {
        Assert.assertEquals(AlleleUtils.getAlleleCounts(alleles, countNoCall), expected);

    }

    @DataProvider(name = "allelesToFreqs")
    public Object[][] getAlleleFreqsData() {
        return new Object[][] {
                // only one allele
                {Stream.of(refA, refA, refA), Collections.singletonMap(refA, 1.0), true},
                {Stream.of(altT, altT), Collections.singletonMap(altT, 1.0), true},
                {Stream.of(altC, altC, altC, altC), Collections.singletonMap(altC, 1.0), true},
                {Stream.of(refA, refA, refA, refA), Collections.singletonMap(refA, 1.0), false},
                {Stream.of(altT, altT, altT), Collections.singletonMap(altT, 1.0), false},
                {Stream.of(altC, altC), Collections.singletonMap(altC, 1.0), false},
                // with NO_CALL
                {Stream.of(refA, refA, Allele.NO_CALL), new HashMap<Allele, Double>(2) {{
                    put(refA, 2d / 3);
                    put(Allele.NO_CALL, 1d / 3);
                }}, true},
                {Stream.of(altT, Allele.NO_CALL), new HashMap<Allele, Double>(2) {{
                    put(altT, 0.5);
                    put(Allele.NO_CALL, 0.5);
                }}, true},
                {Stream.of(refA, Allele.NO_CALL), Collections.singletonMap(refA, 1.0), false},
                {Stream.of(altC, altC, Allele.NO_CALL), Collections.singletonMap(altC, 1.0), false},
                // several alleles
                {Stream.of(refA, refA, refA, altT, altT, altC, Allele.NO_CALL),
                        new HashMap<Allele, Double>(4) {{
                            put(refA, 3d / 7);
                            put(altT, 2d / 7);
                            put(altC, 1d / 7);
                            put(Allele.NO_CALL, 1d / 7);
                        }}, true},
                {Stream.of(refA, refA, refA, altT, altT, altC, Allele.NO_CALL),
                        new HashMap<Allele, Double>(3) {{
                            put(refA, 3d / 6);
                            put(altT, 2d / 6);
                            put(altC, 1d / 6);
                        }}, false}
        };
    }

    @Test(dataProvider = "allelesToFreqs")
    public void testGetAlleleFrequencies(final Stream<Allele> alleles,
            final Map<Allele, Double> expectedFrequencies, final boolean countNoCall)
            throws Exception {
        Assert.assertEquals(AlleleUtils.getAlleleFrequencies(alleles, countNoCall),
                expectedFrequencies);

    }

}
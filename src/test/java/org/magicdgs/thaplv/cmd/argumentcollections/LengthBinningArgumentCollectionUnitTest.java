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

package org.magicdgs.thaplv.cmd.argumentcollections;

import org.magicdgs.thaplv.utils.test.BaseTest;

import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class LengthBinningArgumentCollectionUnitTest extends BaseTest {


    @DataProvider(name = "badArguments")
    public Object[][] getBadArguments() {
        return new Object[][] {
                {-1, 1, 1},
                {2, 1, 1},
                {0, 1, -1}
        };
    }

    @Test(dataProvider = "badArguments", expectedExceptions = IllegalArgumentException.class)
    public void testBadArguments(final int min, final Integer max, final int binDistance) {
        new LengthBinningArgumentCollection(min, max, binDistance);
    }

    @DataProvider(name = "correctArguments")
    public Object[][] getBinningArguments() {
        return new Object[][] {
                {0, 10, 1},
                {10, null, 10},
                {0, 10, 10},
                {10, null, 2}
        };
    }

    @Test(dataProvider = "correctArguments")
    public void testBinningArguments(final int min, final Integer max, final int binDistance)
            throws Exception {
        new LengthBinningArgumentCollection(min, max, binDistance);
    }

}
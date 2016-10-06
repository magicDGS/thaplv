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

import org.broadinstitute.hellbender.exceptions.UserException;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Optional;
import java.util.concurrent.ThreadPoolExecutor;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class MultiThreadComputationArgumentCollectionUnitTest extends BaseTest {

    @DataProvider(name = "badArguments")
    public Object[][] getBadArguments() {
        return new Object[][] {{-1}, {-50}};
    }

    @Test(dataProvider = "badArguments", expectedExceptions = IllegalArgumentException.class)
    public void testBadArguments(final int bufferSize) {
        new MultiThreadComputationArgumentCollection(bufferSize);
    }

    @DataProvider(name = "smallBufferForThreads")
    public Object[][] getSmallBuffer() {
        return new Object[][] {
                {1, 2},
                {2, 3}
        };
    }

    @Test(dataProvider = "smallBufferForThreads")
    public void testSmallBuffer(final int bufferSize, final int threads)
            throws Exception {
        final MultiThreadComputationArgumentCollection args =
                new MultiThreadComputationArgumentCollection(bufferSize);
        Assert.assertFalse(args.getFixedThreadPoolWithBoundedQueue().isPresent());
        args.nThreads = threads;
        final Optional<ThreadPoolExecutor> executor = args.getFixedThreadPoolWithBoundedQueue();
        Assert.assertTrue(executor.isPresent());
        Assert.assertEquals(executor.get().getCorePoolSize(), threads);
        Assert.assertEquals(executor.get().getQueue().remainingCapacity(), threads);
    }

    @DataProvider(name = "hugeBufferForThreads")
    public Object[][] getHugeBuffer() {
        return new Object[][] {
                {100, 2},
                {200, 3}
        };
    }

    @Test(dataProvider = "hugeBufferForThreads")
    public void testHugeBuffer(final int bufferSize, final int threads)
            throws Exception {
        final MultiThreadComputationArgumentCollection args =
                new MultiThreadComputationArgumentCollection(bufferSize);
        Assert.assertFalse(args.getFixedThreadPoolWithBoundedQueue().isPresent());
        args.nThreads = threads;
        final Optional<ThreadPoolExecutor> executor = args.getFixedThreadPoolWithBoundedQueue();
        Assert.assertTrue(executor.isPresent());
        Assert.assertEquals(executor.get().getCorePoolSize(), threads);
        Assert.assertEquals(executor.get().getQueue().remainingCapacity(), bufferSize);
    }

    @DataProvider
    public Object[][] badThreads() {
        return new Object[][] {{-1}, {0}};
    }

    @Test(dataProvider = "badThreads", expectedExceptions = UserException.BadArgumentValue.class)
    public void testBadUserArgument(final int badThreads) throws Exception {
        final MultiThreadComputationArgumentCollection args =
                new MultiThreadComputationArgumentCollection();
        args.nThreads = badThreads;
        args.getFixedThreadPoolWithBoundedQueue();
    }


}
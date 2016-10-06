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

import org.magicdgs.thaplv.utils.concurrent.ExecutorsFactory;

import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollectionDefinition;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Optional;
import java.util.concurrent.ThreadPoolExecutor;

/**
 * Argument collection for multi-thread computation. It contains factory methods to create thread
 * pools and executors.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class MultiThreadComputationArgumentCollection implements ArgumentCollectionDefinition {
    private static final long serialVersionUID = 1L;

    /** Number of threads to use in the computation. Default value is not using multi-thread. */
    @Argument(fullName = "threads", shortName = "nt", doc = "Number of threads to use in the computation.", optional = true)
    public int nThreads = 1;

    // buffer size for the queue
    private final int defaultBufferSize;

    /** Creates a multi-thread argument collection with buffer-size for blocking queues. */
    public MultiThreadComputationArgumentCollection(final int defaultBufferSize) {
        Utils.validateArg(defaultBufferSize > 0, () -> "invalid buffer-size: " + defaultBufferSize);
        this.defaultBufferSize = defaultBufferSize;
    }

    /** Creates a multi-thread argument collection with default buffer-size for blocking queues. */
    public MultiThreadComputationArgumentCollection() {
        this(500_000);
    }

    /** @throws UserException.BadArgumentValue if they are incorrect */
    private void validateUserArgs() {
        if (nThreads <= 0) {
            throw new UserException.BadArgumentValue("threads", String.valueOf(nThreads),
                    "threads should be positive");
        }
    }

    /**
     * Creates a new fixed thread pool with bounded queue using user-provided arguments.
     *
     * @return empty if {@code nThreads == 1}; optional with executor otherwise
     */
    public Optional<ThreadPoolExecutor> getFixedThreadPoolWithBoundedQueue() {
        validateUserArgs();
        if (nThreads == 1) {
            return Optional.empty();
        }
        // don't waste resources when generating threads
        // if the user provides more than the buffer size it's his own risk
        final int bufferSize = Math.max(defaultBufferSize, nThreads);
        return Optional.of((ThreadPoolExecutor) ExecutorsFactory
                .newFixedThreadPoolWithBoundedQueue(nThreads, bufferSize));
    }

}

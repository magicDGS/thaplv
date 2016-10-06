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

package org.magicdgs.thaplv.utils.concurrent;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

/**
 * Factory class for create executors.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class ExecutorsFactory {

    /** Cannot be instantiated. */
    private ExecutorsFactory() {}

    /**
     * Creates a fixed thread pool {@link java.util.concurrent.Executors#newFixedThreadPool}
     * with {@code nThreads} and a bounded working queue where rejected jobs will be added
     * once the buffer size is not full.
     *
     * @param nThreads   the number of threads in the pool
     * @param bufferSize the buffer size for the working queue
     *
     * @return the newly created thread pool
     *
     * @throws IllegalArgumentException for invalid arguments
     */
    // based on http://stackoverflow.com/questions/11568821/how-not-to-overwhelm-java-executorservice-task-queue
    public static ExecutorService newFixedThreadPoolWithBoundedQueue(int nThreads, int bufferSize) {
        Utils.validateArg(nThreads > 0, () -> "invalid number of threads: " + nThreads);
        Utils.validateArg(bufferSize >= 1, () -> "invalid buffer size: " + bufferSize);
        // this is like Executors.newFixedThreadPool() but with an ArrayBlockingQueue queue instead of
        // LinkedBlockingQueue and with a rejected execution handler that put any job that is rejected
        // in the working queue for the executor.
        return new ThreadPoolExecutor(nThreads, nThreads, 0L, TimeUnit.MILLISECONDS,
                new ArrayBlockingQueue<>(bufferSize), Executors.defaultThreadFactory(),
                (r, executor) -> {
                    try {
                        // this will block if the queue is full as opposed to throwing
                        executor.getQueue().put(r);
                    } catch (InterruptedException e) {
                        new InterruptedException(e.getMessage());
                    }
                });
    }

}

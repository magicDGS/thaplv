/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2015 Daniel Gómez-Sánchez
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

package org.magicdgs.thaplv.utils.concurrent;

import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.RejectedExecutionHandler;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

/**
 * Class to provide several executors
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class OtherExecutors {

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
     * @throws IllegalArgumentException if {@code nThreads <= 0}
     */
    public static ExecutorService newFixedThreadPoolWithBoundedQueue(int nThreads, int bufferSize) {
        final BlockingQueue<Runnable> queue = new ArrayBlockingQueue<Runnable>(bufferSize);
        return new ThreadPoolExecutor(nThreads, nThreads, 0L, TimeUnit.MILLISECONDS, queue,
                Executors
                        .defaultThreadFactory(), new RejectedToQueueHandler());
    }

    /**
     * This handler put any job that is rejected in the working queue for the executor.
     * The main use of this class should be for use with a bounded queue.
     *
     * Based on http://stackoverflow.com/questions/11568821/how-not-to-overwhelm-java-executorservice-task-queue
     */
    private static class RejectedToQueueHandler implements RejectedExecutionHandler {
        @Override
        public void rejectedExecution(Runnable r, ThreadPoolExecutor executor) {
            // this will block if the queue is full as opposed to throwing
            try {
                executor.getQueue().put(r);
            } catch (InterruptedException e) {
                throw new RuntimeException(e);
            }
        }
    }

}

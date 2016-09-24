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

import org.magicdgs.thaplv.cmd.ThaplvArgumentDefinitions;

import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollectionDefinition;
import org.broadinstitute.hellbender.exceptions.UserException;

/**
 * Argument collection for encapsulate sliding window parameters.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class SlidingWindowArgumentCollection implements ArgumentCollectionDefinition {
    private static final long serialVersionUID = 1L;

    /** The window size. */
    @Argument(fullName = ThaplvArgumentDefinitions.WINDOW_SIZE_LONG, shortName = ThaplvArgumentDefinitions.WINDOW_SIZE_SHORT, doc = "Window size for compute IBD regions (in bp)", optional = true)
    public int windowSize;

    /** The step size. */
    @Argument(fullName = ThaplvArgumentDefinitions.WINDOW_STEP_LONG, shortName = ThaplvArgumentDefinitions.WINDOW_STEP_SHORT, doc = "Step size for compute IBD regions (in bp)", optional = true)
    public int stepSize;

    /**
     * Constructor for default window and step size.
     */
    public SlidingWindowArgumentCollection(final int defaultWindowSize, final int defaultStepSize) {
        this.windowSize = defaultWindowSize;
        this.stepSize = defaultStepSize;
        try {
            validateArguments();
        } catch (UserException e) {
            throw new IllegalArgumentException("Bad default value: " + e.getMessage());
        }
    }

    /**
     * Validate the arguments. Should be called before using {@link #windowSize} or
     * {@link #stepSize} to avoid user bad input.
     */
    public void validateArguments() {
        if (windowSize < 1) {
            throw new UserException.BadArgumentValue(
                    "--" + ThaplvArgumentDefinitions.WINDOW_SIZE_LONG, String.valueOf(windowSize),
                    "should be a positive integer");
        } else if (stepSize < 1) {
            throw new UserException.BadArgumentValue(
                    "--" + ThaplvArgumentDefinitions.WINDOW_STEP_LONG, String.valueOf(stepSize),
                    "should be a positive integer");
        } else if (stepSize > windowSize) {
            throw new UserException.BadArgumentValue(
                    "--" + ThaplvArgumentDefinitions.WINDOW_STEP_LONG, String.valueOf(stepSize),
                    "should be smaller or equals to --" + ThaplvArgumentDefinitions.WINDOW_SIZE_LONG
            );
        }
    }
}

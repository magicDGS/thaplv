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

import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollectionDefinition;
import org.broadinstitute.hellbender.exceptions.UserException;

/**
 * Argument collection for tools that need binning statistics by distance between variants.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class LengthBinningArgumentCollection implements ArgumentCollectionDefinition {
    private static final long serialVersionUID = 1L;

    /**
     * Minimum distance between SNPs (inclusive). It should be an integer between 0 and
     * {@link Integer#MAX_VALUE}.
     */
    @Argument(fullName = "minimum-distance", shortName = "mindist", doc = "Minimum distance (in bp, inclusive) between SNPs to perform the computation. Set to 0 to do not have a threshold for minimum distance.", optional = true)
    public int min;

    /**
     * Maximum distance between SNPs (inclusive). It should be an integer between 1 and
     * {@link Integer#MAX_VALUE}. If set to {@code nuill}, it won't have a threshold.
     */
    @Argument(fullName = "maximum-distance", shortName = "maxdist", doc = "Maximum distance (in bp, inclusive) between SNPs to perform the computation. Set to null to do not have a threshold for maximum distance.", optional = true)
    public Integer max;

    /** The distance between SNPs to bin them together. */
    @Argument(fullName = "bin-distance", shortName = "bindist", doc = "Distance (in bp) between pairs to be considered in the same bin.", optional = true)
    public int binDistance;

    /**
     * New argument collection with set default values.
     *
     * @param min         default minimum-distance.
     * @param max         default maximum-distance.
     * @param binDistance default bin-distance.
     */
    public LengthBinningArgumentCollection(final int min, final Integer max,
            final int binDistance) {
        this.min = min;
        this.max = max;
        this.binDistance = binDistance;
        try {
            validateArgs();
        } catch (UserException e) {
            // rethrow as illegal argument exception
            throw new IllegalArgumentException(e.getMessage(), e);
        }
    }


    /**
     * Log warnings and info about the parameters.
     *
     * @param logger the logger to use.
     */
    public void logWarnings(final Logger logger) {
        if (max == null && min == 0) {
            logger.warn(
                    "All variants within the chromosome will be used. This could cause a memory overload.");
        } else {
            logger.info(
                    "Only variants with a distance between them in the range [{},{}] bp will be considered.",
                    min, (max == null) ? "chromosome lenght" : max);
        }
    }

    /**
     * Validate the arguments provided by the user.
     *
     * @throws UserException.BadArgumentValue if the arguments are not correct
     */
    public void validateArgs() {
        if (max != null && max < 1) {
            throw new UserException.BadArgumentValue("maximum-distance", String.valueOf(max),
                    "should be positive");
        } else if (max != null && min > max) {
            throw new UserException.BadArgumentValue("minimum-distance", String.valueOf(min),
                    String.format("should be smaller than 'maximum-distance' (%s)", max));
        } else if (min < 0) {
            throw new UserException.BadArgumentValue("minimum-distance", String.valueOf(min),
                    "should be positive or 0");
        }
        if (binDistance < 1) {
            throw new UserException.BadArgumentValue("bin-distance", String.valueOf(binDistance),
                    "should be positive");
        }
    }

    /**
     * Returns {@code true} if distance between two positions are bellow {@link #max};
     * {@code false} otherwise.
     */
    public boolean bellowMaximumDistance(final int position1, final int position2) {
        return (max == null || (position1 <= position2 + max));
    }

    /**
     * Returns {@code true} if distance between two positions exceed {@link #min};
     * {@code false} otherwise.
     */
    public boolean exceedMinimumDistance(final int position1, final int position2) {
        return min == 0 || Math.abs(position1 - position2) >= min;
    }

}

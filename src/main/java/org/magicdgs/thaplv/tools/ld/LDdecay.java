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

package org.magicdgs.thaplv.tools.ld;

import org.magicdgs.thaplv.cmd.ThaplvArgumentDefinitions;
import org.magicdgs.thaplv.cmd.programgroups.AlphaProgramGroup;
import org.magicdgs.thaplv.haplotypes.filters.HaplotypeFilterLibrary;
import org.magicdgs.thaplv.haplotypes.filters.NumberOfMissingFilter;
import org.magicdgs.thaplv.tools.ld.engine.QueueLD;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.HaploidWalker;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.VariantFilter;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.IOException;

/**
 * Class to compute LD in bins in a concurrent way
 * (1.5 million variants)
 * TODO: implement minor allele frequency
 * TODO: unit tests
 * TODO: better documentation
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
@CommandLineProgramProperties(oneLineSummary = "Compute and bin different linkage disequilibrium statistics",
        // TODO: better description
        summary = "Compute and bin different linkage disequilibrium statistics. Median and quantiles are an approximation",
        programGroup = AlphaProgramGroup.class)
public final class LDdecay extends HaploidWalker {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output prefix for LD results.", optional = false)
    public String outputPrefix;

    @Argument(fullName = "minimum-distance", shortName = "mindist", doc = "Minimum distance (in bp, inclusive) between SNPs to compute LD.", optional = true)
    public int min = 0;

    // TODO: is this inclusive or exclusive???
    // TODO: change this to set to null for compute all the chromosome
    @Argument(fullName = "maximum-distance", shortName = "maxdist", doc = "Maximum distance (in bp) between SNPs to compute LD. If set to < 1, compute all the SNPs within the chromosome for minimum distance.", optional = true)
    public int max = 10_000;

    @Argument(fullName = "bin-distance", shortName = "bindist", doc = "Distance (in bp) between pairs to be considered in the same bin.", optional = true)
    public int binDistance = 1_000;

    @Argument(fullName = "chi-square", shortName = "chi", doc = "Chi-square quantile to assess the significance of max. correlation and compute LD.", optional = true)
    public double chiSqrQuantile = 0.95;

    // TODO: change this to set to null for use only no-missing data
    // TODO: change the name
    @Argument(fullName = "minimum-samples", shortName = "mins", doc = "Minimum number of samples available to compute LD. Setting to a negative number only use no-missing data.", optional = true)
    public int minSamples = -1; // TODO: change default value

    // TODO: change logic, always remove them except requested
    @Argument(fullName = "remove-singletons", shortName = "rms", doc = "Remove pairs where the minor variant is a singleton.", optional = true)
    public boolean rmSingletons = false;

    // TODO: add option for minimum allele frequency
    // @Argument(fullName = "minimum-allele-frequency", shortName = "maf", doc = "Minimum allele frequency for a marker to compute LD.", optional = true)
    // public Double maf = null;

    // TODO: add option for only output binned results
    // TODO: probably much better the other way around -> option for output all pairs
    // @Argument(fullName = "only-bin", shortName = "bin", doc = "Output only binned results.", optional = true)
    // public boolean onyBin = false;

    @Argument(fullName = ThaplvArgumentDefinitions.NUMBER_OF_THREAD_LONG, shortName = ThaplvArgumentDefinitions.NUMBER_OF_THREAD_SHORT, doc = "Number of threads to use in the computation.", optional = true)
    public int nThreads = 1;

    @Override
    protected boolean requiresOutputPloidy() {
        return false;
    }

    @Override
    protected boolean allowsDontCheck() {
        return true;
    }

    // this is the queue for store the ld results
    private QueueLD queue;


    @Override
    protected VariantFilter makeVariantFilter() {
        // this is already removing invariant sites
        VariantFilter filter = HaplotypeFilterLibrary.BIALLELIC_FILTER;
        if (rmSingletons) { // TODO: this should change with the parameter
            filter = filter.and(HaplotypeFilterLibrary.NO_SINGLETON_FILTER);
        }
        if (minSamples == -1) { // TODO: this should change with the parameter
            // TODO: I don't know if this is working as expected
            filter = filter.and(new NumberOfMissingFilter(minSamples));
        }
        // TODO: the maf threshold should be here
        return filter;
    }

    // TODO: change messages and validation of args
    @Override
    public void onTraversalStart() {
        if (chiSqrQuantile < 0 || chiSqrQuantile >= 1) {
            throw new UserException.BadArgumentValue("Chi-square quantile should be in the range (0, 1)");
        }
        if (nThreads < 1) {
            throw new UserException.BadArgumentValue("Threads cannot be less than 1");
        }
        if (minSamples < 1) {
            throw new UserException.BadArgumentValue("Minimum samples cannot be smaller than 1");
        }
        if (binDistance < 1) {
            throw new UserException.BadArgumentValue("Bin length cannot be smaller than 1");
        }
        if (max > 0 && min > max) {
            throw new UserException.BadArgumentValue(
                    "Minimum distance between markers cannot be bigger than maximum distance");
        } else if (min < 0) {
            throw new UserException.BadArgumentValue(
                    "Minimum distance between markers cannot be negative");
        }

        // warnings
        if (chiSqrQuantile == 0) {
            logger.warn(
                    "Computation for r2 will be performed for all maximum r2 independently of its value (except for 0) because the chi-square quantile is set to 0");
        }
        if (max < 1 && min == 0) {
            logger.info("All variants within the chromosome will be used to compute LD");
        } else {
            logger.info(
                    "Only variants with a distance between them in the range [{},{}] bp will be considered.",
                    min, (max < 1) ? "chromosome lenght" : max);
        }


        queue = new QueueLD(outputPrefix, binDistance, min, max, minSamples,
                rmSingletons, chiSqrQuantile, nThreads, logger);
    }

    @Override
    public void apply(VariantContext variant, ReadsContext readsContext,
            ReferenceContext referenceContext, FeatureContext featureContext) {
        try {
            // TODO: debugging the variants in RAM was done before here
            // TODO: implement it inside the queue will be better
            queue.add(variant);
        } catch (IOException e) {
            // TODO: change exception handling
            throw new UserException(e.getMessage(), e);
        }
    }
}

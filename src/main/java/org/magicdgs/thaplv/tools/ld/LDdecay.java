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

import org.magicdgs.thaplv.cmd.argumentcollections.LengthBinningArgumentCollection;
import org.magicdgs.thaplv.cmd.argumentcollections.MultiThreadComputationArgumentCollection;
import org.magicdgs.thaplv.cmd.programgroups.AlphaProgramGroup;
import org.magicdgs.thaplv.haplotypes.filters.HaplotypeFilterLibrary;
import org.magicdgs.thaplv.haplotypes.filters.NumberOfMissingFilter;
import org.magicdgs.thaplv.tools.ld.engine.LDdecayOutput;
import org.magicdgs.thaplv.tools.ld.engine.QueueLD;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.HaploidWalker;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.VariantFilter;
import org.broadinstitute.hellbender.exceptions.UserException;

/**
 * Computes linkage disequilibrium statistics (based on Pearson's correlation r<sup>2</sup>)
 * binning
 * on the fly to avoid storing results. In addition, the computation is filtered by the maximum
 * r<sup>2</sup> that can be reached between pairs.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
@CommandLineProgramProperties(oneLineSummary = "Computes binned linkage disequilibrium statistics.",
        summary =
                "Computes and bin different linkage disequilibrium statistics based on Pearson's correlation. "
                        + "It computes average and standard deviation for each bin, and approximation for the median and quantiles. "
                        + "The statistics will be only computed if the maximum correlation could reach a significant value based on the chi-square distribution. "
                        + "In addition, only biallelic SNPs that pass the imposed filters (after removing missing data) wil be used.",
        programGroup = AlphaProgramGroup.class)
public final class LDdecay extends HaploidWalker {

    @VisibleForTesting
    static final String CHI_SQR_QUANTILE_ARGNAME = "chi-square";
    @VisibleForTesting
    static final String MINIMUM_SAMPLES_ARGNAME = "minimum-samples";
    @VisibleForTesting
    static final String INCLUDE_SINGLETONS_ARGNAME = "include-singletons";


    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output prefix for LD results.", optional = false)
    public String outputPrefix;

    @ArgumentCollection(doc = "Binning parameters")
    public LengthBinningArgumentCollection lengthBinningArgumentCollection =
            new LengthBinningArgumentCollection(0, 10_000, 1_000);

    @ArgumentCollection
    public MultiThreadComputationArgumentCollection multiThreadArgumentCollection =
            new MultiThreadComputationArgumentCollection();

    @Argument(fullName = CHI_SQR_QUANTILE_ARGNAME, doc = "Chi-square quantile to assess the significance of max. correlation and compute LD.", optional = true)
    public double chiSqrQuantile = 0.95;

    @Argument(fullName = MINIMUM_SAMPLES_ARGNAME, doc = "Minimum number of samples available to compute LD. Setting to null only use no-missing data.", optional = true)
    public Integer minSamples = null;

    @Argument(fullName = INCLUDE_SINGLETONS_ARGNAME, doc = "Include pairs where the minor variant is a singleton if it fullfit the requirements.", optional = true)
    public boolean includeSingletons = false;

    @Override
    protected boolean requiresOutputPloidy() {
        return false;
    }

    @Override
    protected boolean allowsCheckOnly() {
        return true;
    }

    // this is the queue for store the ld results
    private QueueLD queue;

    // number of variants processed
    private int nVariants = 0;

    // this filter is important to avoid computation/loading in the queue the variants
    @Override
    protected VariantFilter makeVariantFilter() {
        // this is already removing invariant sites
        VariantFilter filter = HaplotypeFilterLibrary.BIALLELIC_FILTER;
        if (includeSingletons) {
            logger.info("Singletons will be included.");
            logger.warn("Including singletons in the analysis may lead to spurious results.");
        } else {
            logger.info("Singletons will be excluded.");
            filter = filter.and(HaplotypeFilterLibrary.NO_SINGLETON_FILTER);
        }
        logger.info("Variants with less than {} samples with missing genotypes will be excluded.",
                minSamples);
        filter = filter.and(new NumberOfMissingFilter(minSamples));
        return filter;
    }

    @Override
    public void onTraversalStart() {
        // check and change the samples arguments
        final int nSamples = getHeaderForVariants().getNGenotypeSamples();
        if (minSamples == null) {
            minSamples = nSamples;
        } else {
            minSamples = Math.min(minSamples, nSamples);
        }

        lengthBinningArgumentCollection.logWarnings(logger);

        // warnings
        if (chiSqrQuantile == 0) {
            logger.warn(
                    "All variant pairs will be included in the analysis because chi-square quantile is set to 0.");
        }

        // create queue
        final LDdecayOutput output =
                new LDdecayOutput(outputPrefix, lengthBinningArgumentCollection.binDistance);
        queue = new QueueLD(output, lengthBinningArgumentCollection, minSamples, !includeSingletons,
                chiSqrQuantile, multiThreadArgumentCollection);
    }

    /** Overrides to validate the arguments. */
    protected String[] customCommandLineValidation() {
        if (chiSqrQuantile < 0 || chiSqrQuantile >= 1) {
            throw new UserException.BadArgumentValue(CHI_SQR_QUANTILE_ARGNAME,
                    String.valueOf(chiSqrQuantile), "should be in the range (0, 1)");
        }
        if (minSamples < 1) {
            throw new UserException.BadArgumentValue(MINIMUM_SAMPLES_ARGNAME,
                    String.valueOf(minSamples), "should be a positive integer");
        }
        lengthBinningArgumentCollection.validateArgs();
        return super.customCommandLineValidation();
    }

    @Override
    public void apply(VariantContext variant, ReadsContext readsContext,
            ReferenceContext referenceContext, FeatureContext featureContext) {
        queue.add(variant);
        if (++nVariants % 1000 == 0) {
            logger.debug("Number of variants in RAM: {}; Number of pairs computed: {}.",
                    () -> queue.variantsInRam(), () -> queue.computedPairs());

        }
    }

    /** This return the number of pairs computed. */
    @Override
    public Object onTraversalSuccess() {
        queue.finalizeQueue();
        final int computedPairs = queue.computedPairs();
        logger.info("Computed linkage disequilibrium statistics on {} pairs out of {}.",
                () -> computedPairs, () -> queue.addedPairs());
        return computedPairs;
    }

    @Override
    public void closeTool() {
        if (queue != null) {
            queue.close();
        }
    }
}

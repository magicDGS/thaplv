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

package org.magicdgs.thaplv.tools.analyses;

import org.magicdgs.thaplv.haplotypes.AlleleVector;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.HaploidWalker;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.util.ArrayList;
import java.util.List;

/**
 * Tool for compute H12, H1 and H2 for haplotypes.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class ComputeH12 extends HaploidWalker {

    @Argument(fullName = "snpWindow", shortName = "w", doc = "Number of SNPs to consider in each non-overlapping window.", optional = true)
    public int snpWindow = 400;

    @Override
    protected boolean requiresOutputPloidy() {
        return false;
    }

    @Override
    protected boolean allowsCheckOnly() {
        return true;
    }

    // accumulated alleles
    protected List<AlleleVector> accumulator;

    @Override
    protected String[] customCommandLineValidation() {
        if (snpWindow <= 0) {
            throw new UserException.BadArgumentValue("snpWindow", String.valueOf(snpWindow),
                    "Number of SNPs in the window should be a positive integer");
        }
        accumulator = new ArrayList<>(snpWindow);
        return super.customCommandLineValidation();
    }

    /** Accumulates in */
    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext,
            final ReferenceContext referenceContext, final FeatureContext featureContext) {
        if (accumulator.size() == snpWindow) {
            apply(accumulator);
            accumulator.clear();
        }
        accumulator.add(new AlleleVector(variant));
    }

    @Override
    public Object onTraversalSuccess() {
        apply(accumulator);
        return null;
    }

    private void apply(final List<AlleleVector> variants) {

        // TODO: compute the statistic
    }

}

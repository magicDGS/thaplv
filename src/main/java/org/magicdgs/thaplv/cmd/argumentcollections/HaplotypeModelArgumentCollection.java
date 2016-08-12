/*
 * Copyright (c) 2016, Daniel Gomez-Sanchez <daniel.gomez.sanchez@hotmail All rights reserved.
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
import org.magicdgs.thaplv.haplotypes.model.HaplotypeModel;
import org.magicdgs.thaplv.haplotypes.model.VariantHaplotypeConverter;

import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollectionDefinition;
import org.broadinstitute.hellbender.exceptions.UserException;

/**
 * Argument collection for haplotype models
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public abstract class HaplotypeModelArgumentCollection implements ArgumentCollectionDefinition {
    private static final long serialVersionUID = 1L;

    private final boolean allowDontCheck;

    @Argument(fullName = ThaplvArgumentDefinitions.HAPLOTYPE_MODEL_LONG, shortName = ThaplvArgumentDefinitions.HAPLOTYPE_MODEL_SHORT, doc = "Haplotype model for convert calls.", common = false, optional = true)
    public HaplotypeModel haplotypeModel = HaplotypeModel.HAPLOID;

    /**
     * Constructor for the argument collection
     *
     * @param allowDontCheck is {@link HaplotypeModel#DONT_CHECK} model allowed?
     */
    public HaplotypeModelArgumentCollection(final boolean allowDontCheck) {
        this.allowDontCheck = allowDontCheck;
    }

    /**
     *  Get the {@link VariantHaplotypeConverter} from the tool, checking for correct parameters
     */
    public final VariantHaplotypeConverter getHaplotypeConverter() {
        if(!allowDontCheck && haplotypeModel == HaplotypeModel.DONT_CHECK) {
            throw new UserException.BadArgumentValue(ThaplvArgumentDefinitions.HAPLOTYPE_MODEL_LONG,
                    haplotypeModel.toString(),
                    "is not allowed in this tool.");
        }
        return haplotypeConverterFromArguments();
    }

    /**
     * Get the haplotype converter from the arguments
     */
    protected abstract VariantHaplotypeConverter haplotypeConverterFromArguments();

    protected void throwDontCheckException() {

    }

}

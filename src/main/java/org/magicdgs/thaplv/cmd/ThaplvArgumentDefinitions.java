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

package org.magicdgs.thaplv.cmd;

/**
 * Argument definitions for thaplv tools
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class ThaplvArgumentDefinitions {

    private ThaplvArgumentDefinitions() {}

    public final static String HAPLOTYPE_MODEL_LONG = "haplotypeModel";
    public final static String HAPLOTYPE_MODEL_SHORT = "H";

    public final static String OUTPUT_PLOIDY_LONG = "outputPloidy";
    public final static String OUTPUT_PLOIDY_SHORT = "ploidy";

    public final static String ONLY_POLYMORPHIC_LONG = "onlyPolymorphic";
    public final static String ONLY_POLYMORPHIC_SHORT = "polymorphic";

    public final static String WINDOW_SIZE_LONG = "widow-size";
    public final static String WINDOW_SIZE_SHORT = "ws";
    public final static String WINDOW_STEP_LONG = "window-step-size";
    public final static String WINDOW_STEP_SHORT = "ss";

    public final static String NUMBER_OF_THREAD_LONG = "threads";
    public final static String NUMBER_OF_THREAD_SHORT = "nt";

}

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

package org.magicdgs.thaplv;

import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.util.ArrayList;
import java.util.List;

/**
 * Main class for thaplv
 *
 * @author Daniel Gomez-Sanchez
 */
public class Main {

    /**
     * The packages we wish to include in our command line.
     */
    protected static List<String> getPackageList() {
        final List<String> packageList = new ArrayList<>();
        // TODO: add some of the bundle tools like IndexFeatureFile
        // packageList.add("org.broadinstitute.hellbender");
        packageList.add("org.magicdgs.thaplv.tools");
        return packageList;
    }

    /**
     * Entry point for thaplv.
     */
    public static void main(final String[] args) {
        try {
            Object result = new Main().instanceMain(args);
            if (result != null) {
                System.out.println("Tool returned:\n" + result);
            }
        } catch (final UserException.CommandLineException e) {
            // the GATK framework already prints the error
            System.exit(1);
        } catch (final UserException e) {
            // this prints the error for user exceptions
            CommandLineProgram.printDecoratedUserExceptionMessage(System.err, e);
            System.exit(2);
        } catch (final Exception e) {
            e.printStackTrace();
            System.exit(3);
        }
    }

    /**
     * This is a necessary hack for use {@linl CommandLineProgramTest}
     */
    public Object instanceMain(final String[] args) {
        return new org.broadinstitute.hellbender.Main()
                .instanceMain(args, getPackageList(), "thaplv");
    }
}

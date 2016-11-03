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

package org.magicdgs.thaplv;

import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Main class for thaplv
 *
 * @author Daniel Gomez-Sanchez
 */
public class Main extends org.broadinstitute.hellbender.Main {

    // TODO: Override when it is in GATK
    protected String getCommandLineName() {
        return "thaplv";
    }

    /** The packages we wish to include in our command line. */
    @Override
    protected List<String> getPackageList() {
        final List<String> packageList = new ArrayList<>();
        packageList.add("org.magicdgs.thaplv.tools");
        return packageList;
    }

    /** The single classes that we wish to include in our command line. */
    @Override
    protected List<Class<? extends CommandLineProgram>> getClassList() {
        // TODO: add some of the bundle tools like IndexFeatureFile
        return Collections.emptyList();
    }

    /** Entry point for thaplv. */
    public static void main(final String[] args) {
        new Main().mainEntry(args);
    }

    /** Entry point for this instance. */
    protected final void mainEntry(final String[] args) {
        try {
            final Object result = instanceMain(args);
            handleResult(result);
            System.exit(0);
        } catch (final UserException.CommandLineException e) {
            //the usage has already been printed so don't print it here.
            System.exit(1);
        } catch (final UserException e) {
            CommandLineProgram.printDecoratedUserExceptionMessage(System.err, e);
            System.exit(2);
        } catch (final Exception e) {
            e.printStackTrace();
            System.exit(3);
        }
    }

    // TODO: remove when in GATK or include special behaviour
    protected void handleResult(final Object result) {
        if (result != null) {
            System.out.println("Tool returned:\n" + result);
        }
    }

    /** Override to include our program name in the command line. */
    @Override
    public Object instanceMain(final String[] args) {
        // TODO: remove this method when it is in GATK
        return instanceMain(args, getPackageList(), getClassList(), getCommandLineName());
    }

}

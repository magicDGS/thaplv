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

package org.broadinstitute.hellbender;

import org.magicdgs.thaplv.utils.test.BaseTest;

import htsjdk.samtools.util.Log;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.logging.BunnyLog;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * This class is a modification of the GATK {@link CommandLineProgramTest} to use our own BaseTest
 * class and do not load statically the GATK resources (not existing), so it is hacking the
 * framework.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class CommandLineProgramTest extends BaseTest {

    /**
     * For testing support.  Given a name of a Main CommandLineProgram and it's arguments, builds
     * the arguments appropriate for calling the
     * program through Main
     *
     * @param args List<String> of command line arguments
     *
     * @return String[] of command line arguments
     */
    public String[] makeCommandLineArgs(final List<String> args) {
        return makeCommandLineArgs(args, getTestedClassName());
    }

    /**
     * For testing support.  Given a name of a Main CommandLineProgram and it's arguments, builds
     * the arguments appropriate for calling the
     * program through Main
     *
     * @param args     List<String> of command line arguments
     * @param toolname name of the tool to test
     *
     * @return String[] of command line arguments
     */
    public String[] makeCommandLineArgs(final List<String> args, final String toolname) {
        List<String> curatedArgs = injectDefaultVerbosity(args);
        final String[] commandLineArgs = new String[curatedArgs.size() + 1];
        commandLineArgs[0] = toolname;
        int i = 1;
        for (final String arg : curatedArgs) {
            commandLineArgs[i++] = arg;
        }
        return commandLineArgs;
    }

    /**
     * Look for --verbosity argument; if not found, supply a default value that minimizes the
     * amount
     * of logging output.
     */
    private List<String> injectDefaultVerbosity(final List<String> args) {

        // global toggle for BunnyLog output.
        BunnyLog.setEnabled(false);

        for (String arg : args) {
            if (arg.equalsIgnoreCase("--" + StandardArgumentDefinitions.VERBOSITY_NAME) || arg
                    .equalsIgnoreCase("-" + StandardArgumentDefinitions.VERBOSITY_NAME)) {
                return args;
            }
        }
        List<String> argsWithVerbosity = new ArrayList<>(args);
        argsWithVerbosity.add("--" + StandardArgumentDefinitions.VERBOSITY_NAME);
        argsWithVerbosity.add(Log.LogLevel.ERROR.name());
        return argsWithVerbosity;
    }

    /**
     * This is using {@link org.magicdgs.thaplv.Main} instead of {@link
     * org.broadinstitute.hellbender.Main}
     */
    public Object runCommandLine(final List<String> args) {
        return new org.magicdgs.thaplv.Main().instanceMain(makeCommandLineArgs(args));
    }

    public Object runCommandLine(final String[] args) {
        return runCommandLine(Arrays.asList(args));
    }

    public Object runCommandLine(final ArgumentsBuilder args) {
        return runCommandLine(args.getArgsList());
    }
}
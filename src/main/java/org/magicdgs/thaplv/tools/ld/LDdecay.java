/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2015 Daniel Gómez-Sánchez
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

package org.magicdgs.thaplv.tools.ld;

import static org.vetmeduni.haplotypetools.utils.Formats.commaFmt;

import org.magicdgs.thaplv.tools.ld.engine.QueueLD;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.vetmeduni.haplotypetools.components.AbstractTool;
import org.vetmeduni.haplotypetools.filters.*;
import org.vetmeduni.haplotypetools.utils.VariantLogger;

import java.io.File;

/**
 * Class to compute LD in bins in a concurrent way
 * TODO: LinkedList is actually faster than MixedQueue and it could handle a whole chromosome easily
 * (1.5 million variants)
 * TODO: implement minor allele frequency
 * TODO: unit tests
 * TODO: better documentation
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class LDdecay extends AbstractTool {
    // default parameters
    protected static int nbins = 1_000;
    protected static int maxDistance = 10_000;
    protected static int minDistance = 0;
    protected static double chiSqrThreshold = 0.95;
    protected static int noMissing = -1;
    // protected static double maf_defaul = 0.05;

    // This tool filter the invariant sites and biallelic
    // TODO: add count filter
    protected static VariantFilterSet filter = new VariantFilterSet() {{
        addFilter(new NonBiallelicFilter());
        addFilter(new InvariantFilter());
    }};


    @Override
    public int run(String[] args) {
        try {
            // Comand line parsing and extracting options
            CommandLine cmd = programParser(args);
            // open the input
            VCFFileReader input = new VCFFileReader(new File(cmd.getOptionValue("i")), false);
            // get the number of samples from the VCFheader
            int numOfSamples = input.getFileHeader().getNGenotypeSamples();
            // obtain the maximum distance to compute LD
            int max = (cmd.hasOption("maxdist")) ? Integer.parseInt(cmd.getOptionValue("maxdist"))
                    : maxDistance;
            // obtain the minimum distance to compute LD
            int min = (cmd.hasOption("mindist")) ? Integer.parseInt(cmd.getOptionValue("mindist"))
                    : minDistance;
            // obtain the bin distance to compute the statistics
            int binDistance =
                    (cmd.hasOption("bindist")) ? Integer.parseInt(cmd.getOptionValue("bindist"))
                            : nbins;
            // obtain the quantile for the chi-square distribution
            double chiSqrQuantile =
                    (cmd.hasOption("chi")) ? Double.parseDouble(cmd.getOptionValue("chi"))
                            : chiSqrThreshold;
            // check if the quantile is correct
            if (chiSqrQuantile < 0 || chiSqrQuantile >= 1) {
                throw new ParseException("Chi-square quantile should be in the range (0, 1)");
            }
            // obtain the minimum counts
            int minSamples = (cmd.hasOption("mins")) ? Integer.parseInt(cmd.getOptionValue("mins"))
                    : noMissing;
            // log if no-missing data is allowed because the minSamples is set to lower than 0
            if (minSamples < 0 || minSamples > numOfSamples) {
                minSamples = numOfSamples;
            }
            // obtain the maf threshold
            //			Double maf = (cmd.hasOption("maf")) ? Double.parseDouble(cmd.getOptionValue("maf")) : maf_defaul;
            //			if(maf < 0 || maf >= 1) {
            //				throw new RuntimeException("MAF threshold should be in the range 0-1");
            //			}
            boolean rmSingletons = cmd.hasOption("rms");
            int nThreads = (cmd.hasOption("t")) ? Integer.parseInt(cmd.getOptionValue("t")) : 1;
            // Generate the queueLD
            QueueLD queue = new QueueLD(cmd.getOptionValue("o"), binDistance, min, max, minSamples,
                    rmSingletons, chiSqrQuantile, nThreads, logger);


            // print the command line
            printCmdLine(args);
            // print some warnings
            logger.warn(
                    "VCF files should contains either haploid individuals or homozygous; other positions will be removed and called as missing");
            logger.warn("Invariant and multi-allelic sites will be filtered");
            // add the minCount filter
            if (minSamples != 0) {
                filter.addFilter(new MissingSamplesCountFilter(minSamples));
            }
            // log if no-missing data is allowed because the minSamples is set to lower than 0
            if (minSamples < 0 || minSamples > numOfSamples) {
                logger.warn(
                        "Minimum number of samples set to lower than 0 or higher than the actual number of samples in the VCF: computing for no-missing data");
            }
            if (chiSqrQuantile == 0) {
                logger.warn(
                        "Computation for r2 will be performed for all maximum r2 independently of its value (except for 0) because the chi-square quantile is set to 0");
            }


            // Generate the iterator
            CloseableIterator<VariantContext> it;
            if (cmd.hasOption("chr")) {
                String chr = cmd.getOptionValue("chr");
                int chrLength = input.getFileHeader().getSequenceDictionary().getSequence(chr)
                        .getSequenceLength();
                it = input.query(chr, 1, chrLength);
                logger.info("Analysing only chromosome ", chr, " (", commaFmt.format(chrLength),
                        " bp)");
            } else {
                it = input.iterator();
                logger.info("Analysing all chromosomes");
            }
            // start a progress logger
            // TODO: change how many times the logger is triggered
            VariantLogger progress = new VariantLogger(logger, 100_000, "Loaded");
            // loggin some information
            logger.info("Computing LD for variants with at least ", minSamples,
                    " genotypes and a maximum r2 into the ", chiSqrQuantile,
                    " chi-square quantile");

            if (max < 1 && min == 0) {
                logger.info("All variants within the chromosome will be used to compute LD");
            } else {
                logger.info("Only variants with a distance between them in the range [", min, ",",
                        (max < 1) ? "chromosome lenght" : max, "] bp will be considered");
            }
            logger.info("Binning results by ", commaFmt.format(binDistance), " bp");

            if (rmSingletons) {
                logger.info("Singletons will be removed");
                filter.addFilter(new SingletonFilter());
            } else {
                logger.info("Singletons will not be removed");
            }

            // start parsing
            while (it.hasNext()) {
                // read the variant
                VariantContext var = it.next();
                // using the genotype filter to reliable call biallelic SNPs even if the ALT column have more than one
                if (!filter.isFilter(var.getGenotypes())) {
                    queue.add(var);
                }
                if (progress.variant(var)) {
                    logger.debug("Variants in RAM: ", commaFmt.format(queue.size()));
                }
            }
            progress.logNumberOfVariantsProcessed();
            logger.debug("Variants in RAM after loading everything: ",
                    commaFmt.format(queue.size()));
            filter.report(logger);
            // close the queue
            queue.close();
        } catch (IllegalArgumentException | ParseException e) {
            logger.debug(e);
            printUsage(e.getMessage());
        } catch (Exception e) {
            logger.debug(e);
            logger.error(e.getMessage());
            return 2;
        }

        return 0;
    }

    // TODO: add the rest of the options
    @Override
    protected Options programOptions() {
        // the input file
        Option input = Option.builder("i")
                .longOpt("input")
                .desc("Input VCF file (indexed).")
                .hasArg()
                .numberOfArgs(1)
                .argName("INPUT.vcf")
                .required()
                .build();
        // the output prefix
        Option output = Option.builder("o")
                .longOpt("output")
                .desc(String.format("Output prefix for LD results"))
                .hasArg()
                .numberOfArgs(1)
                .argName("PREFIX")
                .required()
                .build();
        // the minimum distance to compute LD
        Option minDist = Option.builder("mindist")
                .longOpt("minimum-distance")
                .desc("Minimum distance (inclusive) between SNPs to compute LD in bp. [Default: "
                        + minDistance + "]")
                .hasArg()
                .numberOfArgs(1)
                .argName("INT")
                .required(false)
                .type(int.class)
                .build();
        // the maximum distance to compute LD
        Option maxDist = Option.builder("maxdist")
                .longOpt("maximum-distance")
                .desc("Maximum distance between SNPs to compute LD in bp. If set to < 1, compute for all the SNPs within the chromosome from minimum distance [Default: "
                        + maxDistance + "]")
                .hasArg()
                .numberOfArgs(1)
                .argName("INT")
                .required(false)
                .type(int.class)
                .build();
        // the bin size
        Option binDist = Option.builder("bindist")
                .longOpt("bin-distance")
                .desc("Distance between pairs to be considered in the same bin in bp. [Default: "
                        + nbins + "]")
                .hasArg()
                .numberOfArgs(1)
                .argName("INT")
                .required(false)
                .type(int.class)
                .build();
        // chi-square threshold
        Option chiThr = Option.builder("chi")
                .longOpt("chi-square")
                .desc("Chi-square quantile to assess the significance of r2max and compute LD. [Default: "
                        + chiSqrThreshold + "]")
                .hasArg()
                .numberOfArgs(1)
                .argName("DOUBLE")
                .required(false)
                .type(double.class)
                .build();
        // the minimum number of samples
        Option minSamples = Option.builder("mins")
                .longOpt("minimum-samples")
                .desc("Minimum number of samples available to compute LD. Setting to a negative number only use no-missing data [Default: "
                        + noMissing + "]")
                .hasArg()
                .numberOfArgs(1)
                .argName("INT")
                .required(false)
                .type(int.class)
                .build();
        // no singletons
        Option singletons = Option.builder("rms")
                .longOpt("remove-singletons")
                .desc("Remove pairs where the minor variant is a singleton.")
                .hasArg(false)
                .required(false)
                .build();
        //		Option maf = Option.builder("maf")
        //								.longOpt("minimum-allele-frequency")
        //								.desc("Minimum allele frequency for a marker to compute LD. [Default:"+maf_defaul+"]")
        //								.hasArg()
        //								.numberOfArgs(1)
        //								.argName("double")
        //								.required(false)
        //								.type(double.class)
        //								.build();
        // flag to only output the binned results
        //		Option onlyBin = Option.builder("bin")
        //							   .longOpt("only-bin")
        //							   .desc("Output only binned results.")
        //							   .hasArg(false)
        //							   .required(false)
        //							   .build();
        Option chromosome = Option.builder("chr")
                .longOpt("chromosome")
                .desc("Perform the LD analysis only in this chromosome.")
                .hasArg()
                .numberOfArgs(1)
                .argName("CHROM")
                .required(false)
                .build();
        Option threads = Option.builder("t")
                .longOpt("threads")
                .desc("Number of threads to use in the computation [Default: " + 1 + "]")
                .hasArg()
                .numberOfArgs(1)
                .argName("THREADS")
                .required(false)
                .build();


        Options options = new Options();

        options.addOption(threads);
        options.addOption(chromosome);
        // options.addOption( onlyBin );
        // options.addOption( maf );
        options.addOption(singletons);
        options.addOption(minSamples);
        options.addOption(chiThr);
        options.addOption(maxDist);
        options.addOption(minDist);
        options.addOption(binDist);
        options.addOption(output);
        options.addOption(input);

        return options;
    }
}

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

package org.magicdgs.thaplv.tools.imputation.engine;

import htsjdk.samtools.util.Log;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.vetmeduni.thaplv.haplotypes.differences.DistancePair;
import org.vetmeduni.thaplv.stats.knn.DistanceKNN;
import org.vetmeduni.thaplv.stats.knn.KNNresult;
import org.vetmeduni.thaplv.utils.htsjdk.AlleleUtils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Collector for variants to be imputed and the imputation process
 *
 * @author Daniel Gómez-Sánchez
 */
public class ImputationCollector {

    // raw variants are stored here
    protected final LinkedList<VariantContext> rawVariantList;

    // imputed variants are stored here
    protected final LinkedList<VariantContext> imputedVariants;

    protected final Log logger = Log.getInstance(ImputationCollector.class);

    protected final int windowSize;

    protected final DistanceKNN KNNcomputer;

    protected final boolean addFormat;

    /**
     * The ID for the imputed format in the output VCF
     */
    @Deprecated
    public static final String IMPUTED_FORMAT_ID = "IMP";

    /**
     * The ID for the number of imputed samples in the output VCF
     */
    public static final String IMPUTED_INFO_ID = "IMP";

    /**
     * The ID for the format field for the allele distances for imputation
     */
    public static final String DISTANCE_FORMAT_ID = "KNNDIST";

    /**
     * The ID for the format field for the allele count used for imputation
     */
    public static final String FREQUENCY_FORMAT_ID = "KNNCOUNT";

    /**
     * The ID for the format field for the KNN input (see {@link #KNN_FORMAT_ORDER}).
     */
    public static final String KNN_FORMAT_ID = "KNNINPUT";

    /**
     * The order for the {@link #KNN_FORMAT_ID}
     */
    public static final String[] KNN_FORMAT_ORDER = new String[] {
            "Number of neighbours detected with distance threshold",
            "Number of samples with available information"};

    /**
     * The order of the {@link #IMPUTED_FORMAT_ID} if included
     *
     * @deprecated use the new format for the new implementation
     */
    @Deprecated
    public static final String[] IMPUTED_FORMAT_ORDER_DEPRECATED =
            new String[] {"K-nearest distance",
                    "Number of neighbours with distance lower or equal than K-nearest",
                    "Number of samples with available information",
                    "Imputed allele frequency in neighbours"};

    /**
     * Initialize an imputation collector with the window
     *
     * @param windowSize the window size to use around the variants to impute
     * @param K          the number of neighbours to compute the distance
     * @param addFormat  add the format line to the imputed haplotypes
     */
    public ImputationCollector(int windowSize, int K, boolean addFormat) {
        rawVariantList = new LinkedList<>();
        imputedVariants = new LinkedList<>();
        this.windowSize = windowSize;
        this.KNNcomputer = new DistanceKNN(K);
        this.addFormat = addFormat;
    }

    public ImputationCollector(int windowSize, int K) {
        this(windowSize, K, false);
    }

    /**
     * Add the variant and return a List of already imputed variants
     *
     * @param variant the variant to add
     *
     * @return list with imputed variants (it could be empty)
     */
    public ArrayList<VariantContext> addVariant(VariantContext variant) {
        // final VariantContext toAdd = new VariantContextBuilder(variant).attribute(IMPUTED_FORMAT_ID, 0).make();
        // logger.debug("addVariant(", variant.getContig(), ":", variant.getStart(), ") -> ", variant.getAttribute(IMPUTED_FORMAT_ID));
        final ArrayList<VariantContext> toReturn = new ArrayList<>();
        if (rawVariantList.isEmpty()) {
            // logger.debug("Initializing empty raw list with ", variant.getContig(), ":", variant.getStart());
            rawVariantList.add(variant);
            return toReturn;
        }
        VariantContext notImputed = null;
        // for each of the no imputed variants
        for (int i = imputedVariants.size(); i < rawVariantList.size(); i++) {
            notImputed = rawVariantList.get(i);
            // flush till not imputed
            flushTillVariant(notImputed, toReturn);
            // check if we will need this for imputation
            if (variantsWithinWindow(variant, notImputed)) {
                // logger.debug("Breaking loop at ", notImputed.getContig(), ":", notImputed.getStart());
                break;
            }
            // impute if don't and continue
            // logger.debug("Imputing ", notImputed.getContig(), ":", notImputed.getStart());
            imputeMaybe(notImputed);
        }
        // logger.debug("Raw variants: ", rawVariantList.size(), "; Imputed variants: ", imputedVariants.size());
        // if we already imputed all the raw variants
        if (imputedVariants.size() == rawVariantList.size()) {
            // logger.debug("Flushing imputed");
            // output all the variants
            flushImputed(toReturn);
        } else {
            flushTillVariant(notImputed, toReturn);
        }
        // logger.debug("Raw variants: ", rawVariantList.size(), "; Imputed variants: ", imputedVariants.size(),
        // " (after minimizing)");
        // now we can add the variant
        // logger.debug("Added variant to raw: ", variant.getContig(), ":", variant.getStart());
        rawVariantList.add(variant);
        // and return the final
        return toReturn;
    }

    /**
     * Iterate over the rest of the variants, imputing till the raw variant list is empty or all the
     * variants are
     * imputed
     *
     * @param collector it will contain the rest of the variants already imputed
     */
    public void finalizeImputation(final ArrayList<VariantContext> collector) {
        logger
                .debug("Finalizing. Raw variants: ", rawVariantList.size(), "; Imputed variants: ",
                        imputedVariants.size());
        // if there is no more variants to impute
        if (rawVariantList.size() == imputedVariants.size()) {
            logger.debug("Every variant should be imputed at this point");
            // flush the imputed variants, close the writer and return
            flushImputed(collector);
            logger.debug("Raw variants: ", rawVariantList.size(), "; Imputed variants: ",
                    imputedVariants.size(),
                    collector.size());
            return;
        }
        logger.debug("There are still some variants without imputation: ",
                rawVariantList.size() - imputedVariants.size());
        // get the variant to impute
        final VariantContext toImpute = rawVariantList.get(imputedVariants.size());
        // check if the inputed variants are still within the window
        while (imputedVariants.size() != 0) {
            if (!variantsWithinWindow(toImpute, imputedVariants.peekFirst())) {
                // if not, just output the imputed and remove the from the raw variants
                collector.add(imputedVariants.pollFirst());
                rawVariantList.pollFirst();
            } else {
                break;
            }
        }
        // at this point, in the queue there are only variants that should be used to impute
        imputeMaybe(toImpute);
        // try again to finalize
        finalizeImputation(collector);
    }

    /**
     * Impute a variant if it should be imputed
     *
     * @param variant the variant that could be imputed
     */
    private void imputeMaybe(final VariantContext variant) {
        HashSet<String> samplesToImpute = new HashSet<>();
        HashSet<String> samplesToUse = new HashSet<>();
        for (Genotype geno : variant.getGenotypes()) {
            if (geno.isNoCall() || !geno.isHom()) {
                samplesToImpute.add(geno.getSampleName());
            } else {
                samplesToUse.add(geno.getSampleName());
            }
        }
        if (samplesToImpute.isEmpty()) {
            imputedVariants.add(variant);
            logger.debug("Complete variant at ", variant.getContig(), ":", variant.getStart());
            // printStats(variant, samplesToImpute.size(), "NA", false);
        } else {
            logger.debug("Imputing variant at ", variant.getContig(), ":", variant.getStart(),
                    " with ",
                    samplesToUse.size(), " samples");
            imputedVariants.add(impute(variant, samplesToImpute, samplesToUse));
        }
    }

    /**
     * Impute a variant using {@link #impute(htsjdk.variant.variantcontext.Genotype,
     * org.vetmeduni.thaplv.stats.knn.KNNresult, htsjdk.variant.variantcontext.VariantContext)}
     *
     * @param variant         the variant to impute
     * @param samplesToImpute the samples that should be imputed
     * @param samplesToUse    the samples to compute the KNN
     *
     * @return the imputed variant; the same variant if no variants are around
     */
    private VariantContext impute(final VariantContext variant, Set<String> samplesToImpute,
            Set<String> samplesToUse) {
        GenotypesContext genotypes = GenotypesContext.create(variant.getNSamples());
        final ArrayList<VariantContext> inWindow = getVariantsInWindow(variant);
        if (inWindow.size() == 0) {
            logger.warn("No variants around ", variant.getContig(), ":", variant.getStart(),
                    ". Position won't be imputed");
            return variant;
        }
        int imputed = 0;
        for (Genotype geno : variant.getGenotypes()) {
            if (samplesToImpute.contains(geno.getSampleName())) {
                final KNNresult knn =
                        KNNcomputer.computeKNN(geno.getSampleName(), samplesToUse, inWindow);
                logger.debug(knn, " at ", variant.getContig(), ":", variant.getStart());
                final Genotype imputedGenotype = impute(geno, knn, variant);
                if (imputedGenotype.isCalled()) {
                    imputed++;
                }
                genotypes.add(imputedGenotype);
            } else {
                genotypes.add(geno);
            }
        }
        // create the imputed variant
        return new VariantContextBuilder(variant).genotypes(genotypes)
                .attribute(IMPUTED_INFO_ID, imputed).make();
    }

    /**
     * New implementation for imputation: compute the average distance for each allele for the
     * neighbours and use that
     * allele as the imputed one
     *
     * @param genotype  the genotype to impute
     * @param knnResult result for KNN
     * @param variant   the variant with the information for all samples
     *
     * @return the imputed genotype (could be not imputed if several alleles have the same score)
     */
    private Genotype impute(final Genotype genotype, final KNNresult knnResult,
            final VariantContext variant) {
        // get the neighbours
        DistancePair[] neighbours = knnResult.getAllNeighboursWithKnearestDistance();
        // get the genotypes for the neighbours
        GenotypesContext neighboursGenotypes = variant
                .getGenotypes(Stream.of(neighbours).map(DistancePair::getSample2)
                        .collect(Collectors.toSet()));
        logger.debug("Imputing with ", neighbours.length, "/", knnResult.numberOfComparisons(),
                " haplotypes: ",
                Arrays.toString(neighboursGenotypes.getSampleNames().toArray()));
        HashMap<Allele, Mean> alleleScores = new HashMap<>();
        for (DistancePair n : neighbours) {
            // get the allele for this neighbour
            final Allele a = neighboursGenotypes.get(n.getSample2()).getAllele(0);
            // get the current mean object
            Mean currentMean = alleleScores.putIfAbsent(a, new Mean());
            if (currentMean == null) {
                currentMean = alleleScores.get(a);
            }
            // increment the mean
            currentMean.increment(n.getDistance());
        }
        // get the minimum value
        double max = alleleScores.values().stream().mapToDouble(Mean::getResult).min()
                .orElse(Double.NaN);
        if (Double.isNaN(max)) {
            logger.debug("No maximum value for an allele: ", alleleScores);
            throw new IllegalArgumentException("Something went wrong when imputing");
        }
        final GenotypeBuilder builder = new GenotypeBuilder(genotype);
        final HashMap<Allele, Integer> counts =
                AlleleUtils.computeAlleleCounts(neighboursGenotypes);
        // add the format if requested
        if (addFormat) {
            // adding the distance for the KNN imputation
            final List<Allele> variantAlleles = variant.getAlleles();
            final float[] dist = new float[variantAlleles.size()];
            final int[] c = new int[variantAlleles.size()];
            int i = 0;
            for (Allele al : variantAlleles) {
                final Mean mean = alleleScores.get(al);
                dist[i] = (mean == null) ? Float.NaN : (float) mean.getResult();
                c[i++] = counts.getOrDefault(al, 0);
            }
            // adding also the KNN information
            builder.attribute(DISTANCE_FORMAT_ID, dist).attribute(FREQUENCY_FORMAT_ID, c)
                    .attribute(KNN_FORMAT_ID,
                            new int[] {neighbours.length, knnResult.numberOfComparisons()});
        }
        Set<Allele> bestAlleles =
                alleleScores.entrySet().stream().filter(e -> e.getValue().getResult() == max)
                        .map(Map.Entry::getKey).collect(Collectors.toSet());
        // if there is no allele closer
        if (bestAlleles.size() != 1) {
            logger.debug("Using frequency for resolve tie at ", variant.getContig(), ":",
                    variant.getStart(), " for ",
                    genotype.getSampleName());
            // compute the most frequent alleles
            Set<Allele> frequentAlleles = AlleleUtils.getMostFrequentAlleles(counts);
            logger.debug("Found ", frequentAlleles.size(), " frequent alleles");
            bestAlleles.retainAll(frequentAlleles);
            if (bestAlleles.size() != 1) {
                logger.warn("Best alleles for imputation does not differ in frequency ",
                        genotype.getSampleName(), " for neighbours at ",
                        variant.getContig(), ":", variant.getStart(), ". Sample won't be imputed");
                return builder.make();
            }
        }
        return builder.alleles(
                AlleleUtils.duplicateByPloidy(bestAlleles.iterator().next(), genotype.getPloidy()))
                .make();
    }

    /**
     * Impute a variant
     *
     * @param variant         the variant to impute
     * @param samplesToImpute the samples that should be imputed
     * @param samplesToUse    the samples to compute the KNN
     *
     * @return the imputed variant; the same variant if no variants are around
     *
     * @deprecated implementation using most frequent allele
     */
    @Deprecated
    private VariantContext imputeDeprecated(final VariantContext variant,
            Set<String> samplesToImpute,
            Set<String> samplesToUse) {
        GenotypesContext genotypes = GenotypesContext.create(variant.getNSamples());
        final boolean computeNeighbours = samplesToUse.size() > KNNcomputer.getK();
        // get the variants that should be used to compute distances
        final ArrayList<VariantContext> inWindow;
        final Allele allSamplesAllele;
        final HashMap<Allele, Integer> alleleCounts;
        if (computeNeighbours) {
            allSamplesAllele = null;
            inWindow = getVariantsInWindow(variant);
            // TODO: not imputed variant warning
            if (inWindow.size() == 0) {
                logger.warn("No variants around ", variant.getContig(), ":", variant.getStart(),
                        ". Position won't be imputed");
                return variant;
            }
            alleleCounts = null;
        } else {
            alleleCounts = AlleleUtils.computeAlleleCounts(variant.getGenotypes(samplesToUse));
            final Set<Allele> alleleSet = AlleleUtils.getMostFrequentAlleles(alleleCounts);
            // TODO: change this to add the info to the VCF
            logger
                    .warn("Only ", samplesToUse.size(), " samples found at ", variant.getContig(),
                            ":", variant.getStart(),
                            "; this samples will be the ", KNNcomputer.getK(),
                            "-nearest neighbours independently of the distance for all samples to impute");
            // The allele won't be imputed because there are two alleles
            if (alleleSet.size() != 1) {
                logger.warn("There is no major allele to impute ", variant.getContig(), ":",
                        variant.getStart(),
                        ". Position won't be imputed");
                return variant;
            }
            allSamplesAllele = alleleSet.iterator().next();
            inWindow = null;
        }
        int imputed = 0;
        for (Genotype geno : variant.getGenotypes()) {
            if (samplesToImpute.contains(geno.getSampleName())) {
                if (computeNeighbours) {
                    KNNresult knn =
                            KNNcomputer.computeKNN(geno.getSampleName(), samplesToUse, inWindow);
                    logger.debug(knn, " at ", variant.getContig(), ":", variant.getStart());
                    final Genotype imputedGenotype = imputeDeprecated(geno, knn, variant);
                    if (imputedGenotype.isCalled()) {
                        imputed++;
                    }
                    genotypes.add(imputedGenotype);
                } else {
                    imputed++;
                    genotypes.add(
                            imputeDeprecated(geno, allSamplesAllele, null, samplesToUse.size(),
                                    samplesToUse.size(),
                                    alleleCounts));
                }
            } else {
                genotypes.add(geno);
            }
        }
        // create the imputed variant
        return new VariantContextBuilder(variant).genotypes(genotypes)
                .attribute(IMPUTED_INFO_ID, imputed).make();
    }

    /**
     * Impute using the provided allele (making sure that it is not in half of the samples)
     *
     * @param genotype            the genotype to impute
     * @param allele              the allele to include; {@link htsjdk.variant.variantcontext.Allele#NO_CALL}
     *                            if
     *                            information should be added but it is not imputed
     * @param distance            the distance for the k-nearest neighbour or <code>null</code> if
     *                            not computed
     * @param nNeighbours         number of k-nearest neighbours
     * @param numberOfComparisons number of haplotypes used
     * @param counts              the counts for all the alleles
     *
     * @return a new Genotype with the inputed allele as the genotype
     *
     * @deprecated implementation using most frequent allele
     */
    @Deprecated
    private Genotype imputeDeprecated(final Genotype genotype, final Allele allele, Double distance,
            int nNeighbours,
            int numberOfComparisons, final HashMap<Allele, Integer> counts) {
        GenotypeBuilder builder = new GenotypeBuilder(genotype)
                .alleles(AlleleUtils.duplicateByPloidy(allele, genotype.getPloidy()));
        if (addFormat) {
            final double frequency = (counts == null || Allele.NO_CALL.equals(allele)) ?
                    Double.NaN :
                    counts.get(allele) / counts.values().stream().mapToDouble(Double::valueOf)
                            .sum();
            addAttribute(builder, distance, nNeighbours, numberOfComparisons, frequency);
        }
        return builder.make();
    }

    /**
     * Impute using the provided allele
     *
     * @param genotype  the genotype to impute
     * @param knnResult result for KNN
     * @param variant   the variant with the information for all samples
     *
     * @return a new Genotype with the inputed allele as the genotype
     *
     * @deprecated implementation using most frequent allele
     */
    @Deprecated
    private Genotype imputeDeprecated(final Genotype genotype, final KNNresult knnResult,
            final VariantContext variant) {
        Set<String> neighbours = knnResult.getAllNeighboursNamesWithKnearest();
        logger.debug("Imputing with ", neighbours.size(), "/", knnResult.numberOfComparisons(),
                " haplotypes: ",
                Arrays.toString(neighbours.toArray()));
        final HashMap<Allele, Integer> counts =
                AlleleUtils.computeAlleleCounts(variant.getGenotypes(neighbours));
        final Set<Allele> alleleSet = AlleleUtils.getMostFrequentAlleles(counts);
        if (alleleSet.size() != 1) {
            logger.warn("There is no major allele to impute ", variant.getContig(), ":",
                    variant.getStart(), " for ",
                    genotype.getSampleName(), ". Sample won't be imputed");
            logger.warn("Allele counts: ", counts);
            return imputeDeprecated(genotype, Allele.NO_CALL, knnResult.getkNearestDistance(),
                    neighbours.size(),
                    knnResult.numberOfComparisons(), counts);
        }
        final Allele mostFrequent = alleleSet.iterator().next();
        return imputeDeprecated(genotype, mostFrequent, knnResult.getkNearestDistance(),
                neighbours.size(),
                knnResult.numberOfComparisons(), counts);
    }

    /**
     * Add the attribute in the order described in {@link #IMPUTED_FORMAT_ORDER_DEPRECATED}
     *
     * @param builder             the builder to add the info to
     * @param distance            the k nearest distance
     * @param nNeighbours         the number of neighbours
     * @param numberOfComparisons the number of samples used to KNN
     * @param frequency           the frequency
     *
     * @return the builder with the attribute added
     *
     * @deprecated now the attributes are different
     */
    @Deprecated
    private GenotypeBuilder addAttribute(GenotypeBuilder builder, Number distance,
            Number nNeighbours,
            Number numberOfComparisons, Number frequency) {
        final float[] info = new float[IMPUTED_FORMAT_ORDER_DEPRECATED.length];
        info[0] = distance.floatValue();
        info[1] = nNeighbours.floatValue();
        info[2] = numberOfComparisons.floatValue();
        info[3] = frequency.floatValue();
        return builder.attribute(IMPUTED_FORMAT_ID, info);
    }

    /**
     * Flush the queues till a variant appear and collect into the collector
     *
     * @param stop      the variant where we should stop
     * @param collector the collector to add the imputed variants to
     */
    private void flushTillVariant(final VariantContext stop,
            final ArrayList<VariantContext> collector) {
        // if the not imputed variant is far from the first
        while (!imputedVariants.isEmpty() && !variantsWithinWindow(imputedVariants.peekFirst(),
                stop)) {
            collector.add(imputedVariants.pollFirst());
            rawVariantList.removeFirst();
        }
    }

    /**
     * Flush all the variants that are imputed, both from the raw and the imputed ones
     *
     * @param collector the collector to add the imputed variants to
     *
     * @throws IllegalArgumentException if it is called with different sizes in the lists
     */
    private void flushImputed(ArrayList<VariantContext> collector) {
        if (imputedVariants.size() != rawVariantList.size()) {
            throw new IllegalStateException(
                    "Callign flushImputed() when not all the variants are imputed");
        }
        while (!imputedVariants.isEmpty()) {
            collector.add(imputedVariants.pollFirst());
            rawVariantList.removeFirst();
        }
    }

    /**
     * Compare the position of two variants to check if they are in the same chromosome and in
     * separated as a maximum of
     * {@link #windowSize}
     *
     * @param variant1 the first variant
     * @param variant2 the second variant
     *
     * @return <code>true</code> if they are only {@link #windowSize} bp appart; <code>false</code>
     * otherwise
     */
    private boolean variantsWithinWindow(VariantContext variant1, VariantContext variant2) {
        // only compare if the chromosomes are equal
        if (variant1.getContig().equals(variant2.getContig())) {
            if (Math.abs(variant1.getStart() - variant2.getStart()) <= windowSize) {
                return true;
            }
        }
        return false;
    }

    /**
     * Get the variants in the window. Requires that the list starts with a variant that is in the
     * window
     *
     * @param variant the variant around the list should be get
     *
     * @return the list with the variants around the sample
     */
    private ArrayList<VariantContext> getVariantsInWindow(final VariantContext variant) {
        final ArrayList<VariantContext> toReturn = new ArrayList<>();
        for (VariantContext v : rawVariantList) {
            if (variant.equals(v)) {
                // this variant is not considered for the distance
            } else if (variantsWithinWindow(v, variant)) {
                toReturn.add(v);
            } else {
                break;
            }
        }
        logger
                .debug("Number of variants around ", variant.getContig(), ":", variant.getStart(),
                        " is ", toReturn.size());
        return toReturn;
    }
}

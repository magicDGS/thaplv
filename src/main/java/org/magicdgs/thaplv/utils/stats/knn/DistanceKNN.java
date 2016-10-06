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

package org.magicdgs.thaplv.utils.stats.knn;

import org.magicdgs.thaplv.haplotypes.pairs.DifferencesDistancePair;

import htsjdk.samtools.util.Log;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Simple implementation of K-nearest neighbours for distance between haplotypes
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class DistanceKNN {

    // TODO: this logger should be changed
    private final static Log logger = Log.getInstance(DistanceKNN.class);

    // the number of neighbours
    private final int K;

    /**
     * Constructor
     *
     * @param K number of neighbours
     */
    public DistanceKNN(final int K) {
        this.K = K;
    }

    /**
     * Get the number of K to use in the class
     *
     * @return K
     */
    public int getK() {
        return K;
    }

    public KNNresult computeKNN(final String sampleName, final Set<String> samplesToUse,
            final Collection<VariantContext> variants) {
        // initialize the pairs
        final List<DifferencesDistancePair> pairs = new ArrayList<>();
        for (final String other : samplesToUse) {
            final DifferencesDistancePair p = new DifferencesDistancePair(sampleName, other);
            for (final VariantContext v : variants) {
                p.add(v);
            }
            pairs.add(p);
        }
        return new KNNresult<>(sampleName, K, 0, variants.size(), pairs);
    }

    /**
     * Get the names for the nearest neighbours using a collection of variants and the distance
     * computed by {@link DifferencesDistancePair#getDistance()}
     *
     * @param sampleName   the sample to which the KNN should be computed
     * @param samplesToUse the samples agains the sample will be compared
     * @param variants     list of variants for compute the distance
     * @param returnAll    if <code>true</code> return all the samples with the K-nearest distance
     *                     cutoff; otherwise,
     *                     return only {@link #K} samples
     *
     * @return a set of at least size {@link #K} with the name of the nearest neighbours
     *
     * @deprecated use {@link #computeKNN(String, Set, Collection)} instead
     */
    @Deprecated
    private Set<String> nearestNeighbours(final String sampleName, final Set<String> samplesToUse,
            final Collection<VariantContext> variants, boolean returnAll) {
        // initialize the pairs
        final ArrayList<DifferencesDistancePair> pairs = new ArrayList<>();
        for (final String other : samplesToUse) {
            final DifferencesDistancePair p = new DifferencesDistancePair(sampleName, other);
            for (final VariantContext v : variants) {
                p.add(v);
            }
            pairs.add(p);
        }
        // sort the pairs by the distance
        final DifferencesDistancePair[] sortedByDistance =
                pairs.stream().sorted(Comparator.comparing(DifferencesDistancePair::getDistance))
                        .toArray(DifferencesDistancePair[]::new);
        // if we are returning every sample for a cutoff
        if (returnAll) {
            // get the maximum distance
            final double maxDistance =
                    Arrays.stream(sortedByDistance)
                            .mapToDouble(DifferencesDistancePair::getDistance)
                            .toArray()[K
                            - 1];
            logger.debug("Maximum distance (0-", variants.size(), "] = ", maxDistance);
            return Stream.of(sortedByDistance).filter(p -> p.getDistance() <= maxDistance)
                    .map(DifferencesDistancePair::getSample2).collect(Collectors.toSet());
        }
        // this is not the best approach for imputation, but the best for other issues
        return Stream.of(sortedByDistance).limit(K).map(DifferencesDistancePair::getSample2)
                .collect(Collectors.toSet());
    }

    /**
     * Get the names for the nearest neighbours using a collection of variants
     *
     * @param sampleName   the sample to which the KNN should be computed
     * @param samplesToUse the samples agains the sample will be compared
     * @param variants     list of variants for compute the distance
     *
     * @return a set of size {@link #K} with the name nearest neighbours
     *
     * @deprecated use {@link #computeKNN(String, Set, Collection)} instead
     */
    @Deprecated
    public Set<String> kNearestNeighbours(final String sampleName, final Set<String> samplesToUse,
            final Collection<VariantContext> variants) {
        return nearestNeighbours(sampleName, samplesToUse, variants, false);
    }

    /**
     * Get the names for the all the nearest neighbours using a collection of variants (based on
     * the
     * K-nearest neighbour
     * distance)
     *
     * @param sampleName   the sample to which the KNN should be computed
     * @param samplesToUse the samples agains the sample will be compared
     * @param variants     list of variants for compute the distance
     *
     * @return a set of size {@link #K} or higher with the name nearest neighbours
     *
     * @deprecated use {@link #computeKNN(String, Set, Collection)} instead
     */
    @Deprecated
    public Set<String> allNeighbours(final String sampleName, final Set<String> samplesToUse,
            final Collection<VariantContext> variants) {
        return nearestNeighbours(sampleName, samplesToUse, variants, true);
    }


}

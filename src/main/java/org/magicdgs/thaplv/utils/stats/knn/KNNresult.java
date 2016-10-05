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

import org.magicdgs.thaplv.haplotypes.pairs.DistancePair;

import htsjdk.samtools.util.Log;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Class for storing KNN results
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class KNNresult<DistanceType extends DistancePair> {

    private final static Log logger = Log.getInstance(KNNresult.class);

    // the sample used to compute the neighbours
    protected final String sampleName;

    // the number of neighbours used to compute the maximum distance
    protected final int K;

    // the maximum possible distance
    protected final double maximumPossibleDistance;

    // the minimum possible distance
    protected final double minimumPossibleDistance;

    // the kNearestDistance (initialize on demand)
    private double kNearestDistance;

    // the maximum possible distance
    private final DistancePair[] sortedByDistance;

    KNNresult(final String sampleName, final int K, final double minimumPossibleDistance,
            final double maximumPossibleDistance, final ArrayList<DistanceType> pairs) {
        this.sampleName = sampleName;
        this.K = K;
        this.kNearestDistance = Double.NaN;
        this.minimumPossibleDistance = minimumPossibleDistance;
        this.maximumPossibleDistance = maximumPossibleDistance;
        // TODO: this is unchecked, but safe?
        this.sortedByDistance =
                pairs.stream().sorted(Comparator.comparing(DistanceType::getDistance))
                        .toArray(DistancePair[]::new);
    }

    public final String getSampleName() {
        return sampleName;
    }

    public final int getK() {
        return K;
    }

    public final double getMinimumPossibleDistance() {
        return minimumPossibleDistance;
    }

    public final double getMaximumPossibleDistance() {
        return maximumPossibleDistance;
    }

    public final int numberOfComparisons() {
        return sortedByDistance.length;
    }

    public final double getkNearestDistance() {
        if (Double.isNaN(kNearestDistance)) {
            computeKnearestDistance();
        }
        return kNearestDistance;
    }

    protected void computeKnearestDistance() {
        final int index = (sortedByDistance.length < K) ? sortedByDistance.length - 1 : K - 1;
        kNearestDistance = Arrays.stream(sortedByDistance).mapToDouble(DistancePair::getDistance)
                .toArray()[index];
    }

    public Set<String> getKnearestNeighbours() {
        return Stream.of(sortedByDistance).limit(K).map(DistancePair::getSample2)
                .collect(Collectors.toSet());
    }

    public Set<String> getAllNeighboursNamesWithKnearest() {
        return getFilteredByDistance(sortedByDistance, getkNearestDistance())
                .map(DistancePair::getSample2)
                .collect(Collectors.toSet());
    }

    public DistancePair[] getAllNeighboursWithKnearestDistance() {
        return getFilteredByDistance(sortedByDistance, getkNearestDistance())
                .toArray(DistancePair[]::new);
    }

    @Deprecated
    public Set<String> getKnearestNeighboursOrLess() {
        // logger.warn("Experimental neirgbours computation: K or less");
        double nearestDistance = getkNearestDistance();
        // get all neighbours as a DistancePair array
        DistancePair[] neighbours = getFilteredByDistance(sortedByDistance, nearestDistance)
                .toArray(DistancePair[]::new);
        //logger.debug(Arrays.toString(Arrays.stream(neighbours).map(p -> String.format("%s: %s", p.getSample2(), p.getDistance())).toArray()));
        // logger.debug("Nearest distance = ", nearestDistance, " generates ", neighbours.length, " neighbours");
        // iterate till the distance is
        int iterK = K;
        while (neighbours.length > K) {
            // get the previous distance
            nearestDistance = neighbours[--iterK].getDistance();
            neighbours = getFilteredByDistance(sortedByDistance, nearestDistance)
                    .toArray(DistancePair[]::new);
            // logger.debug("Nearest distance = ", nearestDistance, " generates ", neighbours.length, " neighbours");
            if (iterK == 0) {
                return Collections.emptySet();
            }
        }
        return Arrays.stream(neighbours).map(DistancePair::getSample2).collect(Collectors.toSet());
    }

    private static Stream<DistancePair> getFilteredByDistance(final DistancePair[] pairs,
            final double distance) {
        return Stream.of(pairs).filter(p -> p.getDistance() <= distance);
    }

    @Override
    public String toString() {
        return String.format("%s %s-NN result: %s [%s,%s]", sampleName, K, getkNearestDistance(),
                minimumPossibleDistance, maximumPossibleDistance);
    }
}

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
package org.magicdgs.thaplv.utils;

import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.math3.stat.descriptive.rank.PSquarePercentile;

import java.util.Map;
import java.util.TreeMap;

/**
 * Modified from htsjdk.tribble.util.MathUtils.RunningStat to add quantile calculations with Apache Commons Math
 * Change also the oldMean and oldStdDev fields to fields in the method
 * TODO: change mean/stdev implementation for apache.commons.math3 methods?
 * TODO: better documentation
 * TODO: synchronize in a better way?
 * TODO: unit tests
 */
public class RunningStats {
	private double newMean, newStdDev;
    private PSquarePercentile median;
	private Variance sampleVar;
    private TreeMap<Double, PSquarePercentile> quantiles;

	/**
	 * Default constructor. Only computes mean, standard deviation and median
	 */
    public RunningStats() {
		median = new PSquarePercentile(50.0D);
		sampleVar = new Variance();
		quantiles = new TreeMap<Double, PSquarePercentile>();
    }
    
    /**
     * Constructor with several quantiles (1,100)
     * @param quantiles
     */
    public RunningStats(double... quantiles) {
    	this();
    	for(double quan: quantiles) {
    		if(quan != 50.0D) {
    			this.quantiles.put(quan, new PSquarePercentile(quan));
    		}
    	}
    }
    
    /**
	 * Put a new value to update the values. If NaN is introduced, it is ignored
	 *
     * @param x the value to add to the running mean and variance
	 * @return true if the value is pushed; false otherwise
     */
    public synchronized boolean push(double x) {
		if(!Double.isNaN(x)) {
			median.increment(x);
			// See Knuth TAOCP vol 2, 3rd edition, page 232
			if (median.getN() == 1) {
				newMean = x;
				newStdDev = 0.0;
			} else {
				double oldMean = newMean;
				newMean = oldMean + (x - oldMean) / median.getN();
				newStdDev = newStdDev + (x - oldMean) * (x - newMean);
			}
			sampleVar.increment(x);
			for (PSquarePercentile quan : quantiles.values()) {
				quan.increment(x);
			}
			return true;
		}
		return false;
    }

	/**
	 * Clear all the data in this instance
	 */
    public final synchronized void clear() {
		median.clear();
		for(PSquarePercentile q: quantiles.values()) {
			q.clear();
		}
		newMean = 0;
		newStdDev = 0;
    }
    
    public final synchronized long numDataValues() {
		return median.getN();
    }

	/**
	 * Get the actual mean
	 *
	 * @return the current mean; NaN if just initialized or clear
	 */
    public final synchronized double mean() {
		return (median.getN() > 0) ? newMean : Double.NaN;
    }

	/**
	 * Get the actual variance
	 *
	 * @return the current variance; NaN if just initialized or clear
	 */
    public final synchronized double variance() {
    	return ((median.getN() > 1) ? newStdDev / (median.getN() - 1) : Double.NaN);
    }

	/**
	 * Get the sample variance for this object
	 *
	 * @return	the current sample variance; NaN if just initialized or clear
	 */
	public final synchronized double sampleVariance() {
		return sampleVar.getResult();
	}

	/**
	 * Get the median
	 *
	 * @return the current variance; NaN if just initialized or clear
	 */
    public final synchronized double median() {
    	return median.getResult();
    }

	/**
	 * Have this running stat an specific quantile?
	 *
	 * @param quantile	the quantile to check
	 * @return
	 */
    public final synchronized boolean hasQuantile(double quantile) {
    	return quantiles.containsKey(quantile);
    }

	/**
	 * Get a concrete quantile for this running stat. Throws an error if the quantile is not present
	 *
	 * @param quantile the quantile to obtain
	 * @return	the value of the quantile
	 */
    public final synchronized double getQuantile(double quantile) {
    	if(quantile == 50.0D) {
    		return median.getResult();
    	}
    	if(hasQuantile(quantile)) {
    		return quantiles.get(quantile).getResult();
    	}
    	throw new IllegalArgumentException("Running stats must be initialized with the requested quantile");
    }

	/**
	 * Get the current std deviation
	 *
	 * @return the current std deviation; NaN if just initialized or clear
	 */
    public final synchronized double standardDeviation() {
		return Math.sqrt(variance());
    }

	public final synchronized double sampleStandardDeviation() {
		return Math.sqrt(sampleVariance());
	}

	/**
	 * Get all quantiles except the median
	 *
	 * @return a sorted map with the quantile and the value for it
	 */
	public final synchronized TreeMap<Double, Double> getAllQuantiles() {
		TreeMap<Double, Double> toReturn = new TreeMap<>();
		for(Map.Entry<Double, PSquarePercentile> q: quantiles.entrySet()) {
			toReturn.put(q.getKey(), q.getValue().getResult());
		}
		return toReturn;
	}

	/**
	 * Override method for testing purposes. Don't use it
	 *
	 * @return informative representation for testing
	 */
	public String toString() {
		StringBuilder builder = new StringBuilder();
		builder.append("Stats for "); builder.append(numDataValues()); builder.append(" datapoints. ");
		builder.append("Mean: "); builder.append(newMean); builder.append("; ");
		builder.append("Variance: "); builder.append(variance()); builder.append("; ");
		builder.append("Median: "); builder.append(median.getResult()); builder.append("; ");
		builder.append("Quantiles computed: "); builder.append(quantiles.keySet());
		// TODO: comment this
//		builder.append("Quantiles: ");
//		for(Map.Entry<Double, PSquarePercentile> q: quantiles.entrySet()) {
//			builder.append("[");
//			builder.append(q.getKey());	builder.append("="); builder.append(q.getValue().getResult());
//			builder.append("] ");
//		}
		// TODO: uncomment this
		return builder.toString();
	}
}

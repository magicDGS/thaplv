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

import java.text.DecimalFormat;
import java.util.Map;
import java.util.concurrent.ConcurrentSkipListMap;

/**
 * Class that bin a double based on the length, computing mean, median and quantile
 * TODO: better documentation
 * TODO: getter and setters for the binStats
 * TODO: unit tests??
 * TODO: format bins with no data???
 *
 * @author Daniel Gómez-Sánchez
 */
public class LengthBinning {

	// value that contains the bins
	protected ConcurrentSkipListMap<Integer, RunningStats> binStat;

	// value that contains the bin size
	protected int binLength;

	// which quantiles should be computed?
	protected double[] quantiles;

	protected boolean empty;

	/**
	 * Default constructor. The quantiles to be computed are the 0.01, 0.05, 0.95, 0.99
	 * TODO: it will be better to do the defaults 0.005, 0.025, 0.975, 0.995 to have the 99% and 95% CI?
	 *
	 * @param binLength
	 */
	public LengthBinning(int binLength) {
		this(binLength, 1,5,95,99);
	}

	/**
	 * Constructor with quantiles
	 *
	 * @param binLength
	 * @param quantiles
	 */
	public LengthBinning(int binLength, double... quantiles) {
		if(binLength < 1) {
			throw new IllegalArgumentException("Binning cannot be performed for bins lower than 1 bp");
		}
		this.binStat = new ConcurrentSkipListMap<Integer, RunningStats>();
		this.binLength = binLength;
		this.quantiles = quantiles;
		empty = true;
	}

	public boolean isEmpty() {
		return empty;
	}

	/**
	 * Add a value associated with an appropiate length
	 *
	 * @param length	the length associated with this value
	 * @param value		the value to add to the computation of the statistic
	 */
	public void add(int length, double value) {
		int cBin = toBin(length);
		RunningStats stats = binStat.get(cBin);
		if(stats == null) {
			stats = new RunningStats(quantiles);
			binStat.put(cBin, stats);
		}
		if(stats.push(value)){
			empty = false;
		}
	}

	public String formatHeaderBins() {
		StringBuilder builder = new StringBuilder("bin\tn\tmean\tsd\tmedian");
		for(double q: quantiles) {
			builder.append("\tq");
			builder.append(q/100);
		}
		return builder.toString();
	}

	public String formatBins(DecimalFormat digits) {
		StringBuilder builder = new StringBuilder();
		binStat.entrySet().stream().filter(bin -> bin.getValue().numDataValues() != 0).forEach(bin -> {
			builder.append("\t");
			builder.append(bin.getKey());
			builder.append("\t");
			builder.append(bin.getValue().numDataValues());
			builder.append("\t");
			builder.append(digits.format(bin.getValue().mean()));
			builder.append("\t");
			builder.append(digits.format(bin.getValue().sampleStandardDeviation()));
			builder.append("\t");
			builder.append(digits.format(bin.getValue().median()));
			builder.append("\t");
			for (double quantile : bin.getValue().getAllQuantiles().values()) {
				builder.append(digits.format(quantile));
				builder.append("\t");
			}
			builder.append("\n");
		});
		return builder.toString();
	}

	public String formatBins(String chromosome, DecimalFormat digits) {
		StringBuilder builder = new StringBuilder();
		binStat.entrySet().stream().filter(bin -> bin.getValue().numDataValues() != 0).forEach(bin -> {
			builder.append(chromosome);
			builder.append("\t");
			builder.append(bin.getKey());
			builder.append("\t");
			builder.append(bin.getValue().numDataValues());
			builder.append("\t");
			builder.append(digits.format(bin.getValue().mean()));
			builder.append("\t");
			builder.append(digits.format(bin.getValue().sampleStandardDeviation()));
			builder.append("\t");
			builder.append(digits.format(bin.getValue().median()));
			builder.append("\t");
			for (double quantile : bin.getValue().getAllQuantiles().values()) {
				builder.append(digits.format(quantile));
				builder.append("\t");
			}
			builder.append("\n");
		});
		return builder.toString();
	}

	/**
	 * Get the bin for a concrete value
	 *
	 * @param length
	 * @return
	 */
	private int toBin(int length) {
		int bin = (int) Math.floor((length-1)/binLength);
		bin*= binLength;
		bin+= binLength;
		return bin;
	}

	/**
	 * Override method for testing purposes
	 *
	 * @return
	 */
	public String toString() {
		StringBuilder builder = new StringBuilder();
		for(Map.Entry<Integer, RunningStats> bin: binStat.entrySet()) {
			builder.append("\tBin #");
			builder.append(bin.getKey());
			builder.append("\t-> ");
			builder.append(bin.getValue());
			builder.append("\n");
		}
		return builder.toString();
	}


	public void clear() {
		for(RunningStats stat: binStat.values()) {
			stat.clear();
		}
		empty = true;
	}
}

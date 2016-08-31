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

package org.magicdgs.thaplv.utils.ld;


import org.magicdgs.thaplv.haplotypes.light.SNPpair;
import org.magicdgs.thaplv.utils.LengthBinning;

import java.util.ArrayList;
import java.util.List;

/**
 *	Class for binning LD results
 *
 * @author Daniel Gómez-Sánchez
 */
public class LDbinning {
	// list with the bins to write (for some statistics we don't need binning
	public static final ArrayList<String> statsToBin = new ArrayList<String>(){{add("r2"); add("r2norm"); add("rw_pos"); add("rw_neg");}};
	//
	public static final List<Double> quantiles = new ArrayList<Double>(){{add(1D); add(5D); add(95D); add(99D);}};
	// array that contains the binning in the order: positive rw, negative rw, r2 and r2'
	private LengthBinning[] binningStats;

	/**
	 * Initialize this LD binning with certain length
	 *
	 * @param binLength	the length to bin
	 */
	public LDbinning(int binLength) {
		// init a new Array with the length og the stats to bin
		binningStats = new LengthBinning[statsToBin.size()];
		for(int i = 0; i < statsToBin.size(); i++) {
			binningStats[i] = new LengthBinning(binLength, quantiles.stream().mapToDouble(d -> d).toArray());
		}
	}

	/**
	 * Add to a distance bin a value for the statistic idx
	 *
	 * @param idx
	 * @param distance
	 * @param value
	 */
	private synchronized void addToBin(int idx, int distance, double value) {
		binningStats[idx].add(distance, value);
	}


	/**
	 * Add pair and associated values
	 *
	 * @param pair
	 * @param values
	 */
	public void add(SNPpair pair, double[] values) {
		String chromosomeName = pair.getReferenceA();
		if(!chromosomeName.equals(pair.getReferenceB())) {
			// TODO: uncomment when implemented for several chromosomes
			// chromosomeName = String.join("-", chromosomeName, pair.getReferenceB());
			throw new IllegalArgumentException("Output for pairs in different chromosomes not implemented");
		}
		add(pair.getDistance(), values);
	}

	/**
	 * Add a new value
	 *
	 * @param distance
	 * @param values
	 */
	private void add(int distance, double[] values) {
		if(values.length != LDfunctions.rStatisticsOrder.length) {
			throw new IllegalArgumentException("Only "+ LDfunctions.rStatisticsOrder.length+" statistics could be written in the output file");
		}
		// iterate over the values
		for(int i = 0; i < LDfunctions.rStatisticsOrder.length; i++) {
			// update the binnings if this value is binned
			int idx = statsToBin.indexOf(LDfunctions.rStatisticsOrder[i]);
			// if it have an index, added
			if(idx != -1) {
				// Only this needs to be synchronized
				addToBin(idx, distance, values[i]);
			} else {
				// check if it is asked to separate positive and negative values
				int posIndex = statsToBin.indexOf(String.format("%s_pos", LDfunctions.rStatisticsOrder[i]));
				int negIndex = statsToBin.indexOf(String.format("%s_neg", LDfunctions.rStatisticsOrder[i]));
				// bin according to the sign of the value
				// TODO: is this binning correct? Where to put 0? Options:
				// TODO: to negative (like current), to positive, to both or ignore it
				if(posIndex != -1 && values[i] > 0) {
					addToBin(posIndex, distance, values[i]);
				} else if(negIndex != -1 && values[i] <= 0) {
					addToBin(negIndex, distance, values[i]);
				}

			}
		}
	}

	public LengthBinning getBin(String stat) {
		int idx = statsToBin.indexOf(stat);
		if(idx < 0) {
			throw new IllegalArgumentException("Stat is not part of this LDbinning");
		}
		return getBin(idx);
	}

	public LengthBinning getBin(int idx) {
		return binningStats[idx];
	}

	/**
	 * Clear the stats, setting to 0
	 */
	public void clear() {
		for(int i = 0; i < statsToBin.size(); i++) {
			binningStats[i].clear();
		}
	}

}

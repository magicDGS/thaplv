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

package org.magicdgs.thaplv.haplotypes.light;

import htsjdk.variant.variantcontext.VariantContext;

import java.util.stream.IntStream;

/**
 * This class contain information for a pair of SNPs to efficiently compute LD
 * A and B are polarized (A > a and B > b); missing data is not stored
 *
 * Both variants could be in different chromosomes, but computations for them are not allowed yet
 *
 * @author Daniel Gómez-Sánchez
 */
public class SNPpair {

	private String referenceA;
	private int positionA;
	private String referenceB;
	private int positionB;

	/**
	 * Counts for AB, Ab and ab, respectively
	 * ab can be computed with totalCounts-sum(counts)
	 */
	private int[] counts;

	/**
	 * Total counts (without missing)
	 */
	private int totalCounts;

	/**
	 * Constructor from one light-weight genotype and one variant from htsjdk
	 *
	 * @param genotype1	the first genotype
	 * @param variant2  the second genotype
	 */
	public SNPpair(LightGenotype genotype1, VariantContext variant2) {
		this(genotype1, new LightGenotype(variant2));
	}

	/**
	 * Constructor from two variant from htsjdk
	 *
	 * @param variant1	the first genotype
	 * @param variant2  the second genotype
	 */
	public SNPpair(VariantContext variant1, VariantContext variant2) {
		this(new LightGenotype(variant1), variant2);
	}

	/**
	 * Constructor from two light-weight genotypes
	 *
	 * @param genotype1	the first genotype
	 * @param genotype2  the second genotype
	 */
	public SNPpair(LightGenotype genotype1, LightGenotype genotype2) {
		if(genotype1.size() != genotype2.size()) {
			throw new IllegalArgumentException("Cannot initialize a SNPpair with genoptypes with differen lengths");
		}
		int AB = 0;
		int Ab = 0;
		int aB = 0;
		totalCounts = 0;
		referenceA = genotype1.getContig();
		positionA = genotype1.getPosition();
		referenceB = genotype2.getContig();
		positionB = genotype2.getPosition();
		for(int i = 0; i < genotype1.size(); i++) {
			LightGenotype.SNP snp1 = genotype1.getGenotypeAt(i);
			LightGenotype.SNP snp2 = genotype2.getGenotypeAt(i);
			if(snp1 == LightGenotype.SNP.N || snp2 == LightGenotype.SNP.N) {
				// to don't store the total counts
				continue;
			} else if(snp1 == snp2) {
				if(snp1 == LightGenotype.SNP.A) {
					AB++;
				}
			} else {
				if(snp1 == LightGenotype.SNP.A) {
					Ab++;
				} else {
					aB++;
				}
			}
			totalCounts++;
		}
		counts = new int[]{AB, Ab, aB};
		polarizeAlleles();
	}

	/**
	 * Polarize alleles in terms of pA > pa and pB > pb
	 */
	private void polarizeAlleles() {
		// counts are in the order AB, Ab, aB
		if(getMajorAC_A() < getMinorAC_A()) {
			// change Ab <=> ab
			counts[1] = getCount_ab();
			// change AB <=> aB
			swap(0,2);
		}
		if(getMajorAC_B() < getMinorAC_B()) {
			// change Ab <=> AB
			swap(0, 1);
			// chage ab <=> aB
			counts[2] = getCount_ab();
		}
	}

	/**
	 * Swap between two counts. Useful for polarize.
	 * Indexes are from 0 to 2 in the order: AB, Ab, aB
	 *
	 * @param index1 the first genotype to swap
	 * @param index2 the second genotype to swap
	 */
	private void swap(int index1, int index2) {
		int temp = counts[index1];
		counts[index1] = counts[index2];
		counts[index2] = temp;
	}

	/**
	 * Check if at least one of the SNPs is a singleton
	 *
	 * @return true if A, B or both are singletons; false otherwise
	 */
	public boolean oneIsSingleton() {
		return getMinorAC_A() == 1 || getMinorAC_B() == 1;
	}

	/**
	 * Get the reference for the first variant (A)
	 *
	 * @return	the reference
	 */
	public String getReferenceA() {
		return referenceA;
	}

	/**
	 * Get the reference for the second variant (B)
	 *
	 * @return the reference
	 */
	public String getReferenceB() {
		return referenceB;
	}

	/**
	 * Get the position for the firs variant (A)
	 *
	 * @return the position
	 */
	public int getPositionA() {
		return positionA;
	}

	/**
	 * Get the position for the firs variant (B)
	 *
	 * @return the position
	 */
	public int getPositionB() {
		return positionB;
	}

	/**
	 * Get the distance between the two SNPs
	 * TODO: throw an error when they are in different chromosomes
	 *
	 * @return	the distance between the 2 SNPs
	 */
	public int getDistance() {
		return Math.abs(positionA-positionB);
	}

	/**
	 * Check if the first variant (A) is invariant
	 *
	 * @return	true if is invariant; false otherwise
	 */
	public boolean isInvariant_A() {
		return getMajorAC_A() == 0 || getMinorAC_A() == 0;
	}

	/**
	 * Check if the second variant (B) is invariant
	 *
	 * @return true if is invariant; false otherwise
	 */
	public boolean isInvariant_B() {
		return getMajorAC_B() == 0 || getMinorAC_B() == 0;
	}

	/**
	 * Get the total counts
	 *
	 * @return	the number of available "haplotypes"
	 */
	public int getTotalCounts() {
		return totalCounts;
	}

	/**
	 * Get the major allele count for the first variant (A)
	 *
	 * @return	the major allele count
	 */
	public int getMajorAC_A() {
		// counts are in the order AB, Ab, aB, ab
		return getCount_AB() + getCount_Ab();
	}

	/**
	 * Get the major allele count for the second variant (B)
	 *
	 * @return	the major allele count
	 */
	public int getMajorAC_B() {
		// counts are in the order AB, Ab, aB, ab
		return getCount_AB() + getCount_aB();
	}

	/**
	 * Get the minor allele count for the first variant (A)
	 *
	 * @return	the minor allele count
	 */
	public int getMinorAC_A() {
		// counts are in the order AB, Ab, aB, ab
		return getCount_aB() + getCount_ab();
	}

	/**
	 * Get the minor allele count for the second variant (B)
	 *
	 * @return	the minor allele count
	 */
	public int getMinorAC_B() {
		// counts are in the order AB, Ab, aB, ab
		return getCount_Ab() + getCount_ab();
	}

	/**
	 * Get the major allele frequency for the first variant (A)
	 *
	 * @return	the major allele frequency
	 */
	public double getMajorAF_A() {
		// counts are in the order AB, Ab, aB, ab
		return getMajorAC_A()/(double) getTotalCounts();
	}

	/**
	 * Get the major allele frequency for the second variant (B)
	 *
	 * @return	the major allele frequency
	 */
	public double getMajorAF_B() {
		// counts are in the order AB, Ab, aB, ab
		return getMajorAC_B()/(double) getTotalCounts();
	}

	/**
	 * Get the minor allele frequency for the first variant (A)
	 *
	 * @return	the minor allele frequency
	 */
	public double getMinorAF_A() {
		// counts are in the order AB, Ab, aB, ab
		return getMinorAC_A()/(double) getTotalCounts();
	}

	/**
	 * Get the minor allele frequency for the second variant (B)
	 *
	 * @return	the minor allele frequency
	 */
	public double getMinorAF_B() {
		// counts are in the order AB, Ab, aB, ab
		return getMinorAC_B()/(double) getTotalCounts();
	}

	/**
	 * Get the count for the "haplotype" AB
	 *
	 * @return the counts
	 */
	public int getCount_AB() {
		return counts[0];
	}

	/**
	 * Get the count for the "haplotype" ab
	 *
	 * @return the counts
	 */
	public int getCount_ab() {
		return totalCounts - IntStream.of(counts).sum();
	}

	/**
	 * Get the count for the "haplotype" Ab
	 *
	 * @return the counts
	 */
	public int getCount_Ab() {
		return counts[1];
	}

	/**
	 * Get the count for the "haplotype" ab
	 *
	 * @return the counts
	 */
	public int getCount_aB() {
		return counts[2];
	}

	/**
	 * Get the frequency for the "haplotype" AB
	 *
	 * @return the frequency
	 */
	public double getFreq_AB() {
		return getCount_AB()/(double) getTotalCounts();
	}

	/**
	 * Get the frequency for the "haplotype" ab
	 *
	 * @return the frequency
	 */
	public double getFreq_ab() {
		return getCount_ab()/(double) getTotalCounts();
	}

	/**
	 * Get the frequency for the "haplotype" Ab
	 *
	 * @return the frequency
	 */
	public double getFreq_Ab() {
		return getCount_Ab()/(double) getTotalCounts();
	}

	/**
	 * Get the frequency for the "haplotype" aB
	 *
	 * @return the frequency
	 */
	public double getFreq_aB() {
		return getCount_aB()/(double) getTotalCounts();
	}

	/**
	 * Representation of the SNP pair with ref:position
	 *
	 * @return	representation of the SNPpair
	 */
	public String toString() {
		return String.format("%s:%s-%s:%s", referenceA, positionA, referenceB, positionB);
	}
}

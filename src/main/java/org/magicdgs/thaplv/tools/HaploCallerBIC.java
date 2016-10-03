package org.magicdgs.thaplv.tools;

import org.magicdgs.thaplv.cmd.ThaplvArgumentDefinitions;
import org.magicdgs.thaplv.cmd.programgroups.AlphaProgramGroup;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import org.apache.commons.math3.stat.interval.BinomialConfidenceInterval;
import org.apache.commons.math3.stat.interval.ClopperPearsonInterval;
import org.apache.commons.math3.stat.interval.ConfidenceInterval;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.LocusWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.IntStream;

/**
 * Class for calling haplotypes in a similar way as described in
 * <a href=http://onlinelibrary.wiley.com/wol1/doi/10.1111/mec.12594/full>Kapun et al. 2014</a>.
 *
 * This does not apply any threshold by now.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
// TODO: add filters using the coverage distribution of the chromosomes
@CommandLineProgramProperties(summary = "TODO", oneLineSummary = "TODO", programGroup = AlphaProgramGroup.class)
public class HaploCallerBIC extends LocusWalker {

    // Default value is set for maintain compatibility with diploid software
    @Argument(fullName = ThaplvArgumentDefinitions.OUTPUT_PLOIDY_LONG, shortName = ThaplvArgumentDefinitions.OUTPUT_PLOIDY_SHORT, doc = "Ploidy for the output haplotype. Default value is used for maintain compatibility with other software.", common = false, optional = true)
    public int ploidy = 2;

    private SAMFileHeader header;
    private Set<String> samples;

    // TODO: provide an option to compute the confidence interval with other methods?
    // TODO: this is suppose to be the exact one.
    private static final BinomialConfidenceInterval BIC = new ClopperPearsonInterval();

    // TODO: provide this as a parameter?? will be better to adjust
    private static final double CONFIDENCE_LEVEL = 0.9;

    @Override
    public void onTraversalStart() {
        header = getHeaderForReads();
        samples = ReadUtils.getSamplesFromHeader(header);
    }

    @Override
    public void apply(final AlignmentContext alignmentContext,
            final ReferenceContext referenceContext, final FeatureContext featureContext) {
        // TODO: filter pileup and remove the overlaps
        final ReadPileup pileup = alignmentContext.getBasePileup();
        // get the set of alleles in the pileup
        final Allele refAllele = Allele.create(referenceContext.getBase(), true);
        final Set<Allele> alleles = new TreeSet<>();
        alleles.add(refAllele);
        pileup.forEach(p -> {
            final Allele al = Allele.create(p.getBase());
            if (!al.equals(refAllele, true)) {
                alleles.add(al);
            }
        });
        // get the genotypes
        final int refIndex = BaseUtils.simpleBaseToBaseIndex(referenceContext.getBase());
        final List<Genotype> genotypes = new ArrayList<>(samples.size());
        final Map<String, ReadPileup> stratified = pileup.splitBySample(header, null);
        for (final String s : samples) {
            final ReadPileup samPileup = stratified.getOrDefault(s, null);
            if (samPileup == null) {
                genotypes.add(GenotypeBuilder.createMissing(s, ploidy));
            } else {
                final GenotypeBuilder builder = new GenotypeBuilder(s).phased(true);
                // TODO: allele counts should be sub-sample to the ref allele and the second major count
                final int[] alleleCounts = samPileup.getBaseCounts();
                // computing the total depth
                final int depth = IntStream.of(alleleCounts).sum();
                builder.DP(depth);
                // creates the confidence interval
                final ConfidenceInterval interval =
                        BIC.createInterval(depth, alleleCounts[refIndex], CONFIDENCE_LEVEL);
                if (IntStream.of(alleleCounts).min().orElse(0) / (double) depth > interval
                        .getLowerBound()) {
                    // TODO: extract the filter to a constant and create the VCFHeaderLine
                    builder.filter("BICFilter");
                    // TODO: this does not fit the data, so add a filter to the genotype
                }
                // TODO: add the genotype to the builder and the AD tag
                // TODO: add other filters regarding coverage
            }
        }
    }

}

package org.magicdgs.thaplv.utils.stats;

import org.magicdgs.thaplv.utils.test.BaseTest;

import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Map;

/**
 * This test is using a liberal tolerance for the statistics computed with
 * {@link org.apache.commons.math3.stat.descriptive.rank.PSquarePercentile}.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class RunningStatsUnitTest extends BaseTest {

    private static final double QUANTILE_TOLERANCE = 1e-5;

    @Test
    public void test1MStats() throws Exception {
        // TODO: there are some code repetition in this test that could be solved with data providers
        // creating an instance
        final RunningStats stats = new RunningStats(1, 90);
        // assert that it is empty
        Assert.assertEquals(stats.numDataValues(), 0);
        // checking if the quatiles are properly set
        Assert.assertTrue(stats.hasQuantile(1));
        Assert.assertTrue(stats.hasQuantile(90));
        Assert.assertFalse(stats.hasQuantile(5));
        // empty stats returns NaN and empty quantiles
        Assert.assertEquals(stats.mean(), Double.NaN, DEFAULT_TOLERANCE);
        Assert.assertEquals(stats.median(), Double.NaN, QUANTILE_TOLERANCE);
        Assert.assertEquals(stats.variance(), Double.NaN, DEFAULT_TOLERANCE);
        Assert.assertEquals(stats.standardDeviation(), Double.NaN, DEFAULT_TOLERANCE);
        Assert.assertEquals(stats.sampleVariance(), Double.NaN, DEFAULT_TOLERANCE);
        Assert.assertEquals(stats.sampleStandardDeviation(), Double.NaN, DEFAULT_TOLERANCE);
        Assert.assertEquals(stats.getQuantile(1), Double.NaN, QUANTILE_TOLERANCE);
        Assert.assertEquals(stats.getQuantile(90), Double.NaN, QUANTILE_TOLERANCE);
        Assert.assertThrows(IllegalArgumentException.class, () -> stats.getQuantile(5));
        int num = 0;
        for (double i = 0d; i < 1d; i += 0.000001d) {
            stats.push(i);
            // assert than we are incrementing the number of data
            Assert.assertEquals(stats.numDataValues(), ++num);
        }
        Assert.assertEquals(stats.mean(), 0.50, DEFAULT_TOLERANCE);
        Assert.assertEquals(stats.median(), 0.50, QUANTILE_TOLERANCE);
        Assert.assertEquals(stats.variance(), 0.08333358, DEFAULT_TOLERANCE);
        Assert.assertEquals(stats.standardDeviation(), 0.2886756, DEFAULT_TOLERANCE);
        Assert.assertEquals(stats.sampleVariance(), 0.08333358, DEFAULT_TOLERANCE);
        Assert.assertEquals(stats.sampleStandardDeviation(), 0.2886756, DEFAULT_TOLERANCE);
        Assert.assertEquals(stats.getQuantile(1), 0.01, QUANTILE_TOLERANCE);
        Assert.assertEquals(stats.getQuantile(90), 0.90, QUANTILE_TOLERANCE);
        Assert.assertThrows(IllegalArgumentException.class, () -> stats.getQuantile(5));
        final Map<Double, Double> quantiles = stats.getAllQuantiles();
        Assert.assertEquals(stats.getQuantile(1), quantiles.get(1d), QUANTILE_TOLERANCE);
        Assert.assertEquals(stats.getQuantile(90), quantiles.get(90d), QUANTILE_TOLERANCE);
        // testing if clear works
        stats.clear();
        Assert.assertNotEquals(stats.numDataValues(), 0);
    }
}
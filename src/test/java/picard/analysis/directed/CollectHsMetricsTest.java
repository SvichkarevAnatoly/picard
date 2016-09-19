package picard.analysis.directed;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.Histogram;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;

public class CollectHsMetricsTest extends CommandLineProgramTest {
    private final static File TEST_DIR = new File("testdata/picard/analysis/directed/CollectHsMetrics");

    @Override
    public String getCommandLineProgramName() {
        return CollectHsMetrics.class.getSimpleName();
    }

    @DataProvider(name = "collectHsMetricsDataProvider")
    public Object[][] targetedIntervalDataProvider() {
        final String referenceFile = TEST_DIR + "/chrM.fasta";
        final String intervals = TEST_DIR + "/chrM.interval_list";
        final String twoSmallIntervals = TEST_DIR + "/two-small.interval_list";

        return new Object[][] {
                // test that all bases (read 2) with base quality 1 are filtered out
                // TODO: base quality minimum should not be 1 - pick something higher
                {TEST_DIR + "/lowbaseq.sam",    intervals, 1, 10, true,  2, 200, 0.5, 0.0, 0.50, 0.0,   1000},
                // test that read 2 (with mapping quality 1) is filtered out with minimum mapping quality 2
                {TEST_DIR + "/lowmapq.sam",     intervals, 2, 0, true,  2, 202, 0,   0.0, 0.505, 0.0,   1000},
                // test that we clip overlapping bases
                {TEST_DIR + "/overlapping.sam", intervals, 0, 0, true,  2, 202, 0,   0.5, 0.505, 0.505, 1000},
                // test that we do not clip overlapping bases
                {TEST_DIR + "/overlapping.sam", intervals, 0, 0, false, 2, 202, 0,   0.0, 0.505, 0.505, 1000},
                // 0 bin test TODO: better title
                {TEST_DIR + "/single-short-read.sam", twoSmallIntervals, 20, 20, true, 1, 10, 0.0, 0.0, 0.5, 0.0, 1000 }

        };
    }

    @Test(dataProvider = "collectHsMetricsDataProvider")
    public void runCollectTargetedMetricsTest(final String input,
                                              final String targetIntervals,
                                              final int minimumMappingQuality,
                                              final int minimumBaseQuality,
                                              final boolean clipOverlappingReads,
                                              final int totalReads,
                                              final int pfUqBasesAligned,
                                              final double pctExcBaseq,
                                              final double pctExcOverlap,
                                              final double pctTargetBases1x,
                                              final double pctTargetBases2x,
                                              final int sampleSize) throws IOException {

        final File outfile = File.createTempFile("CollectHsMetrics", ".hs_metrics", TEST_DIR);
        outfile.deleteOnExit();

        final String[] args = new String[] {
                "TARGET_INTERVALS=" + targetIntervals,
                "BAIT_INTERVALS=" + targetIntervals,
                "INPUT=" + input,
                "OUTPUT=" + outfile,
                "MINIMUM_MAPPING_QUALITY=" + minimumMappingQuality,
                "MINIMUM_BASE_QUALITY=" + minimumBaseQuality,
                "CLIP_OVERLAPPING_READS=" + clipOverlappingReads,
                "SAMPLE_SIZE=" + sampleSize
        };

        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<HsMetrics, Comparable<?>> output = new MetricsFile<HsMetrics, Comparable<?>>();
        output.read(new FileReader(outfile));

        for (final HsMetrics metrics : output.getMetrics()) {
            // overlap
            Assert.assertEquals(metrics.TOTAL_READS, totalReads);
            Assert.assertEquals(metrics.PF_UQ_BASES_ALIGNED, pfUqBasesAligned);
            Assert.assertEquals(metrics.PCT_EXC_BASEQ, pctExcBaseq);
            Assert.assertEquals(metrics.PCT_EXC_OVERLAP, pctExcOverlap);
            Assert.assertEquals(metrics.PCT_TARGET_BASES_1X, pctTargetBases1x);
            Assert.assertEquals(metrics.PCT_TARGET_BASES_2X, pctTargetBases2x);
        }
    }

    @Test
    public void testCoverageHistogram() throws IOException {
        final String input = TEST_DIR + "/single-short-read.sam";
        final String targetIntervals = TEST_DIR + "/two-small.interval_list";
        final int minimumMappingQuality = 20;
        final int minimumBaseQuality = 20;
        final boolean clipOverlappingReads = true;
        final int sampleSize = 10;

        final File outfile = File.createTempFile("testCoverageHistogram", ".hs_metrics", TEST_DIR);
        outfile.deleteOnExit();

        final String[] args = new String[] {
                "TARGET_INTERVALS=" + targetIntervals,
                "BAIT_INTERVALS=" + targetIntervals,
                "INPUT=" + input,
                "OUTPUT=" + outfile,
                "MINIMUM_MAPPING_QUALITY=" + minimumMappingQuality,
                "MINIMUM_BASE_QUALITY=" + minimumBaseQuality,
                "CLIP_OVERLAPPING_READS=" + clipOverlappingReads,
                "SAMPLE_SIZE=" + sampleSize
        };

        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<HsMetrics, Integer> output = new MetricsFile<>();
        output.read(new FileReader(outfile));
        final Histogram<Integer> coverageHistogram = output.getAllHistograms().get(0);
        Assert.assertEquals(coverageHistogram.get(0).getValue(), 10.0);
        Assert.assertEquals(coverageHistogram.get(1).getValue(), 10.0);
    }
}

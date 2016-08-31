package picard.fingerprint;

import htsjdk.samtools.metrics.MetricsFile;
import org.testng.Assert;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.vcf.SamTestUtils;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;

/**
 * Tests for CrosscheckReadgroupFingerprints
 */
public class CrosscheckReadGroupFingerprintsTest {

    private final static File TEST_DIR = new File("testdata/picard/fingerprint/");
    private final static File HAPLOTYPE_MAP = new File(TEST_DIR, "Homo_sapiens_assembly19.haplotype_database.subset.txt");

    static private final File NA12891_r1_sam = new File(TEST_DIR, "NA12891.over.fingerprints.r1.sam");
    static private final File NA12891_r2_sam = new File(TEST_DIR, "NA12891.over.fingerprints.r2.sam");

    //this is a copy of a previous one, but with a different sample name
    static private final File NA12891_named_NA12892_r1_sam = new File(TEST_DIR, "NA12891_named_NA12892.over.fingerprints.r1.sam");

    static private final File NA12892_r1_sam = new File(TEST_DIR, "NA12892.over.fingerprints.r1.sam");
    static private final File NA12892_r2_sam = new File(TEST_DIR, "NA12892.over.fingerprints.r2.sam");

    static private File NA12891_r1, NA12891_r2, NA12891_named_NA12892_r1, NA12892_r1, NA12892_r2;

    static private final int NA12891_r1_RGs = 27;
    static private final int NA12891_r2_RGs = 26;
    static private final int NA12892_r1_RGs = 25;
    static private final int NA12892_r2_RGs = 26;

    @BeforeTest
    public void setup() throws IOException {
        NA12891_r1 = SamTestUtils.createIndexedBam(NA12891_r1_sam, NA12891_r1_sam);
        NA12891_r2 = SamTestUtils.createIndexedBam(NA12891_r2_sam, NA12891_r2_sam);
        NA12891_named_NA12892_r1 = SamTestUtils.createIndexedBam(NA12891_named_NA12892_r1_sam, NA12891_named_NA12892_r1_sam);
        NA12892_r1 = SamTestUtils.createIndexedBam(NA12892_r1_sam, NA12892_r1_sam);
        NA12892_r2 = SamTestUtils.createIndexedBam(NA12892_r2_sam, NA12892_r2_sam);
    }

    @DataProvider(name = "bamFilesRGs")
    public Object[][] bamFilesRGs() {
        return new Object[][] {
                new Object[]{NA12891_r1, NA12891_r2, false, 0, (NA12891_r1_RGs + NA12891_r2_RGs) * (NA12891_r1_RGs + NA12891_r2_RGs - 1) / 2},
                new Object[]{NA12891_r1, NA12892_r1, false, 0, (NA12891_r1_RGs + NA12892_r1_RGs) * (NA12891_r1_RGs + NA12892_r1_RGs - 1) / 2},
                new Object[]{NA12891_r1, NA12892_r2, false, 0, (NA12891_r1_RGs + NA12892_r2_RGs) * (NA12891_r1_RGs + NA12892_r2_RGs - 1) / 2},
                new Object[]{NA12892_r1, NA12892_r2, false, 0, (NA12892_r1_RGs + NA12892_r2_RGs) * (NA12892_r1_RGs + NA12892_r2_RGs - 1) / 2},
                new Object[]{NA12892_r2, NA12891_r2, false, 0, (NA12892_r2_RGs + NA12891_r2_RGs) * (NA12892_r2_RGs + NA12891_r2_RGs - 1) / 2},
                new Object[]{NA12892_r2, NA12891_r1, false, 0, (NA12892_r2_RGs + NA12891_r1_RGs) * (NA12892_r2_RGs + NA12891_r1_RGs - 1) / 2},
                new Object[]{NA12891_r1, NA12891_r2, true, 0, (NA12891_r1_RGs + NA12891_r2_RGs) * (NA12891_r1_RGs + NA12891_r2_RGs - 1) / 2},
                new Object[]{NA12891_r1, NA12892_r1, true, 1, (NA12891_r1_RGs + NA12892_r1_RGs) * (NA12891_r1_RGs + NA12892_r1_RGs - 1) / 2},
                new Object[]{NA12891_r1, NA12892_r2, true, 1, (NA12891_r1_RGs + NA12892_r2_RGs) * (NA12891_r1_RGs + NA12892_r2_RGs - 1) / 2},
                new Object[]{NA12892_r1, NA12892_r2, true, 0, (NA12892_r1_RGs + NA12892_r2_RGs) * (NA12892_r1_RGs + NA12892_r2_RGs - 1) / 2},
                new Object[]{NA12892_r2, NA12891_r2, true, 1, (NA12892_r2_RGs + NA12891_r2_RGs) * (NA12892_r2_RGs + NA12891_r2_RGs - 1) / 2},
                new Object[]{NA12892_r2, NA12891_r1, true, 1, (NA12892_r2_RGs + NA12891_r1_RGs) * (NA12892_r2_RGs + NA12891_r1_RGs - 1) / 2}
        };
    }

    @Test(dataProvider = "bamFilesRGs")
    public void testCrossCheckRGs(final File file1, final File file2, final boolean expectAllMatch, final int expectedRetVal, final int expectedNMetrics) throws IOException {

        File metrics = File.createTempFile("Fingerprinting", "NA1291.RG.crosscheck_metrics");
        metrics.deleteOnExit();

        final CrosscheckReadGroupFingerprints crossChecker = new CrosscheckReadGroupFingerprints();
        final String[] args = new String[]{
                "INPUT=" + file1.getAbsolutePath(),
                "INPUT=" + file2.getAbsolutePath(),
                "OUTPUT=" + metrics.getAbsolutePath(),
                "HAPLOTYPE_MAP=" + HAPLOTYPE_MAP,
                "LOD_THRESHOLD=" + -2.0,
                "EXPECT_ALL_GROUPS_TO_MATCH=" + expectAllMatch
        };
        Assert.assertEquals(crossChecker.instanceMain(args), expectedRetVal);
        final MetricsFile<CrosscheckMetric, Comparable<?>> metricsOutput = new MetricsFile<>();
        metricsOutput.read(new FileReader(metrics));

        Assert.assertFalse(metricsOutput.getMetrics().stream()
                .anyMatch(m -> m.DATA_TYPE != CrosscheckMetric.DataType.READGROUP));

        if (expectAllMatch) {
            Assert.assertTrue(metricsOutput.getMetrics().stream()
                    .allMatch(m -> m.RESULT == CrosscheckMetric.FingerprintResult.INCONCLUSIVE ||
                            m.RESULT.isMatch() == m.LEFT_SAMPLE.equals(m.RIGHT_SAMPLE)));
        } else if (expectedRetVal == 0) {
            Assert.assertTrue(metricsOutput.getMetrics().stream()
                    .allMatch(m -> m.RESULT == CrosscheckMetric.FingerprintResult.INCONCLUSIVE ||
                            m.RESULT.isExpected()));
        } else {
            Assert.assertTrue(metricsOutput.getMetrics().stream()
                    .anyMatch(m -> !m.RESULT.isExpected()));
        }

        Assert.assertEquals(metricsOutput.getMetrics().size(), expectedNMetrics);
    }

    @DataProvider(name = "bamFilesLBs")
    public Object[][] bamFilesLBs() {

        return new Object[][]{
                new Object[]{NA12891_r1, NA12891_r2, 0},
                new Object[]{NA12891_r1, NA12892_r1, 0},
                new Object[]{NA12892_r2, NA12891_r2, 0},
                new Object[]{NA12892_r2, NA12891_r1, 0},
                new Object[]{NA12891_r1, NA12891_named_NA12892_r1_sam, 1}, //error since expected match but found a mismatch
                new Object[]{NA12892_r1, NA12891_named_NA12892_r1_sam, 1}, //error since expected mismatch but found a match
        };
    }

    @Test(dataProvider = "bamFilesLBs")
    public void testCrossCheckLBs(final File file1, final File file2, final int expectedRetVal) throws IOException {
        File metrics = File.createTempFile("Fingerprinting", "NA1291.LB.crosscheck_metrics");
        metrics.deleteOnExit();

        final CrosscheckReadGroupFingerprints crossChecker = new CrosscheckReadGroupFingerprints();
        final String[] args = new String[]{
                "INPUT=" + file1.getAbsolutePath(),
                "INPUT=" + file2.getAbsolutePath(),
                "OUTPUT=" + metrics.getAbsolutePath(),
                "HAPLOTYPE_MAP=" + HAPLOTYPE_MAP,
                "LOD_THRESHOLD=" + -1.0,
                "CROSSCHECK_LIBRARIES=true"
        };
        Assert.assertEquals(crossChecker.instanceMain(args), expectedRetVal);
        final MetricsFile<CrosscheckMetric, Comparable<?>> metricsOutput = new MetricsFile<>();
        metricsOutput.read(new FileReader(metrics));

        Assert.assertFalse(metricsOutput.getMetrics().stream()
                .anyMatch(m -> m.DATA_TYPE != CrosscheckMetric.DataType.LIBRARY));

        if (expectedRetVal == 0) {
            Assert.assertTrue(metricsOutput.getMetrics().stream()
                    .allMatch(m -> m.RESULT == CrosscheckMetric.FingerprintResult.INCONCLUSIVE ||
                            m.RESULT.isExpected()));
        } else {
            Assert.assertTrue(metricsOutput.getMetrics().stream()
                    .anyMatch(m -> !m.RESULT.isExpected()));
        }
        // make sure that some things are conclusive:
        Assert.assertTrue(metricsOutput.getMetrics().stream()
                .anyMatch(m -> m.RESULT != CrosscheckMetric.FingerprintResult.INCONCLUSIVE));
        Assert.assertEquals(metricsOutput.getMetrics().size(), 1);
    }

    @DataProvider(name = "bamFilesSMs")
    public Object[][] bamFilesSMs() {

        return new Object[][]{
                new Object[]{NA12891_r1, NA12891_r2, 0, 1},
                new Object[]{NA12891_r1, NA12892_r1, 0, 2},
                new Object[]{NA12892_r2, NA12891_r2, 0, 2},
                new Object[]{NA12892_r2, NA12891_named_NA12892_r1, 0, 1}, // no error since only one sample in aggregate
                new Object[]{NA12891_r2, NA12891_named_NA12892_r1, 1, 2}, //expected match
        };
    }

    @Test(dataProvider = "bamFilesSMs")
    public void testCrossCheckSMs(final File file1, final File file2, final int expectedRetVal, final int numberOfSamples) throws IOException {
        File metrics = File.createTempFile("Fingerprinting", "NA1291.SM.crosscheck_metrics");
        metrics.deleteOnExit();

        final CrosscheckReadGroupFingerprints crossChecker = new CrosscheckReadGroupFingerprints();
        final String[] args = new String[]{
                "INPUT=" + file1.getAbsolutePath(),
                "INPUT=" + file2.getAbsolutePath(),
                "OUTPUT=" + metrics.getAbsolutePath(),
                "HAPLOTYPE_MAP=" + HAPLOTYPE_MAP,
                "LOD_THRESHOLD=" + -1.0,
                "CROSSCHECK_SAMPLES=true"
        };
        Assert.assertEquals(crossChecker.instanceMain(args), expectedRetVal);
        final MetricsFile<CrosscheckMetric, Comparable<?>> metricsOutput = new MetricsFile<>();
        metricsOutput.read(new FileReader(metrics));

        Assert.assertFalse(metricsOutput.getMetrics().stream()
                .anyMatch(m -> m.DATA_TYPE != CrosscheckMetric.DataType.SAMPLE));

        if (expectedRetVal == 0) {
            Assert.assertTrue(metricsOutput.getMetrics().stream()
                    .allMatch(m -> m.RESULT == CrosscheckMetric.FingerprintResult.INCONCLUSIVE ||
                            m.RESULT.isExpected()));
        } else {
            Assert.assertTrue(metricsOutput.getMetrics().stream()
                    .anyMatch(m -> !m.RESULT.isExpected()));
        }

        // make sure that some things are conclusive:
        if (numberOfSamples > 1) {
            Assert.assertTrue(metricsOutput.getMetrics().stream()
                    .anyMatch(m -> m.RESULT != CrosscheckMetric.FingerprintResult.INCONCLUSIVE));
        }
        Assert.assertEquals(metricsOutput.getMetrics().size(), numberOfSamples * (numberOfSamples - 1) / 2);
    }
}
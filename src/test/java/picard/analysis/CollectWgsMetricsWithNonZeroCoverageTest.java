package picard.analysis;

import htsjdk.samtools.*;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;


import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import static org.testng.Assert.*;

public class CollectWgsMetricsWithNonZeroCoverageTest extends CommandLineProgramTest {
    private final static File REF_DICT_DIR = new File("testdata/picard/sam/CollectGcBiasMetrics/");
    private final static File TEST_DIR = new File("testdata/picard/sam/");
    private final File referenceDict = new File(REF_DICT_DIR, "MSmallHeader.dict");
    private File tempSamFile;
    private File outfile;

    private final static int READ_PAIR_DISTANCE = 99;
    private final static String SAMPLE = "TestSample1";
    private final static String READ_GROUP_ID = "TestReadGroup1";
    private final static String PLATFORM = "ILLUMINA";
    private final static String LIBRARY = "TestLibrary1";

    public String getCommandLineProgramName() {
        return CollectWgsMetricsWithNonZeroCoverage.class.getSimpleName();
    }

    @Test
    public void testExcludedBases() throws IOException {
        final File reference = new File("testdata/picard/quality/chrM.reference.fasta");
        final File testSamFile = File.createTempFile("CollectWgsMetrics", ".bam", TEST_DIR);
        testSamFile.deleteOnExit();

        /**
         *  Our test SAM looks as follows:
         *
         *   ----------   <- reads with great base qualities (60) ->  ----------
         *   ----------                                               ----------
         *   ----------                                               ----------
         *   **********   <- reads with poor base qualities (10) ->   **********
         *   **********                                               **********
         *   **********                                               **********
         *
         *  We exclude half of the bases because they are low quality.
         *  We do not exceed the coverage cap (3), thus none of the bases is excluded as such.
         *
         */

        final SAMRecordSetBuilder setBuilder = createTestSAMBuilder(reference);
        setBuilder.setReadLength(10);
        for (int i = 0; i < 3; i++){
            setBuilder.addPair("GreatBQRead:" + i, 0, 1, 30, false, false, "10M", "10M", false, true, 60);
        }

        for (int i = 0; i < 3; i++){
            setBuilder.addPair("PoorBQRead:" + i, 0, 1, 30, false, false, "10M", "10M", false, true, 10);
        }

        final SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(setBuilder.getHeader(), false, testSamFile);

        for (final SAMRecord record : setBuilder) {
            writer.addAlignment(record);
        }

        writer.close();

        final File outfile = File.createTempFile("testExcludedBases-metrics", ".txt");
        outfile.deleteOnExit();

        final File chartOutFile = File.createTempFile("testExcludedBases",".pdf");
        chartOutFile.deleteOnExit();


        final String[] args = new String[] {
                "INPUT="  + testSamFile.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "REFERENCE_SEQUENCE=" + reference.getAbsolutePath(),
                "INCLUDE_BQ_HISTOGRAM=true",
                "COVERAGE_CAP=3",
                "CHART_OUTPUT=" + chartOutFile.getAbsolutePath()
        };

        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<CollectWgsMetrics.WgsMetrics, Integer> output = new MetricsFile<>();
        output.read(new FileReader(outfile));

        final CollectWgsMetrics.WgsMetrics metrics = output.getMetrics().get(0);
        final CollectWgsMetrics.WgsMetrics nonZeroMetrics = output.getMetrics().get(1);


        // Some metrics should not change between with and without zero
        Assert.assertEquals(nonZeroMetrics.PCT_EXC_BASEQ, metrics.PCT_EXC_BASEQ);
        Assert.assertEquals(nonZeroMetrics.PCT_EXC_CAPPED, metrics.PCT_EXC_CAPPED);

        // Other metrics change when we ignore the zero depth bin
        Assert.assertEquals(nonZeroMetrics.GENOME_TERRITORY, 20);
        Assert.assertEquals(nonZeroMetrics.MEAN_COVERAGE, 3.0);
    }

    // TODO: put this in a util class to be shared with CollectWgsMetricsTest.java
    private SAMRecordSetBuilder createTestSAMBuilder(final File reference){
        final SAMFileHeader header = new SAMFileHeader();

        //Check that dictionary file is readable and then set header dictionary
        try {
            header.setSequenceDictionary(SAMSequenceDictionaryExtractor.extractDictionary(reference));
            header.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        } catch (final SAMException e) {
            e.printStackTrace();
        }

        //Set readGroupRecord
        final SAMReadGroupRecord readGroupRecord = new SAMReadGroupRecord(READ_GROUP_ID);
        readGroupRecord.setSample(SAMPLE);
        readGroupRecord.setPlatform(PLATFORM);
        readGroupRecord.setLibrary(LIBRARY);
        readGroupRecord.setPlatformUnit(READ_GROUP_ID);
        header.addReadGroup(readGroupRecord);

        final SAMRecordSetBuilder setBuilder = new SAMRecordSetBuilder(true, SAMFileHeader.SortOrder.coordinate, true, 100);
        setBuilder.setReadGroup(readGroupRecord);
        setBuilder.setUseNmFlag(true);
        setBuilder.setHeader(header);

        return(setBuilder);
    }


}
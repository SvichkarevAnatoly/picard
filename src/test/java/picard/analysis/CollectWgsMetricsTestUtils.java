package picard.analysis;


import htsjdk.samtools.*;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

import java.io.File;
import java.io.IOError;
import java.io.IOException;


/**
 * Contains util methods for CollectWgsMetricsTest, CollectWgsMetricsWithNonZeroCoverageTest
 */

public class CollectWgsMetricsTestUtils {

    /**
     *
     * @param reference
     * @param readGroupId
     * @param sample
     * @param platform
     * @param library
     * @return
     */
    protected static SAMRecordSetBuilder createTestSAMBuilder(final File reference,
                                                              final String readGroupId,
                                                              final String sample,
                                                              final String platform,
                                                              final String library){
        final SAMFileHeader header = new SAMFileHeader();

        //Check that dictionary file is readable and then set header dictionary
        try {
            header.setSequenceDictionary(SAMSequenceDictionaryExtractor.extractDictionary(reference));
            header.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        } catch (final SAMException e) {
            e.printStackTrace();
        }

        //Set readGroupRecord
        final SAMReadGroupRecord readGroupRecord = new SAMReadGroupRecord(readGroupId);
        readGroupRecord.setSample(sample);
        readGroupRecord.setPlatform(platform);
        readGroupRecord.setLibrary(library);
        readGroupRecord.setPlatformUnit(readGroupId);
        header.addReadGroup(readGroupRecord);

        final SAMRecordSetBuilder setBuilder = new SAMRecordSetBuilder(true, SAMFileHeader.SortOrder.coordinate, true, 100);
        setBuilder.setReadGroup(readGroupRecord);
        setBuilder.setUseNmFlag(true);
        setBuilder.setHeader(header);

        return(setBuilder);
    }


    /**
     * Template code for creating a custom SAM file for testing. Modify to suit your needs.
     *
     */
    private static void createTestSAM() throws IOException {
        final File testDir = new File("testdata/picard/analysis/directed/CollectHsMetrics/");
        final File reference = new File("testdata/picard/quality/chrM.reference.fasta");
        final String readGroupId = "TestReadGroup";
        final String sample = "TestSample";
        final String platform = "Illumina";
        final String library = "TestLibrary";
        final int numReads = 1;

        File samFile = File.createTempFile("TestSam", ".bam", testDir);
        samFile.deleteOnExit();

        final SAMRecordSetBuilder setBuilder = createTestSAMBuilder(reference, readGroupId, sample, platform, library);
        setBuilder.setReadLength(100);

        for (int i = 0; i < numReads; i++){
            setBuilder.addPair("MediocreBaseQ" + i, 0, 1, 200, false, false, "100M", "100M", false, true, 20);
        }

        for (int i = 0; i < numReads; i++){
            setBuilder.addPair("PoorBaseQ" + i, 0, 1, 200, false, false, "100M", "100M", false, true, 2);
        }

        final SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(setBuilder.getHeader(), false, samFile);

        for (final SAMRecord record : setBuilder) {
            writer.addAlignment(record);
        }

        writer.close();
        String message = "Place a breakpoint here";
    }

    /**
     *
     */
    public static void main(String[] args){
        try { createTestSAM(); } catch(IOException e) { ; }
    }

}

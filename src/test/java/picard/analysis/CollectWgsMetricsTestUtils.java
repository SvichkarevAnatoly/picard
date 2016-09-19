package picard.analysis;


import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecordSetBuilder;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

import java.io.File;

public class CollectWgsMetricsTestUtils {
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

}

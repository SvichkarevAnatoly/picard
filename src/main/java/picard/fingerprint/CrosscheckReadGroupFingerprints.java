
/*
 * The MIT License
 *
 * Copyright (c) 2010 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package picard.fingerprint;

import htsjdk.samtools.BamFileIoUtils;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.Fingerprinting;

import java.io.File;
import java.io.OutputStreamWriter;
import java.util.*;
import java.util.concurrent.TimeUnit;
import java.util.function.Function;
import java.util.stream.Collectors;

import static java.util.stream.Collectors.toSet;
import static picard.fingerprint.CrosscheckMetric.FingerprintResult.*;
import static picard.fingerprint.CrosscheckMetric.FingerprintResult;
import static picard.fingerprint.FingerprintChecker.*;

/**
 * Program to check that all (read-)groups within the set of BAM files appear to come from the same
 * individual. Can br used to cross-check libraries or samples too.
 *
 * @author Tim Fennell
 * @author Yossi Farjoun
 */
@CommandLineProgramProperties(
        usage = "Checks if all (read-)groups within a set of BAM files appear to come from the same individual.",
        usageShort = "Checks if all (read-)groups appear to come from the same individual.",
        programGroup = Fingerprinting.class
)
public class CrosscheckReadGroupFingerprints extends CommandLineProgram {

    @Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME,
            doc = "One or more input BAM files (or lists of BAM files) to compare fingerprints for.")
    public List<File> INPUT;

    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, optional = true,
            doc = "Optional output file to write metrics to. Default is to write to stdout.")
    public File OUTPUT;

    @Option(shortName = "H", doc = "The file of haplotype data to use to pick SNPs to fingerprint")
    public File HAPLOTYPE_MAP;

    @Option(shortName = "LOD",
            doc = "If any two groups (with the same sample name) match with a LOD score lower than the threshold " +
                    "the program will exit with a non-zero code to indicate error." +
                    " Program will also exit with an error if it finds two groups with different sample name that" +
                    "match with a LOD score greater than -LOD_THRESHOLD." +
                    "\n\n" +
                    "LOD score 0 means equal likelihood" +
                    "that the groups match vs. come from different individuals, negative LOD scores mean N logs more likely " +
                    "that the groups are from different individuals, and positive numbers mean N logs more likely that" +
                    " the groups are from the sample individual. ")
    public double LOD_THRESHOLD = 0;

    @Option(doc = "Roll fingerprints up to the library level and compare libraries against each other", mutex = {"CROSSCHECK_SAMPLES"})
    public boolean CROSSCHECK_LIBRARIES = false;

    @Option(doc = "Roll fingerprints up to the sample level and compare samples against each other.", mutex = {"CROSSCHECK_LIBRARIES"})
    public boolean CROSSCHECK_SAMPLES = false;

    @Option(doc = "The number of threads to use to process BAM files and generate Fingerprints.")
    public int NUM_THREADS = 1;

    @Option(doc = "Allow the use of duplicate reads in performing the comparison. Can be useful when duplicate " +
            "marking has been overly aggressive and coverage is low.")
    public boolean ALLOW_DUPLICATE_READS = false;

    @Option(doc = "If true then only read groups that do not relate to each other as expected will have their LODs reported.")
    public boolean OUTPUT_ERRORS_ONLY = false;

    @Option(doc = "The rate at which a heterozygous genotype in a normal sample turns into a homozygous (via loss of heterozygosity) " +
            "in the tumor (model assumes independent events, so this needs to be larger than reality).", optional = true)
    public double LOSS_OF_HET_RATE = 0.5;

    @Option(doc = "Expect all groups' fingerprints to match, irrespective of their sample names.  By default (with this value set to " +
            "false), groups (readgroups, libraries or samples) with different sample names are expected to mismatch, and those with the " +
            "same sample name are expected to match. ")
    public boolean EXPECT_ALL_GROUPS_TO_MATCH = false;

    @Option(doc = "When one or more mismatches between read groups are detected, exit with this value instead of 0.")
    public int EXIT_CODE_WHEN_MISMATCH = 1;

    private final Log log = Log.getInstance(CrosscheckReadGroupFingerprints.class);

    // the rate at which we expect a given genotype to change for an individual (exceedingly small)
    final static private double GENOTYPING_ERROR_RATE = 1e-6;

    @Override
    protected int doWork() {
        // Check inputs
        INPUT.forEach(IOUtil::assertFileIsReadable);
        IOUtil.assertFileIsReadable(HAPLOTYPE_MAP);
        if (OUTPUT != null) IOUtil.assertFileIsWritable(OUTPUT);

        final HaplotypeMap map = new HaplotypeMap(HAPLOTYPE_MAP);
        final FingerprintChecker checker = new FingerprintChecker(map);

        checker.setAllowDuplicateReads(ALLOW_DUPLICATE_READS);

        log.info("Done checking input files, moving onto fingerprinting files.");

        List<File> unrolledFiles = IOUtil.unrollFiles(INPUT, BamFileIoUtils.BAM_FILE_EXTENSION, IOUtil.SAM_FILE_EXTENSION);
        final Map<FingerprintIdDetails, Fingerprint> fpMap = checker.fingerprintSamFiles(unrolledFiles, NUM_THREADS, 1, TimeUnit.DAYS);

        log.info("Finished generating fingerprints from BAM files, moving on to cross-checking.");

        final List<CrosscheckMetric> metrics = new ArrayList<>();

        final int numUnexpected;
        if (CROSSCHECK_SAMPLES) {
            numUnexpected = crossCheckGrouped(fpMap, metrics, entry -> entry.getKey().sample, CrosscheckMetric.DataType.SAMPLE);
        } else if (CROSSCHECK_LIBRARIES) {
            numUnexpected = crossCheckGrouped(fpMap, metrics, entry -> entry.getKey().sample + "::" + entry.getKey().library, CrosscheckMetric.DataType.LIBRARY);
        } else {
            numUnexpected = crossCheckFingerprints(fpMap, CrosscheckMetric.DataType.READGROUP, metrics);
        }

        MetricsFile<CrosscheckMetric, ?> metricsFile = getMetricsFile();

        metricsFile.addAllMetrics(metrics);
        if (OUTPUT != null) {
            metricsFile.write(OUTPUT);
        } else {
            metricsFile.write(new OutputStreamWriter(System.out));
        }

        if (numUnexpected > 0) {
            log.info("WARNING: At least two read groups did not relate as expected.");
            return EXIT_CODE_WHEN_MISMATCH;
        } else {
            log.info("All read groups related as expected.");
            return 0;
        }
    }

    /**
     * Method that combines the fingerprint evidence across all the read groups for the same sample
     * and then produces a matrix of LOD scores for comparing every sample with every other sample.
     */
    private int crossCheckGrouped(final Map<FingerprintIdDetails, Fingerprint> fingerprints, final List<CrosscheckMetric> metrics,
                                  final Function<Map.Entry<FingerprintIdDetails, Fingerprint>, String> by,
                                  CrosscheckMetric.DataType type) {

        //collect the various entries according to the grouping "by"
        final Map<String, List<Map.Entry<FingerprintIdDetails, Fingerprint>>> collection =
                fingerprints.entrySet()
                        .stream()
                        .collect(Collectors.groupingBy(by));

        Map<FingerprintIdDetails, Fingerprint> fingerprintsByLibrary = collection.entrySet().stream()
                .collect(Collectors.toMap(
                        entry -> {
        //merge the keys (unequal values are eliminated).

                            final FingerprintIdDetails finalId = new FingerprintIdDetails();
                            entry.getValue().stream().forEach(id -> finalId.merge(id.getKey()));
                            return finalId;

                        }, entry -> {
        //merge the values by merging the fingerprints.

                            final FingerprintIdDetails firstDetail = entry.getValue().get(0).getKey();
                            //use the "by" function to determine the "info" part of the fingerprint
                            final Fingerprint sampleFp = new Fingerprint(firstDetail.sample, null, by.apply(entry.getValue().get(0)));
                            entry.getValue().stream().map(Map.Entry::getValue).collect(toSet()).forEach(sampleFp::merge);
                            return sampleFp;

                        }));

        return crossCheckFingerprints(fingerprintsByLibrary, type, metrics);
    }

    /**
     * Method that pairwise checks every pair of read groups and reports a LOD score for the two read groups
     * coming from the same sample.
     */
    private int crossCheckFingerprints(final Map<FingerprintIdDetails, Fingerprint> fingerprints, final CrosscheckMetric.DataType type, final List<CrosscheckMetric> metrics) {
        int unexpectedResults = 0;

        final List<FingerprintIdDetails> fingerprintIdDetails = new ArrayList<>(fingerprints.keySet());

        for (int i = 0; i < fingerprintIdDetails.size(); i++) {
            final FingerprintIdDetails lhsRg = fingerprintIdDetails.get(i);
            for (int j = i + 1; j < fingerprintIdDetails.size(); j++) {
                final FingerprintIdDetails rhsRg = fingerprintIdDetails.get(j);
                final boolean expectedToMatch = EXPECT_ALL_GROUPS_TO_MATCH || lhsRg.sample.equals(rhsRg.sample);

                final MatchResults results = FingerprintChecker.calculateMatchResults(fingerprints.get(lhsRg), fingerprints.get(rhsRg), GENOTYPING_ERROR_RATE, LOSS_OF_HET_RATE);

                final CrosscheckMetric.FingerprintResult result = getMatchResults(expectedToMatch, results);

                if (!OUTPUT_ERRORS_ONLY || !result.isExpected()) {
                    metrics.add(getMatchDetails(result, results, lhsRg, rhsRg, type));
                }
                if (result != INCONCLUSIVE && !result.isExpected()) unexpectedResults++;
            }
        }

        return unexpectedResults;
    }

    /**
     * Generates tab delimited string containing details about a possible match between fingerprints on two different SAMReadGroupRecords
     *
     * @param matchResult    String describing the match type.
     * @param results        MatchResults object
     * @param leftPuDetails  left hand side FingerprintIdDetails
     * @param rightPuDetails right hand side FingerprintIdDetails
     * @return tab delimited string containing details about a possible match
     */
    private CrosscheckMetric getMatchDetails(final FingerprintResult matchResult,
                                             final MatchResults results,
                                             final FingerprintIdDetails leftPuDetails,
                                             final FingerprintIdDetails rightPuDetails,
                                             final CrosscheckMetric.DataType type) {
        final CrosscheckMetric metric = new CrosscheckMetric();

        metric.RESULT = matchResult;
        metric.LOD_SCORE = results.getLOD();
        metric.LOD_SCORE_TUMOR_NORMAL = results.getLodTN();
        metric.LOD_SCORE_NORMAL_TUMOR = results.getLodNT();
        metric.DATA_TYPE = type;

        metric.LEFT_RUN_BARCODE = leftPuDetails.molecularBarcode;
        metric.LEFT_LANE = leftPuDetails.runLane;
        metric.LEFT_MOLECULAR_BARCODE_SEQUENCE = leftPuDetails.runBarcode;
        metric.LEFT_LIBRARY = leftPuDetails.library;
        metric.LEFT_SAMPLE = leftPuDetails.sample;

        metric.RIGHT_RUN_BARCODE = rightPuDetails.runBarcode;
        metric.RIGHT_LANE = rightPuDetails.runLane;
        metric.RIGHT_MOLECULAR_BARCODE_SEQUENCE = rightPuDetails.runBarcode;
        metric.RIGHT_LIBRARY = rightPuDetails.library;
        metric.RIGHT_SAMPLE = rightPuDetails.sample;

        return metric;
    }

    private FingerprintResult getMatchResults(boolean expectedToMatch, MatchResults results) {
        if (expectedToMatch) {
            if (results.getLOD() < LOD_THRESHOLD) {
                return UNEXPECTED_MISMATCH;
            } else if (results.getLOD() > -LOD_THRESHOLD) {
                return EXPECTED_MATCH;
            } else {
                return INCONCLUSIVE;
            }
        } else {
            if (results.getLOD() > -LOD_THRESHOLD) {
                return UNEXPECTED_MATCH;
            } else if (results.getLOD() < LOD_THRESHOLD) {
                return EXPECTED_MISMATCH;
            } else {
                return INCONCLUSIVE;
            }
        }
    }
}
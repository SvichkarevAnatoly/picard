package picard.fingerprint;

import htsjdk.samtools.metrics.MetricBase;

/**
 * A Metric class to hold the result of crosschecking fingerprints.
 * The same metric will be used for crosschecking Readgroups, libraries or samples.
 */
public class CrosscheckMetric extends MetricBase {

    // An enum representing whether the result of the fingerprinting was expected and whether it was a match.
    public enum FingerprintResult {
        EXPECTED_MATCH(true, true),
        EXPECTED_MISMATCH(true, false),
        UNEXPECTED_MATCH(false, true),
        UNEXPECTED_MISMATCH(false, false),
        INCONCLUSIVE(null, null);

        final private Boolean isExpected;
        final private Boolean isMatch;

        FingerprintResult(Boolean isExpected, Boolean isMatch){
            this.isExpected = isExpected;
            this.isMatch = isMatch;
        }

        public Boolean isExpected() {
            return isExpected;
        }

        public Boolean isMatch() {
            return isMatch;
        }
    }

    public enum DataType {
        SAMPLE,
        LIBRARY,
        READGROUP
    }

    public FingerprintResult RESULT;
    public DataType DATA_TYPE;
    public Double LOD_SCORE;
    public Double LOD_SCORE_TUMOR_NORMAL;
    public Double LOD_SCORE_NORMAL_TUMOR;
    public String LEFT_RUN_BARCODE="";
    public Integer LEFT_LANE = -1;
    public String LEFT_MOLECULAR_BARCODE_SEQUENCE ="";
    public String LEFT_LIBRARY = "";
    public String LEFT_SAMPLE;
    public String RIGHT_RUN_BARCODE ="";
    public Integer RIGHT_LANE = -1;
    public String RIGHT_MOLECULAR_BARCODE_SEQUENCE ="";
    public String RIGHT_LIBRARY ="";
    public String RIGHT_SAMPLE;


}

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

import htsjdk.samtools.metrics.MetricBase;

/**
 * A Metric class to hold the result of crosschecking fingerprints.
 * The same metric will be used for crosschecking Readgroups, libraries or samples.
 *
 * @author Yossi Farjoun
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

        FingerprintResult(Boolean isExpected, Boolean isMatch) {
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
    public String LEFT_RUN_BARCODE = "";
    public Integer LEFT_LANE = -1;
    public String LEFT_MOLECULAR_BARCODE_SEQUENCE = "";
    public String LEFT_LIBRARY = "";
    public String LEFT_SAMPLE;
    public String RIGHT_RUN_BARCODE = "";
    public Integer RIGHT_LANE = -1;
    public String RIGHT_MOLECULAR_BARCODE_SEQUENCE = "";
    public String RIGHT_LIBRARY = "";
    public String RIGHT_SAMPLE;
}

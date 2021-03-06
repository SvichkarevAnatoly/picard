/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
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

package picard.analysis.directed;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.IntervalList;
import picard.analysis.MetricAccumulationLevel;
import picard.sam.DuplicationMetrics;

import java.io.File;
import java.util.List;
import java.util.Set;

/**
 * Calculates HS metrics for a given SAM or BAM file. Requires the input of a list of
 * target intervals and a list of bait intervals. Can be invoked either on an entire
 * iterator of SAMRecords or be passed SAMRecords one at a time.
 *
 * @author Jonathan Burke
 */
public class HsMetricCollector extends TargetMetricsCollector<HsMetrics> {

    public HsMetricCollector(final Set<MetricAccumulationLevel> accumulationLevels,
                             final List<SAMReadGroupRecord> samRgRecords,
                             final ReferenceSequenceFile refFile,
                             final File perTargetCoverage,
                             final File perBaseCoverage,
                             final IntervalList targetIntervals,
                             final IntervalList probeIntervals,
                             final String probeSetName,
                             final int nearProbeDistance,
                             final int minimumMappingQuality,
                             final int minimumBaseQuality,
                             final boolean clipOverlappingReads,
                             final int coverageCap,
                             final int sampleSize) {
        super(accumulationLevels, samRgRecords, refFile, perTargetCoverage, perBaseCoverage, targetIntervals, probeIntervals, probeSetName, nearProbeDistance, minimumMappingQuality, minimumBaseQuality, clipOverlappingReads, coverageCap, sampleSize);
    }

    public HsMetricCollector(final Set<MetricAccumulationLevel> accumulationLevels,
                             final List<SAMReadGroupRecord> samRgRecords,
                             final ReferenceSequenceFile refFile,
                             final File perTargetCoverage,
                             final File perBaseCoverage,
                             final IntervalList targetIntervals,
                             final IntervalList probeIntervals,
                             final String probeSetName,
                             final int nearProbeDistance,
                             final int minimumMappingQuality,
                             final int minimumBaseQuality,
                             final boolean clipOverlappingReads,
                             final boolean noSideEffects,
                             final int coverageCap,
                             final int sampleSize) {
        super(accumulationLevels, samRgRecords, refFile, perTargetCoverage, perBaseCoverage, targetIntervals, probeIntervals, probeSetName, nearProbeDistance, minimumMappingQuality, minimumBaseQuality, clipOverlappingReads, noSideEffects, coverageCap, sampleSize);
    }

    @Override
    public HsMetrics convertMetric(final TargetMetrics targetMetrics) {
        final HsMetrics hsMetrics = new HsMetrics();
        TargetMetricsCollector.reflectiveCopy(targetMetrics, hsMetrics,
                new String[]{"PROBE_SET", "PROBE_TERRITORY", "ON_PROBE_BASES", "NEAR_PROBE_BASES", "OFF_PROBE_BASES", "PCT_OFF_PROBE", "ON_PROBE_VS_SELECTED", "MEAN_PROBE_COVERAGE"},
                new String[]{"BAIT_SET",  "BAIT_TERRITORY",  "ON_BAIT_BASES",  "NEAR_BAIT_BASES",  "OFF_BAIT_BASES",  "PCT_OFF_BAIT",  "ON_BAIT_VS_SELECTED",  "MEAN_BAIT_COVERAGE"}
        );

        hsMetrics.BAIT_DESIGN_EFFICIENCY = (double) hsMetrics.TARGET_TERRITORY / (double) hsMetrics.BAIT_TERRITORY;
        hsMetrics.PCT_USABLE_BASES_ON_BAIT   = hsMetrics.ON_BAIT_BASES   / (double) targetMetrics.PF_BASES;
        hsMetrics.PCT_USABLE_BASES_ON_TARGET = hsMetrics.ON_TARGET_BASES / (double) targetMetrics.PF_BASES;
        hsMetrics.HS_LIBRARY_SIZE = DuplicationMetrics.estimateLibrarySize(targetMetrics.PF_SELECTED_PAIRS, targetMetrics.PF_SELECTED_UNIQUE_PAIRS);

        //need HSLIBRARY_SIZE
        hsMetrics.HS_PENALTY_10X = calculateHsPenalty(hsMetrics.HS_LIBRARY_SIZE, targetMetrics, 10);
        hsMetrics.HS_PENALTY_20X = calculateHsPenalty(hsMetrics.HS_LIBRARY_SIZE, targetMetrics, 20);
	    hsMetrics.HS_PENALTY_30X = calculateHsPenalty(hsMetrics.HS_LIBRARY_SIZE, targetMetrics, 30);
	    hsMetrics.HS_PENALTY_40X = calculateHsPenalty(hsMetrics.HS_LIBRARY_SIZE, targetMetrics, 40);
	    hsMetrics.HS_PENALTY_50X = calculateHsPenalty(hsMetrics.HS_LIBRARY_SIZE, targetMetrics, 50);
	    hsMetrics.HS_PENALTY_100X = calculateHsPenalty(hsMetrics.HS_LIBRARY_SIZE, targetMetrics, 100);
        return hsMetrics;
    }

    /**
     * Attempts to calculate the HS penalty incurred by the library in order to get 80%
     * of target bases (in non-zero-covered targets) to a specific target coverage (e.g. 20X).
     *
     * @param coverageGoal the desired coverage target (e.g. 20X)
     * @return the hs penalty - a multiplier that tells if you want, e.g. 20X coverage, then you will
     *         need to produce this many PF aligned bases per target bases in your design!
     */
    private double calculateHsPenalty(final Long librarySize, final TargetMetrics targetMetrics, final int coverageGoal) {
        if (librarySize == null) return 0;

        final double meanCoverage  = targetMetrics.ON_TARGET_FROM_PAIR_BASES / (double) targetMetrics.TARGET_TERRITORY;
        final double fold80        = targetMetrics.FOLD_80_BASE_PENALTY;
        final long pairs           = targetMetrics.PF_SELECTED_PAIRS;
        final long uniquePairs     = targetMetrics.PF_SELECTED_UNIQUE_PAIRS;
        final double onTargetPct   = (double) targetMetrics.ON_TARGET_BASES / (double) targetMetrics.PF_UQ_BASES_ALIGNED;

        final double uniquePairGoalMultiplier = (coverageGoal / meanCoverage) * fold80;
        double pairMultiplier = uniquePairGoalMultiplier;
        double increment = 1;
        boolean goingUp = uniquePairGoalMultiplier >= 1;
        double finalPairMultiplier = -1;

        // Converge "pairMultiplier" to the number that gives us a uniquePairMultiplier equal
        // to the coverage multiplier we desire.  If we can't get there with 1000X coverage,
        // we're not going to get there!
        for (int i=0; i<10000; ++i) {
            final double uniquePairMultiplier = DuplicationMetrics.estimateRoi(librarySize, pairMultiplier, pairs, uniquePairs);

            if (Math.abs(uniquePairMultiplier - uniquePairGoalMultiplier) / uniquePairGoalMultiplier <= 0.001) {
                finalPairMultiplier  = pairMultiplier;
                break;
            }
            else if ((uniquePairMultiplier > uniquePairGoalMultiplier && goingUp) ||
                    (uniquePairMultiplier < uniquePairGoalMultiplier && !goingUp)){
                increment /= 2;
                goingUp = !goingUp;
            }

            pairMultiplier += (goingUp ? increment : -increment);
        }

        if (finalPairMultiplier == -1) {
            return -1;
        }
        else {
            final double uniqueFraction = (uniquePairs * uniquePairGoalMultiplier) / (pairs * finalPairMultiplier);
            return (1 / uniqueFraction) * fold80 * (1 / onTargetPct);
        }
    }

}

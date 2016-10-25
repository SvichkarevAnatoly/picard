package picard.sam.markduplicates.util;

import htsjdk.samtools.util.SequenceUtil;

public class ReadSequence {
    byte[] read;
    long hash;

    public byte[] getRead() {
        return read;
    }

    public long getLongHashCode() {
        return hash;
    }

    public void setRead(byte[] read) {
        this.read = read;
    }

    public void initHashCode(final int seedLength) {
        hash = SequenceUtil.getReadBaseHashCode(read[0]);
        for (int i = 1; i < seedLength; i++) {
            hash = (hash << 4) + SequenceUtil.getReadBaseHashCode(read[i]);
        }
    }

    public void setHash(long hash) {
        this.hash = hash;
    }
}

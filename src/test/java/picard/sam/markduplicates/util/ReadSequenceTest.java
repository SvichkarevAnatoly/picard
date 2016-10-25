package picard.sam.markduplicates.util;

import org.testng.Assert;
import org.testng.annotations.Test;

public class ReadSequenceTest {
    public static final int SEED_LENGTH = 5;

    @Test
    public void hashOneReadBaseACGT() {
        final ReadSequence rs = new ReadSequence();
        rs.setRead("A".getBytes());
        rs.initHashCode(1);
        Assert.assertEquals(rs.getLongHashCode(), 1);

        rs.setRead("C".getBytes());
        rs.initHashCode(1);
        Assert.assertEquals(rs.getLongHashCode(), 2);

        rs.setRead("G".getBytes());
        rs.initHashCode(1);
        Assert.assertEquals(rs.getLongHashCode(), 4);

        rs.setRead("T".getBytes());
        rs.initHashCode(1);
        Assert.assertEquals(rs.getLongHashCode(), 8);
    }

    @Test
    public void hashOneReadBaseLowerCase_acgt() {
        final ReadSequence rs = new ReadSequence();
        rs.setRead("a".getBytes());
        rs.initHashCode(1);
        Assert.assertEquals(rs.getLongHashCode(), 1);

        rs.setRead("c".getBytes());
        rs.initHashCode(1);
        Assert.assertEquals(rs.getLongHashCode(), 2);

        rs.setRead("g".getBytes());
        rs.initHashCode(1);
        Assert.assertEquals(rs.getLongHashCode(), 4);

        rs.setRead("t".getBytes());
        rs.initHashCode(1);
        Assert.assertEquals(rs.getLongHashCode(), 8);
    }

    @Test
    public void hashOneReadBaseIUPAC() {
        final ReadSequence rs = new ReadSequence();
        rs.setRead("M".getBytes());
        rs.initHashCode(1);
        Assert.assertEquals(rs.getLongHashCode(), 0b11); // A & C
        Assert.assertEquals(rs.getLongHashCode(), 3); // 0b11 == 3

        rs.setRead("V".getBytes());
        rs.initHashCode(1);
        Assert.assertEquals(rs.getLongHashCode(), 0b111); // A & C & G
    }
}

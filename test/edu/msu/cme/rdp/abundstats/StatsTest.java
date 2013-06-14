/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package edu.msu.cme.rdp.abundstats;

import edu.msu.cme.rdp.abundstats.Sample;
import edu.msu.cme.rdp.abundstats.calculators.Sorensen;
import edu.msu.cme.rdp.abundstats.calculators.Jaccard;
import edu.msu.cme.rdp.abundstats.AbundStatsCalculator.UV;
import java.util.Arrays;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author fishjord
 */
public class StatsTest {


    private static final Sample s1 = new Sample("Sample1", Arrays.asList(21, 20, 43,  0,  23, 11,  1, 92, 0, 12, 15,  1,  2, 21, 5,  2));
    private static final Sample s2 = new Sample("Sample2", Arrays.asList( 0, 40, 20, 30, 100, 19, 15,  1, 1,  1,  1, 25, 15,  2, 2, 16));

    private static final float expectedJaccard = 0.9615f;
    private static final float expectedSorensen = 0.9804f;
    private static final float expectedU = 1.0000f;
    private static final float expectedV = 0.9615f;

    public StatsTest() {
    }

    /**
     * Test of calculate method, of class Jaccard.
     */
    @Test
    public void testCalcJaccard() {
        assertEquals("Jaccard index doesn't match expected!", new Jaccard().calculate(s1, s2), expectedJaccard, 0.0001);
        assertEquals("Sorensen index doesn't match expected!", new Sorensen().calculate(s1, s2), expectedSorensen, 0.0001);

        UV uv = AbundStatsCalculator.calculateUV(s1, s2);
        assertEquals("U doesn't match expected", uv.U, expectedU, 0.0001);
        assertEquals("V doesn't match expected", uv.V, expectedV, 0.0001);
    }

}
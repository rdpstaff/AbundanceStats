/*
 * Copyright (C) 2012 Michigan State University <rdpstaff at msu.edu>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.msu.cme.rdp.abundstats;

import java.io.PrintStream;
import java.util.List;

/**
 *
 * @author fishjord
 */
public class AbundStatUtils {

    public static float[][] computeSymMatrix(List<Sample> samples, AbundStatsCalculator calc) {
        float[][] ret = new float[samples.size()][samples.size()];

        for(int s1 = 0;s1 < samples.size();s1++) {
            ret[s1][s1] = 0;
            for(int s2 = s1 + 1;s2 < samples.size();s2++) {
                ret[s1][s2] = ret[s2][s1] = calc.calculate(samples.get(s1), samples.get(s2));
            }
        }

        return ret;
    }

    public static void writeMatrix(PrintStream out, List<Sample> samples, float[][] matrix) {

        out.print("\t");
        for(Sample s1 : samples) {
            out.print(s1.getSampleName() + "\t");
        }
        out.println();

        for(int s1 = 0;s1 < samples.size();s1++) {

            out.print(samples.get(s1).getSampleName() + "\t");
            for(int s2 = 0;s2 < samples.size();s2++) {

                out.print(matrix[s1][s2] + "\t");
            }
            out.println();
        }
    }
}

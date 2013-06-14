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

import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author fishjord
 */
public abstract class AbundStatsCalculator {

    protected static class UV {

        public float U = 0, V = 0;
    }

    public static UV calculateUV(Sample s1, Sample s2) {
        if (s1.getAbundances().size() != s2.getAbundances().size()) {
            throw new IllegalArgumentException("Abundance data in samples one and two must have the same number of species");
        }

        UV ret = new UV();

        List<Integer> X = s1.getAbundances();
        List<Integer> Y = s2.getAbundances();
        float f1Plus = 0; // f1+
        float f2Plus = 0; // f2+
        float fPlus1 = 0; // f+1
        float fPlus2 = 0; // f+2

        float n = 0;
        float m = 0;
        List<Integer> d12 = new ArrayList();

        //Combining a couple of steps
        //First check for D1, D2, and D12 (number of non 0 in sample 1, sample 2, and both respectively)
        //Additionally compute f(+)1/2(+) at the same time
        for (int index = 0; index < X.size(); index++) {
            int xi = X.get(index);
            int yi = Y.get(index);

            if (xi != 0 && yi != 0) {
                d12.add(index);

                if (xi == 1) {
                    f1Plus++;
                } else if (xi == 2) {
                    f2Plus++;
                }

                if (yi == 1) {
                    fPlus1++;
                } else if (yi == 2) {
                    fPlus2++;
                }
            }

            n += xi;
            m += yi;
        }


        if (f1Plus == 0) {
            f1Plus = 1;
        }
        if (f2Plus == 0) {
            f2Plus = 1;
        }
        if (fPlus1 == 0) {
            fPlus1 = 1;
        }
        if (fPlus2 == 0) {
            fPlus2 = 1;
        }


        float upt1 = 0;
        float upt2 = 0;
        for (int index : d12) {
            upt1 += X.get(index) / n;

            if (Y.get(index) == 1) {
                upt2 += X.get(index) / n;
            }
        }

        ret.U = upt1 + (((m - 1) / m) * (fPlus1 / (2 * fPlus2)) * upt2);
        ret.U = Math.min(ret.U, 1);

        float vpt1 = 0;
        float vpt2 = 0;
        for (int index : d12) {
            vpt1 += Y.get(index) / m;

            if (X.get(index) == 1) {
                vpt2 += Y.get(index) / m;
            }
        }

        ret.V = vpt1 + (((n - 1) / n) * (f1Plus / (2 * f2Plus)) * vpt2);
        ret.V = Math.min(ret.V, 1);

        //System.out.println("d1= " + n + " d2= " + m + " d12=" + d12.size());
        //System.out.println("f1+= " + f1Plus + " f2Plus= " + f2Plus + " f+1=" + fPlus1 + " f+2=" + fPlus2);
        //System.out.println("U= " + ret.U + " V=" + ret.V);

        return ret;
    }

    abstract public float calculate(Sample s1, Sample s2);

    abstract public String getCalcName();
}

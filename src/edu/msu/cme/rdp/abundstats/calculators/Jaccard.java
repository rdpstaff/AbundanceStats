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

package edu.msu.cme.rdp.abundstats.calculators;

import edu.msu.cme.rdp.abundstats.Sample;
import edu.msu.cme.rdp.abundstats.AbundStatsCalculator;

/**
 *
 * @author fishjord
 */
public class Jaccard extends AbundStatsCalculator{

    private boolean distance = false;
    
    public Jaccard() {}
    public Jaccard(boolean distance) { this.distance = distance; }

    public float calculate(Sample sample1, Sample sample2) {
        UV uv = calculateUV(sample1, sample2);

        if(uv.U == 0 && uv.V == 0) {
            return (distance)? 1 : 0;
        }

        float sim = Math.min((uv.U * uv.V) / (uv.U + uv.V - uv.U * uv.V), 1);
        return (distance)? 1 - sim : sim;
    }

    public String getCalcName() {
        return "jaccard";
    }
}

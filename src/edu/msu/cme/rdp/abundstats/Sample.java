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
public class Sample {
    private String sampleName;
    private List<Integer> abundances;

    public Sample(String sampleName, List<Integer> abundances) {
        this.sampleName = sampleName;
        this.abundances = abundances;
    }

    public Sample(String sampleName) {
        this.sampleName = sampleName;
        this.abundances = new ArrayList();
    }

    public List<Integer> getAbundances() {
        return abundances;
    }

    public String getSampleName() {
        return sampleName;
    }

    public void addSpecies(int abundence) {
        abundances.add(abundence);
    }
}

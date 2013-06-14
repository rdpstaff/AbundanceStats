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
package edu.msu.cme.rdp.rarefaction;

/**
 *
 * @author fishjord
 */
public class RarefactionResult {

    public double distance;
    public double[] eArray;
    public double[] vArray;

    public RarefactionResult(double distance, double[] eArray, double[] vArray) {
        this.distance = distance;
        this.eArray = eArray;
        this.vArray = vArray;
    }

    public double getDistance() {
        return distance;
    }

    public double[] geteArray() {
        return eArray;
    }

    public double[] getvArray() {
        return vArray;
    }
}

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

import edu.msu.cme.pyro.cluster.io.RDPClustParser.ClusterSample;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.io.File;
import java.io.PrintStream;
import java.io.IOException;

/**
 *
 * @author fishjord
 */
public class RarefactionWriter {

    public static final int MAX_PLOT_POINTS = 1000;
    public static final DecimalFormat dformat = new DecimalFormat("0.00");

    public static int[] getIndiciesToPlot(int numSeqs) {
        int[] ret = new int[MAX_PLOT_POINTS];

        int index = 0;
        for (int plotIndex = 0; plotIndex < numSeqs; plotIndex += (int) Math.ceil(numSeqs / (double) MAX_PLOT_POINTS)) {
            ret[index++] = plotIndex;
        }

        ret = Arrays.copyOfRange(ret, 0, index);

        return ret;
    }

    public static void writeRarefactionResults(Map<ClusterSample, List<RarefactionResult>> rarefactionResults, File outputDirectory) throws IOException {
        if (!outputDirectory.isDirectory()) {
            throw new IOException(outputDirectory + " isn't a directory!");
        }

        for (ClusterSample sample : rarefactionResults.keySet()) {
            writeRarefactionResult(rarefactionResults.get(sample), sample.getSeqs(), new File(outputDirectory, sample.getName() + "_rarefaction.txt"));
        }
    }

    public static void writeRarefactionResult(List<RarefactionResult> results, int numSeqs,  File outputFile) throws IOException {
        PrintStream out = new PrintStream(outputFile);

        for (RarefactionResult val : results) {
            out.print("\t" + dformat.format(val.distance));
        }
        for (RarefactionResult val : results) {
            out.print("\t" + dformat.format(val.distance) + "U\t" + dformat.format(val.distance) + "L");
        }
        out.print("\n");

        for (int plotIndex : getIndiciesToPlot(numSeqs)) {
            out.print(plotIndex);
            // output only the rarefaction values, easy for chart
            for (RarefactionResult val : results) {
                double E = val.eArray[plotIndex];
                out.print("\t" + dformat.format(E));
            }
            //output upper and lower
            for (RarefactionResult val : results) {
                double E = val.eArray[plotIndex];
                double V = val.vArray[plotIndex];
                double conf = Rarefaction.getConf(V);
                out.print("\t" + dformat.format(E + conf) + "\t" + dformat.format(E - conf));
            }
            out.println();
        }

        out.close();
    }
}

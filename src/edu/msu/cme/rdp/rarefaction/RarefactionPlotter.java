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

import java.util.Map;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.ArrayList;
import edu.msu.cme.pyro.cluster.io.RDPClustParser;
import edu.msu.cme.pyro.cluster.io.RDPClustParser.ClusterSample;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

/**
 *
 * @author fishjord
 */
public class RarefactionPlotter {

    private static final double[] interestedPoints = {0.001f, 0.031f, 0.051f, 0.101f};

    public static void plotRarefaction(Map<ClusterSample, List<RarefactionResult>> rarefactionResults, File outputDirectory) throws IOException {
        if (!outputDirectory.isDirectory()) {
            throw new IOException(outputDirectory + " isn't a directory!");
        }

        for (ClusterSample sample : rarefactionResults.keySet()) {
            List<RarefactionResult> resultsToPlot = new ArrayList();

            int currCutoff = 0;
            RarefactionResult lastResult = null;
            for (RarefactionResult result : rarefactionResults.get(sample)) {
                if (result.getDistance() <= interestedPoints[currCutoff]) {
                    lastResult = result;
                } else {
                    if (lastResult != null) {
                        resultsToPlot.add(lastResult);

                        if (++currCutoff >= interestedPoints.length) {
                            break;
                        }
                    }
                }
            }

            plotToFile(resultsToPlot, sample.getSeqs(), sample.getName() + " Rarefaction", new File(outputDirectory, sample.getName() + "_rarefaction.png"));
        }

    }

    private static void plotToFile(List<RarefactionResult> resultsToPlot, int numSeqs, String title, File outputFile) throws IOException {
        XYSeriesCollection dataset = new XYSeriesCollection();

        for (RarefactionResult result : resultsToPlot) {
            XYSeries series = new XYSeries(RarefactionWriter.dformat.format(result.distance) + " distance");

            for (int plotIndex : RarefactionWriter.getIndiciesToPlot(numSeqs)) {
                series.add(result.geteArray()[plotIndex], plotIndex);
            }

            dataset.addSeries(series);
        }

        JFreeChart chart = ChartFactory.createXYLineChart(title, "", "", dataset, PlotOrientation.HORIZONTAL, true, false, false);
        ChartUtilities.saveChartAsPNG(outputFile, chart, 1280, 1024);
    }

    public static void main(String[] args) throws Exception {

        args = new String[]{"/scratch/test_files/rarefaction/FUSION26.comp.clust", "rarefaction_test.txt"};

        RDPClustParser parser = new RDPClustParser(new File(args[0]));

        Map<ClusterSample, List<RarefactionResult>> resultsMap = Rarefaction.doRarefaction(parser);

        plotRarefaction(resultsMap, new File("."));
    }
}

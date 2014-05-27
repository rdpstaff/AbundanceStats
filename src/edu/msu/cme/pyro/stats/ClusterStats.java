/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.msu.cme.pyro.stats;

import edu.msu.cme.pyro.cluster.io.ClusterToBiom;
import edu.msu.cme.pyro.cluster.io.RDPClustParser;
import edu.msu.cme.pyro.cluster.io.RDPClustParser.ClusterSample;
import edu.msu.cme.pyro.cluster.io.RDPClustParser.Cutoff;
import edu.msu.cme.pyro.cluster.utils.Cluster;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

/**
 *
 * @author fishjord
 */
public class ClusterStats {

    public static void writeStats(File statsFile, File statsImg, Map<String, Integer> otuMap, int totalSeqs) throws IOException {
        PrintStream out = new PrintStream(statsFile);
        out.println("Total sequences\t" + totalSeqs);
        out.println();

        out.println("Cutoff\tClusters(OTUs)");
        XYSeries cutoffToSeqSeries = new XYSeries("Sequences");
        for (String key : otuMap.keySet()) {
            out.println(key + "\t" + otuMap.get(key));
            cutoffToSeqSeries.add(otuMap.get(key), Float.valueOf(key));
        }
        out.close();

        JFreeChart chart = ChartFactory.createScatterPlot("Clusters vs Distance", "Number of Clusters", "Distance Cutoff", new XYSeriesCollection(cutoffToSeqSeries), PlotOrientation.HORIZONTAL, false, true, false);
        ChartUtilities.writeChartAsPNG(new PrintStream(statsImg), chart, 640, 640);
    }

    public static void writeArityStats(File statsFile, File statsImg, Map<String, Map<Integer, Integer>> cutoffArityMap, Map<String, Integer> numClusters, int numSamples) throws IOException {
        PrintStream out = new PrintStream(statsFile);
        DefaultCategoryDataset dataset = new DefaultCategoryDataset();

        out.print("cutoff");
        for (int arity = 1; arity <= numSamples; arity++) {
            out.print("\t" + arity);
        }
        out.println("\tnum_clusters");

        for (String cutoff : cutoffArityMap.keySet()) {
            out.print(cutoff);
            double totalClusters = numClusters.get(cutoff);
            for (Integer arity = 1; arity <= numSamples; arity++) {
                int clustCount = cutoffArityMap.get(cutoff).get(arity);
                out.print("\t" + clustCount);
                dataset.addValue((clustCount / totalClusters) * 100, arity, cutoff);
            }
            out.println((int) totalClusters);
        }
        out.close();

        JFreeChart chart = ChartFactory.createStackedBarChart("", "Distance Cutoff", "% of Total Clusters", dataset, PlotOrientation.HORIZONTAL, true, false, false);

        ChartUtilities.writeChartAsPNG(new PrintStream(statsImg), chart, 1280, 1280);
    }

    public static void writeStats(File outDir, List<File> clusterFiles) throws IOException {
        writeStats(outDir, clusterFiles, true);
    }
    
    public static void writeStats(File outDir, List<File> clusterFiles, boolean outputBiom) throws IOException {

        for (File clusterFile : clusterFiles) {
            RDPClustParser parser = new RDPClustParser(clusterFile);

            Map<ClusterSample, Map<String, Integer>> otuCounts = new LinkedHashMap();
            Cutoff cutoff;

            for (ClusterSample cs : parser.getClusterSamples()) {
                otuCounts.put(cs, new LinkedHashMap());
            }

            while ((cutoff = parser.readNextCutoff()) != null) {
                for (ClusterSample sample : parser.getClusterSamples()) {
                    int otus = 0;
                    for (Cluster c : cutoff.getClusters().get(sample.getName())) {
                        if (c.getNumberOfSeqs() > 0) {
                            otus++;
                        }
                    }
                    otuCounts.get(sample).put(cutoff.getCutoff(), otus);
                }
                // if outBiom is true
                if ( outputBiom){                
                    File outFile = new File(outDir, clusterFile.getName() + "_" + cutoff.getCutoff() + ".biom");
                    PrintStream out = new PrintStream(outFile);
                    try {
                        ClusterToBiom.writeCutoff(cutoff, out);
                    } finally {
                        out.close();
                    }                
                }
            }

            for (ClusterSample sample : parser.getClusterSamples()) {
                File statOutFile = new File(outDir, clusterFile.getName() + "_" + sample.getName() + "_otus.txt");
                File otuChartOutFile = new File(outDir, clusterFile.getName() + "_" + sample.getName() + "_otus.png");
                int totalSeqs = sample.getSeqs();
                ClusterStats.writeStats(statOutFile, otuChartOutFile, otuCounts.get(sample), totalSeqs);
            }
                    
            parser.close();
        }
    }

    

    public static void main(String[] args) throws Exception {
        if (args.length < 2) {
            System.out.println("Usage: ClusterStats <resultDir> <cluster_file>...");
            return;
        }

        File outDir = new File(args[0]);
        List<File> clustFiles = new ArrayList();
        for(int index = 1;index < args.length;index++) {
            clustFiles.add(new File(args[index]));
        }
        
        writeStats(outDir, clustFiles);
    }
}

package edu.msu.cme.rdp.rarefaction;

import edu.msu.cme.pyro.cluster.io.RDPClustParser;
import edu.msu.cme.pyro.cluster.io.RDPClustParser.ClusterSample;
import edu.msu.cme.pyro.cluster.io.RDPClustParser.Cutoff;
import edu.msu.cme.pyro.cluster.utils.Cluster;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.HashSet;
import java.util.Map;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.Iterator;

public class Rarefaction {

    /**
     * So this class is an extreme optimization on an inner for loop
     *
     * Originally this was representated by a Map<Integer, Double>, but it would result in
     * N^3 * 3 function calls in the calcV function per cutoff, so just think of this as a
     * map, where usedPValIndicies is a list of keys, and pvalues are the values (there can
     * and often are unused values)
     *
     * usedPValIndicies is a list of the cluster sizes and sums there of
     */
    private static class PValues {

        int[] usedPValIndices;
        double[] pvalues;
    }

    /**
     * This class like PValues is an extreme optimization on a Map<Integer, Integer>
     * in this case however there is no key/value list, instead there are two parallel arrays
     * so clusterSizes[i] is the size of an arbitrary cluster and clustersOfSize[i] is the
     * number of clusters of that size
     *
     * Also we keep track of the largest cluster, since we'll need it to initialize the pvalues
     * (largest number of possible pvals is 2 * largest cluster size + 1), and it's practically free
     * to calculate that value while mapping the cluster sizes
     */
    private static class ClusterCounts {

        int[] clusterSizes;
        int[] clustersOfSize;
        int largestCluster = 0;
    }

    public static Map<ClusterSample, List<RarefactionResult>> doRarefaction(RDPClustParser parser) throws IOException {

        Map<ClusterSample, List<RarefactionResult>> resultsMap = new HashMap();
        for (ClusterSample clusterSample : parser.getClusterSamples()) {
            resultsMap.put(clusterSample, new ArrayList());
        }

        Cutoff c;
        while((c = parser.readNextCutoff()) != null) {
            Map<ClusterSample, RarefactionResult> result = cutoffRarefaction(parser.getClusterSamples(), c);

            for (ClusterSample clusterSample : result.keySet()) {
                resultsMap.get(clusterSample).add(result.get(clusterSample));
            }
        }

        return resultsMap;
    }

    private static Map<ClusterSample, RarefactionResult> cutoffRarefaction(List<ClusterSample> samples, Cutoff cutoff) {

        double distance = Double.valueOf(cutoff.getCutoff());
        Map<ClusterSample, RarefactionResult> ret = new HashMap();

        for (ClusterSample sample : samples) {
            int numSeqsInSample = sample.getSeqs();

            double[] eValues = new double[numSeqsInSample + 1];
            double[] vValues = new double[numSeqsInSample + 1];

            List<Cluster> clusters = cutoff.getClusters().get(sample.getName());
            for (Iterator<Cluster> iterator = clusters.listIterator(); iterator.hasNext();) {
                Cluster c = iterator.next();
                if (c.getNumberOfSeqs() == 0) {
                    iterator.remove();
                }
            }

            ClusterCounts clusterSizeCounts = mapClusterSizes(clusters);
            PValues pvals = initializePValues(clusterSizeCounts);

            long p = 0, e = 0, v = 0;
            for (int i = 1; i <= numSeqsInSample; i++) {
                updatePValues(i, numSeqsInSample, pvals);
                eValues[i] = calcE(clusterSizeCounts, pvals);
                vValues[i] = calcV(clusterSizeCounts, pvals, clusters);
            }

            ret.put(sample, new RarefactionResult(distance, eValues, vValues));
        }

        return ret;
    }

    /**
     * This function maps the size of clusters to the numbers of clusters of that size
     *
     * The result is a hella optimized class
     *
     * @see ClsuterCounts
     *
     * @param oneSampleCutoffClusters
     * @return
     */
    private static ClusterCounts mapClusterSizes(List<Cluster> oneSampleCutoffClusters) {
        Map<Integer, Integer> clusterSizeCounts = new HashMap();

        ClusterCounts ret = new ClusterCounts();
        for (Cluster cluster : oneSampleCutoffClusters) {
            Integer count = clusterSizeCounts.get(cluster.getNumberOfSeqs());
            if (count == null) {
                count = 0;
            }

            if (cluster.getNumberOfSeqs() > ret.largestCluster) {
                ret.largestCluster = cluster.getNumberOfSeqs();
            }

            clusterSizeCounts.put(cluster.getNumberOfSeqs(), count + 1);
        }

        ret.clusterSizes = new int[clusterSizeCounts.size()];
        ret.clustersOfSize = new int[clusterSizeCounts.size()];

        int index = 0;
        for (Entry<Integer, Integer> entry : clusterSizeCounts.entrySet()) {
            ret.clusterSizes[index] = entry.getKey();
            ret.clustersOfSize[index] = entry.getValue();

            index++;
        }

        return ret;
    }

    private static PValues initializePValues(ClusterCounts clusterSizeCounts) {
        PValues ret = new PValues();

        ret.pvalues = new double[clusterSizeCounts.largestCluster * 2 + 1];
        Set<Integer> usedPValues = new HashSet();

        for (int i = 0; i < clusterSizeCounts.clusterSizes.length; i++) {
            int countI = clusterSizeCounts.clusterSizes[i];

            usedPValues.add(countI);
            ret.pvalues[countI] = 1f;
            for (int j = i; j < clusterSizeCounts.clusterSizes.length; j++) {
                int countK = clusterSizeCounts.clusterSizes[j] + countI;

                usedPValues.add(countK);
                ret.pvalues[countK] = 1f;
            }
        }

        ret.usedPValIndices = new int[usedPValues.size()];
        int index = 0;
        for (int usedPValue : usedPValues) {
            ret.usedPValIndices[index++] = usedPValue;
        }

        return ret;
    }

    private static void updatePValues(double n, int totalSeqCount, PValues pvals) {
        double denominator = 1.0d / (totalSeqCount - n + 1.0d);
        for (int H : pvals.usedPValIndices) {
            double pval = pvals.pvalues[H];
            pval *= (totalSeqCount - (double) H - n + 1.0d) * denominator;
            pvals.pvalues[H] = pval;
        }
    }

    private static double calcE(ClusterCounts clusterSizeToCounts, PValues pvals) {
        double E = 0;

        for (int index = 0; index < clusterSizeToCounts.clusterSizes.length; index++) {
            int clusterSize = clusterSizeToCounts.clusterSizes[index];
            int clustersOfSize = clusterSizeToCounts.clustersOfSize[index];
            E += clustersOfSize * (1.0d - pvals.pvalues[clusterSize]);
        }

        return E;
    }

    private static double calcV(ClusterCounts clusterSizeToCounts, PValues pvals, List<Cluster> clusters) {
        double V = 0;

        for (int i = 0; i < clusterSizeToCounts.clusterSizes.length; i++) {
            int sizeI = clusterSizeToCounts.clusterSizes[i];
            int countI = clusterSizeToCounts.clustersOfSize[i];

            for (int j = 0; j < clusterSizeToCounts.clusterSizes.length; j++) {
                int sizeJ = clusterSizeToCounts.clusterSizes[j];
                int countJ = clusterSizeToCounts.clustersOfSize[j];

                if (sizeI != sizeJ) {
                    V += countI * countJ * (pvals.pvalues[sizeI + sizeJ] - pvals.pvalues[sizeI] * pvals.pvalues[sizeJ]);
                } else if (sizeI > 1) { // use to remove possibility that count_i > 1/2N because pValues[>N] may go to infinity
                    V += countI * (countI - 1) * (pvals.pvalues[sizeI + sizeJ] - pvals.pvalues[sizeI] * pvals.pvalues[sizeJ]);
                }
            }
        }

        for (Cluster cluster : clusters) {
            double pval = pvals.pvalues[cluster.getNumberOfSeqs()];
            V += pval * (1 - pval);
        }

        return V;
    }

    public static double getConf(double V) {
        double conf = 0.0;
        if (V > 0) {  // can be very small negative due to rounding error
            conf = Math.sqrt(V) * 1.96d;
        }
        return conf;

    }

    public static void main(String[] args) throws IOException {
        String usage = "Usage:java Rarefaction in.clust output_dir [plot]";
        if (args.length != 2 && args.length != 3) {
            throw new IllegalArgumentException(usage);
        }

        File clusterFile = new File(args[0]);
        File outputDir = new File(args[1]);

        if (!outputDir.exists()) {
            if (!outputDir.mkdirs()) {
                throw new IOException("Failed to make output directory " + outputDir);
            }
        } else if (!outputDir.exists()) {
        }

        boolean plot = args.length == 3;

        RDPClustParser parser = new RDPClustParser(new File(args[0]));
        Map<ClusterSample, List<RarefactionResult>> resultsMap = doRarefaction(parser);
        parser.close();
        RarefactionWriter.writeRarefactionResults(resultsMap, outputDir);

        if (plot) {
            RarefactionPlotter.plotRarefaction(resultsMap, outputDir);
        }
    }
}

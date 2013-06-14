package edu.msu.cme.rdp.abundstats.cli;

import edu.msu.cme.pyro.cluster.io.RDPClustParser;
import edu.msu.cme.pyro.cluster.utils.Cluster;
import java.io.*;
import java.text.DecimalFormat;
import java.util.*;

/**
 *
 * @author siddiq15 class to produce "R" compatible format from the cluster
 * file.
 */
public class RFormatter {

    private static final DecimalFormat DECIMAL_FORMAT = new DecimalFormat("#.##");

    /**
     * This method writes data in tab delimited tabular format in the following
     * manner that is compatible to R software. output a matrix file with first
     * row contains the clusterNO and each row represent one sample The first
     * row contains the OTU names (tab delimited) The remaining rows contain the
     * count for each sample, tab delimited OTU_1 OTU_2 ... sample1 count count
     * ... sample2 count count ...
     *
     */
    public static boolean createTabulatedFormat(RDPClustParser.Cutoff cutoff, PrintStream stream) throws IOException {

        // find the clusters from all the samples in the same otu
        TreeSet<Integer> allOTUset = new TreeSet<Integer>();
        HashMap<String, HashMap<Integer, Integer>> clusterMap = new HashMap<String, HashMap<Integer, Integer>>(); // sampleName , HashMap (clusterID, count)
        int maxClustID = 0;
        for (String sample : cutoff.getClusters().keySet()) {

            for (Cluster clust : cutoff.getClusters().get(sample)) {
                int clusterID = clust.getId();
                maxClustID = clusterID;
                allOTUset.add(clusterID);

                HashMap<Integer, Integer> countMap = clusterMap.get(sample);

                if (countMap == null) {
                    countMap = new HashMap<Integer, Integer>();
                    clusterMap.put(sample, countMap);
                }
                countMap.put(clusterID, clust.getNumberOfSeqs());
            }
        }

        String otuFormat = "\tOTU_%0" + Integer.toString(maxClustID).length() + "d";
        for (Integer clusterNo : allOTUset) {
            stream.print(String.format(otuFormat, clusterNo));
        }
        stream.println(" ");

        for (String sample : clusterMap.keySet()) {
            // remove "aligned_" or "_trimmed" if present in sample name
            stream.print(sample.replace("aligned_", "").replace("_trimmed", ""));
            HashMap<Integer, Integer> countMap = clusterMap.get(sample);
            for (Integer clusterNo : allOTUset) {

                if (countMap.get(clusterNo) != null) {
                    stream.print("\t" + countMap.get(clusterNo));
                } else {
                    stream.print("\t" + 0);
                }
            }
            stream.println();
        }

        stream.close();
        return true;
    }

    public static boolean createTabulatedFormatForRange(File clusterFile, double distCutoffStart, double distCutoffEnd, File userTempDir) throws IOException {
        boolean distanceFound = false;
        if (distCutoffStart < 0.0 || distCutoffStart > 0.5 || distCutoffStart > distCutoffEnd) {
            throw new IllegalArgumentException("R Format Error: invalid distance cutoff start value");
        }

        if (distCutoffEnd < 0.0 || distCutoffEnd > 0.5 || distCutoffStart > distCutoffEnd) {
            throw new IllegalArgumentException("R Format Error: invalid distance cutoff end value");
        }

        RDPClustParser parser = new RDPClustParser(clusterFile);


        double curDistCutoff = distCutoffStart;
        ArrayList<Double> nonAvailDist = new ArrayList<Double>();
        while (curDistCutoff <= distCutoffEnd) {
            File distfile = new File(userTempDir, "rformat_dist_" + curDistCutoff + ".txt");
            PrintStream writer = new PrintStream(distfile);

            RDPClustParser.Cutoff cutoff = parser.getCutoff(Double.toString(curDistCutoff));
            // you want to update the distanceFound variable to true if there is any distance found, so that the caller of this method would kwow at least a distance was found in the given cutoff range.

            if (cutoff != null) {
                distanceFound = true;
                RFormatter.createTabulatedFormat(cutoff, writer);
            } else {
                // if there is no distance found for the current distance then delete the current distance file (since it's empty now) that was created before the call SpadeInputFormatter.createSimpleFormat(clusterFile, distCutoff, writer).
                distfile.delete();
                // store the distance that was not found in the cluster input file, for writing that to a file after the while loop
                nonAvailDist.add(curDistCutoff);
            }
            curDistCutoff = curDistCutoff + 0.01;
            curDistCutoff = Double.parseDouble(DECIMAL_FORMAT.format(curDistCutoff)); // this is to make sure the double value has 2 digit after deciaml point
        }
        // write the non available distances to a file for the user to let them know.
        if (!nonAvailDist.isEmpty()) {
            File nonAvailDistFile = new File(userTempDir, "rformat_non_available_dist.txt");
            PrintStream nonAvailDistWriter = new PrintStream(nonAvailDistFile);
            nonAvailDistWriter.println("You requested to generate R input format with distance cutoff start value = " + distCutoffStart + ", and the end value = " + distCutoffEnd + ".");
            nonAvailDistWriter.println("Warning: The following distances are not found in your cluster file '" + clusterFile.getName() + "'.");
            for (Double dist : nonAvailDist) {
                nonAvailDistWriter.println(dist);
            }
            nonAvailDistWriter.flush();
            nonAvailDistWriter.close();
        }
        parser.close();
        return distanceFound;
    }

    public static void main(String[] args) throws IOException {
        if (args.length != 4) {
            throw new IllegalArgumentException("Usage: clusterFile outdir startDist endDist");
        }

        File clusterFile = new File(args[0]);
        File userTempDir = new File(args[1]);
        Double startDist = Double.parseDouble(args[2]);
        Double endDist = Double.parseDouble(args[3]);


        RFormatter.createTabulatedFormatForRange(clusterFile, startDist, endDist, userTempDir);


    }
}

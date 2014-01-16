package edu.msu.cme.rdp.abundstats;

import edu.msu.cme.pyro.cluster.io.RDPClustParser;
import edu.msu.cme.pyro.cluster.io.RDPClustParser.ClusterSample;
import edu.msu.cme.pyro.cluster.io.RDPClustParser.Cutoff;
import edu.msu.cme.pyro.cluster.utils.Cluster;
import java.io.File;
import java.io.PrintWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.Set;
import java.util.Map;
import java.util.List;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;

public class ShannonChao {

    private static final DecimalFormat dFormat = new DecimalFormat();

    static {
        dFormat.setMinimumFractionDigits(3);
        dFormat.setMaximumFractionDigits(5);
    }

    public static class ShannonResult {

        private double h;
        private double varH;
        private double evenness;

        public ShannonResult(double h, double varH, double evenness) {
            this.h = h;
            this.varH = varH;
            this.evenness = evenness;
        }

        public double getEvenness() {
            return evenness;
        }

        public double getH() {
            return h;
        }

        public double getVarH() {
            return varH;
        }
    }

    public static class ChaoResult {

        private double chao;
        private double LCI95;
        private double UCI95;

        public ChaoResult(double chao, double LCI95, double UCI95) {
            this.chao = chao;
            this.LCI95 = LCI95;
            this.UCI95 = UCI95;
        }

        public double getLCI95() {
            return LCI95;
        }

        public double getUCI95() {
            return UCI95;
        }

        public double getChao() {
            return chao;
        }
    }

    public static class ShannonChaoResult {

        private double distance;
        private int numClusters;
        private ChaoResult chaoResult;
        private ShannonResult shannonResult;

        public ShannonChaoResult(double distance, int numClusters, ChaoResult chaoResult, ShannonResult shannonResult) {
            this.distance = distance;
            this.chaoResult = chaoResult;
            this.shannonResult = shannonResult;
            this.numClusters = numClusters;
        }

        public ChaoResult getChaoResult() {
            return chaoResult;
        }

        public double getDistance() {
            return distance;
        }

        public ShannonResult getShannonResult() {
            return shannonResult;
        }

        public int getNumClusters() {
            return numClusters;
        }
    }

    /**
     *
     *
     * @param parser
     * @return Map of sample names to ShannonChaoResults, one for each cutoff in
     * the cluster file
     */
    public static Map<ClusterSample, List<ShannonChaoResult>> process(RDPClustParser parser) throws IOException {
        Map<ClusterSample, List<ShannonChaoResult>> ret = new LinkedHashMap();

        for (ClusterSample sample : parser.getClusterSamples()) {
            ret.put(sample, new ArrayList());
        }

        Cutoff cutoff;
        while ((cutoff = parser.readNextCutoff()) != null) {
            for (ClusterSample sample : parser.getClusterSamples()) {
                List<Cluster> sampleClusters = cutoff.getClusters().get(sample.getName());

                Map<Integer, Integer> clusterSizeToCounts = new LinkedHashMap();
                clusterSizeToCounts.put(1, 0);
                clusterSizeToCounts.put(2, 0);

                for (Iterator<Cluster> iterator = sampleClusters.listIterator(); iterator.hasNext();) {
                    Cluster c = iterator.next();
                    if (c.getNumberOfSeqs() == 0) {
                        iterator.remove();
                    }
                }

                for (Cluster c : sampleClusters) {
                    if (!clusterSizeToCounts.containsKey(c.getNumberOfSeqs())) {
                        clusterSizeToCounts.put(c.getNumberOfSeqs(), 0);
                    }

                    clusterSizeToCounts.put(c.getNumberOfSeqs(), clusterSizeToCounts.get(c.getNumberOfSeqs()) + 1);
                }

                ChaoResult chaoResult = chaoI(sample.getSeqs(), sampleClusters.size(), clusterSizeToCounts.get(1), clusterSizeToCounts.get(2));
                ShannonResult shannonResult = shannon(clusterSizeToCounts, sample.getSeqs(), sampleClusters.size());

                ret.get(sample).add(new ShannonChaoResult(Double.parseDouble(cutoff.getCutoff()), sampleClusters.size(), chaoResult, shannonResult));
            }
        }

        return ret;
    }

    private static ChaoResult chaoI(double totalSequences, double totalClusters, double singles, double doubles) {
        double chao = 0;
        double var = Double.NaN;
        if (singles > 0 || (singles == 0 && doubles == 0)) {
            chao = totalClusters + (singles * (singles - 1) / (2 * (doubles + 1)));
        } else {
            chao = totalClusters + (singles * singles) / (2 * doubles);
        }

        if (singles > 0 && doubles > 0) {
            double v1 = singles * (singles - 1) / (2 * (doubles + 1));
            double v2 = singles * Math.pow((2 * singles - 1), 2) / (4 * Math.pow((doubles + 1), 2));
            double v3 = singles * singles * doubles * Math.pow((singles - 1), 2) / (4 * Math.pow((doubles + 1), 4));
            var = v1 + v2 + v3;

        } else if (singles > 0 && doubles == 0) {

            double v1 = singles * (singles - 1) / 2;
            double v2 = singles * Math.pow((2 * singles - 1), 2) / 4;
            double v3 = Math.pow(singles, 4) / (4 * chao);
            var = v1 + v2 - v3;

        } else if (singles == 0 && doubles > 0) {
            var = totalClusters * Math.exp(-totalSequences / totalClusters) * (1 - Math.exp(-totalSequences / totalClusters));
        }

        double c = Math.exp(1.96 * Math.sqrt(Math.log(1 + (var / Math.pow(chao - totalClusters, 2)))));

        double LCI95 = totalClusters + (chao - totalClusters) / c;
        double UCI95 = totalClusters + (chao - totalClusters) * c;

        return new ChaoResult(chao, LCI95, UCI95);
    }

    private static ShannonResult shannon(Map<Integer, Integer> clusterSizeToCounts, double totalSequences, double totalClusters) {
        double H = 0;
        double accV = 0;

        for (Integer cSize : clusterSizeToCounts.keySet()) {
            double clusterCount = clusterSizeToCounts.get(cSize);
            double clusterSize = cSize;


            H -= clusterCount * (clusterSize / totalSequences * Math.log(clusterSize / totalSequences));
            accV += clusterCount * (clusterSize / totalSequences * Math.pow(Math.log(clusterSize / totalSequences), 2));
        }

        double varH = (accV - Math.pow(H, 2)) / totalSequences + (totalClusters - 1) / (2 * Math.pow(totalSequences, 2));
        double Evenness = H / Math.log(totalClusters);

        return new ShannonResult(H, varH, Evenness);
    }

    public static void computeShannonChao1(List<File> clusterFiles, PrintStream out) throws IOException {
        RDPClustParser parser;
        Map<ClusterSample, List<ShannonChaoResult>> results;

        out.println("sampleID\tdistance\tN\tclusters\tchao\tLCI95\tUCI95\tH'\tvarH\tE");
        Set<String> seenSampleNames = new LinkedHashSet();
        for (File clusterFile : clusterFiles) {

            parser = new RDPClustParser(clusterFile);
            results = ShannonChao.process(parser);
            parser.close();

            for (ClusterSample sample : results.keySet()) {
                if (seenSampleNames.contains(sample.getName())) {
                    throw new IOException("Sample name " + sample.getName() + " appears more than once");
                }

                seenSampleNames.add(sample.getName());

                for (ShannonChaoResult result : results.get(sample)) {
                    out.println(sample.getName() + "\t" + result.getDistance() + "\t" + sample.getSeqs() + "\t" + result.getNumClusters() + "\t"
                            + dFormat.format(result.getChaoResult().getChao()) + "\t"
                            + dFormat.format(result.getChaoResult().getLCI95()) + "\t"
                            + dFormat.format(result.getChaoResult().getUCI95()) + "\t"
                            + dFormat.format(result.getShannonResult().getH()) + "\t"
                            + dFormat.format(result.getShannonResult().getVarH()) + "\t"
                            + dFormat.format(result.getShannonResult().getEvenness()));
                }

            }
        }

    }

    /**
     *
     * @param args arg1 = directory that contains one or more cluster file arg2
     * = output file
     * @throws IOException
     */
    public static void main(String[] args) throws IOException {
        if (args.length != 2) {
            System.err.println("USAGE: Shannon <input dir> <output_file>");
            System.err.println("Where <input dir> contains one or more cluster files to process (cluster files must have .clust extension");
            return;
        }

        File dir = new File(args[0]);
        File outfile = new File(args[1]);
        File[] childFiles;
        if (dir.isDirectory()) {
            childFiles = dir.listFiles();
        } else if (dir.exists()) {
            childFiles = new File[]{dir};
        } else {
            throw new IOException(dir.getAbsolutePath() + " doesn't exist");
        }

        PrintWriter writer = null;
        try {
            writer = new PrintWriter(outfile);
            writer.println("sampleID\tdistance\tN\tclusters\tchao\tLCI95\tUCI95\tH'\tvarH\tE");
            Set<String> seenSampleNames = new LinkedHashSet();

            for (int i = 0; i < childFiles.length; i++) {
                if (!childFiles[i].isDirectory() && childFiles[i].getName().endsWith(".clust")) {
                    RDPClustParser parser = new RDPClustParser(childFiles[i]);
                    Map<ClusterSample, List<ShannonChaoResult>> results = ShannonChao.process(parser);
                    parser.close();

                    for (ClusterSample sample : results.keySet()) {
                        if (seenSampleNames.contains(sample.getName())) {
                            throw new IOException("Sample name " + sample.getName() + " appears more than once");
                        }

                        seenSampleNames.add(sample.getName());

                        for (ShannonChaoResult result : results.get(sample)) {
                            writer.println(sample.getName() + "\t" + result.getDistance() + "\t" + sample.getSeqs() + "\t" + result.getNumClusters() + "\t"
                                    + dFormat.format(result.getChaoResult().getChao()) + "\t"
                                    + dFormat.format(result.getChaoResult().getLCI95()) + "\t"
                                    + dFormat.format(result.getChaoResult().getUCI95()) + "\t"
                                    + dFormat.format(result.getShannonResult().getH()) + "\t"
                                    + dFormat.format(result.getShannonResult().getVarH()) + "\t"
                                    + dFormat.format(result.getShannonResult().getEvenness()));
                        }

                    }
                }
            }
        } finally {
            if (writer != null) {
                writer.flush();
                writer.close();
            }
        }
    }
}

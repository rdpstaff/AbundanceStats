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
package edu.msu.cme.rdp.abundstats.cli;

import edu.msu.cme.pyro.cluster.io.RDPClustParser;
import edu.msu.cme.pyro.cluster.io.RDPClustParser.ClusterSample;
import edu.msu.cme.pyro.cluster.io.RDPClustParser.Cutoff;
import edu.msu.cme.pyro.cluster.utils.Cluster;
import edu.msu.cme.rdp.abundstats.AbundStatUtils;
import edu.msu.cme.rdp.abundstats.AbundStatsCalculator;
import edu.msu.cme.rdp.abundstats.RPlotter;
import edu.msu.cme.rdp.abundstats.Sample;
import edu.msu.cme.rdp.abundstats.calculators.Jaccard;
import edu.msu.cme.rdp.abundstats.calculators.Sorensen;
import java.io.*;
import java.util.List;
import java.util.ArrayList;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.PosixParser;

/**
 *
 * @author fishjord
 */
public class AbundMain {

    private static final String rplotterTemplate = "png(filename=\"%{heatmap_file_name}\", width=1200, height=1200)\n"
            + "col_grad <- c(\"#0000FF\",\"#0202FF\",\"#0505FF\",\"#0707FF\",\"#0A0AFF\",\"#0C0CFF\",\"#0F0FFF\",\"#1111FF\",\"#1414FF\",\"#1616FF\",\n"
            + "\"#1919FF\",\"#1C1CFF\",\"#1E1EFF\",\"#2121FF\",\"#2323FF\",\"#2626FF\",\"#2828FF\",\"#2B2BFF\",\"#2D2DFF\",\"#3030FF\",\"#3232FF\",\n"
            + "\"#3535FF\",\"#3838FF\",\"#3A3AFF\",\"#3D3DFF\",\"#3F3FFF\",\"#4242FF\",\"#4444FF\",\"#4747FF\",\"#4949FF\",\"#4C4CFF\",\"#4F4FFF\",\n"
            + "\"#5151FF\",\"#5454FF\",\"#5656FF\",\"#5959FF\",\"#5B5BFF\",\"#5E5EFF\",\"#6060FF\",\"#6363FF\",\"#6565FF\",\"#6868FF\",\"#6B6BFF\",\n"
            + "\"#6D6DFF\",\"#7070FF\",\"#7272FF\",\"#7575FF\",\"#7777FF\",\"#7A7AFF\",\"#7C7CFF\",\"#7F7FFF\",\"#8282FF\",\"#8484FF\",\"#8787FF\",\n"
            + "\"#8989FF\",\"#8C8CFF\",\"#8E8EFF\",\"#9191FF\",\"#9393FF\",\"#9696FF\",\"#9999FF\",\"#9B9BFF\",\"#9E9EFF\",\"#A0A0FF\",\"#A3A3FF\",\n"
            + "\"#A5A5FF\",\"#A8A8FF\",\"#AAAAFF\",\"#ADADFF\",\"#AFAFFF\",\"#B2B2FF\",\"#B5B5FF\",\"#B7B7FF\",\"#BABAFF\",\"#BCBCFF\",\"#BFBFFF\",\n"
            + "\"#C1C1FF\",\"#C4C4FF\",\"#C6C6FF\",\"#C9C9FF\",\"#CCCCFF\",\"#CECEFF\",\"#D1D1FF\",\"#D3D3FF\",\"#D6D6FF\",\"#D8D8FF\",\"#DBDBFF\",\n"
            + "\"#DDDDFF\",\"#E0E0FF\",\"#E2E2FF\",\"#E5E5FF\",\"#E8E8FF\",\"#EAEAFF\",\"#EDEDFF\",\"#EFEFFF\",\"#F2F2FF\",\"#F4F4FF\",\"#F7F7FF\",\"#F9F9FF\",\"#FFFFFF\")\n"
            + "distFile <- read.table(\"%{plot_input_matrix}\", header= TRUE, sep=\"\", dec=\".\", fill = TRUE, row.names = 1)\n"
            + "matrix <- as.matrix(distFile)\n"
            + "distMatrix <- as.dist(matrix)\n"
            + "hc <- hclust(distMatrix, \"average\")\n"
            + "dend <- as.dendrogram(hc)\n"
            + "hv <- heatmap(matrix, Rowv = dend, Colv = \"Rowv\", col=col_grad, scale=\"none\",\n"
            + "        symm = TRUE, margin=c(7,7),\n"
            + "        xlab = \"xseq\", ylab= \"yseq\",\n"
            + "        main = \"%{plot_title}\")\n"
            + "png(filename=\"%{plot_file_name}\", width=1200, height=1024)\n"
            + "plot(hc, hang=-1, cex=1, xlab = \"\", ylab=\"distance\", main = \"%{plot_title}\", sub=\"UPGMA Clustering\")";
    private static final Options options = new Options();

    static {
        options.addOption("r", "result-dir", true, "Directory to put the result files in (default=.)");
        options.addOption("t", "otu-table", false, "input file is an otu table, not rdp cluster file");
        options.addOption("R", "R-location", true, "Triggers the R plotter subsystem, provide the path to the R command");
        options.addOption("l", "lower-cutoff", true, "Lowest cutoff in the cluster file to compute stats for");
        options.addOption("u", "upper-cutoff", true, "Highest cutoff in the cluster file to compute stats for");
        options.addOption("j", "jaccard", false, "Compute jaccard abundance stats");
        options.addOption("s", "sorensen", false, "Compute sorensen abundance stats");
    }

    private static void processSamples(List<Sample> samples, List<AbundStatsCalculator> statCalcs, File resultDir, String prefix, RPlotter plotter) throws IOException {
        for (AbundStatsCalculator calc : statCalcs) {
            float[][] matrix = AbundStatUtils.computeSymMatrix(samples, calc);
            File matrixOut = new File(resultDir, prefix + calc.getCalcName() + ".txt");
            PrintStream out = new PrintStream(matrixOut);
            AbundStatUtils.writeMatrix(out, samples, matrix);
            out.close();

            if (plotter != null) {
                plotter.generateRPlot(resultDir.getParentFile(), matrixOut, new File(resultDir, prefix + calc.getCalcName()), prefix + calc.getCalcName());
            }
        }

    }

    public static void main(String[] args) throws IOException {
        File inputFile;
        File resultDir = new File(".");
        RPlotter plotter = null;
        boolean isClusterFile = true;
        List<AbundStatsCalculator> statCalcs = new ArrayList();
        double clustCutoffFrom = Double.MIN_VALUE, clustCutoffTo = Double.MAX_VALUE;

        String usage = "Main [options] <cluster file>";
        try {
            CommandLine line = new PosixParser().parse(options, args);

            if (line.hasOption("result-dir")) {
                resultDir = new File(line.getOptionValue("result-dir"));
                if (!resultDir.exists() && !resultDir.mkdirs()) {
                    throw new Exception("Result directory " + resultDir + " does not exist and could not be created");
                }
            }

            if (line.hasOption("R-location")) {
                plotter = new RPlotter();
                plotter.setCommandTemplate(rplotterTemplate);
                plotter.setRPath(line.getOptionValue("R-location"));
                plotter.setOutFileExt(".png");

                if (!new File(plotter.getRPath()).canExecute()) {
                    throw new Exception(plotter.getRPath() + " does not exist or is not exectuable");
                }
            }

            if (line.hasOption("lower-cutoff")) {
                clustCutoffFrom = Double.valueOf(line.getOptionValue("lower-cutoff"));
            }

            if (line.hasOption("upper-cutoff")) {
                clustCutoffTo = Double.valueOf(line.getOptionValue("upper-cutoff"));
            }

            if (line.hasOption("jaccard")) {
                statCalcs.add(new Jaccard(true));
            }

            if (line.hasOption("sorensen")) {
                statCalcs.add(new Sorensen(true));
            }

            if (line.hasOption("otu-table")) {
                isClusterFile = false;
            }

            if (statCalcs.isEmpty()) {
                throw new Exception("Must specify at least one stat to compute (jaccard, sorensen)");
            }

            args = line.getArgs();
            if (args.length != 1) {
                throw new Exception("Unexpected number of command line arguments");
            }

            inputFile = new File(args[0]);

        } catch (Exception e) {
            new HelpFormatter().printHelp(usage, options);
            System.err.println("Error: " + e.getMessage());
            return;
        }

        if (isClusterFile) {
            RDPClustParser parser;
            parser = new RDPClustParser(inputFile);

            try {
                if (parser.getClusterSamples().size() == 1) {
                    throw new IOException("Cluster file must have more than one sample");
                }

                List<Cutoff> cutoffs = parser.getCutoffs(clustCutoffFrom, clustCutoffTo);
                if (cutoffs.isEmpty()) {
                    throw new IOException("No cutoffs in cluster file in range [" + clustCutoffFrom + "-" + clustCutoffTo + "]");
                }

                for (Cutoff cutoff : cutoffs) {
                    List<Sample> samples = new ArrayList();

                    for (ClusterSample clustSample : parser.getClusterSamples()) {
                        Sample s = new Sample(clustSample.getName());
                        for (Cluster clust : cutoff.getClusters().get(clustSample.getName())) {
                            s.addSpecies(clust.getNumberOfSeqs());
                        }
                        samples.add(s);
                    }

                    processSamples(samples, statCalcs, resultDir, cutoff.getCutoff() + "_", plotter);
                }


            } finally {
                parser.close();
            }
        } else {
            List<Sample> samples = new ArrayList();
            BufferedReader reader = new BufferedReader(new FileReader(inputFile));
            String line = reader.readLine();

            if (line == null || line.split("\\s+").length < 2) {
                throw new IOException("Must be 2 or more samples for abundance statistic calculations!");
            }
            int numSamples = line.split("\\s+").length;

            boolean header = true;
            try {
                Integer.valueOf(line.split("\\s+")[0]);
                header = false;
            } catch (Exception e) {
            }

            if (header) {
                for (String s : line.split("\\s+")) {
                    samples.add(new Sample(s));
                }
            } else {
                int sample = 0;
                for (String s : line.split("\\s+")) {
                    samples.add(new Sample("" + sample));
                    samples.get(sample).addSpecies(Integer.valueOf(s));
                    sample++;
                }
            }

            int lineno = 2;
            while ((line = reader.readLine()) != null) {
                if (line.trim().equals("")) {
                    continue;
                }
                int sample = 0;
                if (line.split("\\s+").length != numSamples) {
                    System.err.println("Line number " + lineno + " didn't have the expected number of samples (contained " + line.split("\\s+").length + ", expected " + numSamples + ")");
                }

                for (String s : line.split("\\s+")) {
                    samples.get(sample).addSpecies(Integer.valueOf(s));
                    sample++;
                }

                lineno++;
            }

            processSamples(samples, statCalcs, resultDir, inputFile.getName(), plotter);
        }
    }
}

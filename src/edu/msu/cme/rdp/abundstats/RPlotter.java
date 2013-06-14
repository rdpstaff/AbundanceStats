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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.io.Reader;

/**
 *
 * @author fishjord
 */
public class RPlotter {
    private static final String PLOT_FILE_NAME = "%{plot_file_name}";
    private static final String PLOT_TITLE = "%{plot_title}";
    private static final String PLOT_INPUT = "%{plot_input_matrix}";
    private static final String HEATMAP_FILE_NAME = "%{heatmap_file_name}";

    private String RPath;
    private String commandTemplate;
    private String outFileExt;

    public String getRPath() {
        return RPath;
    }

    public void setRPath(String RPath) {
        this.RPath = RPath;
    }

    public String getCommandTemplate() {
        return commandTemplate;
    }

    public void setCommandTemplate(String commandTemplate) {
        if(commandTemplate.indexOf(PLOT_FILE_NAME) == -1 || commandTemplate.indexOf(PLOT_FILE_NAME) != commandTemplate.lastIndexOf(PLOT_FILE_NAME)) {
            throw new IllegalArgumentException("Command template must contain exactly one occurance of " + PLOT_FILE_NAME + " for inserting the output file name");
        }

        if(commandTemplate.indexOf(PLOT_TITLE) == -1) {
            throw new IllegalArgumentException("Command template must contain exactly one occurance of " + PLOT_TITLE + " for inserting the plot title");
        }

        if(commandTemplate.indexOf(PLOT_INPUT) == -1 || commandTemplate.indexOf(PLOT_INPUT) != commandTemplate.lastIndexOf(PLOT_INPUT)) {
            throw new IllegalArgumentException("Command template must contain exactly one occurance of " + PLOT_INPUT + " for inserting the input matrix file");
        }

        if(commandTemplate.indexOf(HEATMAP_FILE_NAME) == -1 || commandTemplate.indexOf(HEATMAP_FILE_NAME) != commandTemplate.lastIndexOf(HEATMAP_FILE_NAME)) {
            throw new IllegalArgumentException("Command template must contain exactly one occurance of " + HEATMAP_FILE_NAME + " for inserting the output heatmap file");
        }

        this.commandTemplate = commandTemplate;
    }

    public String getOutFileExt() {
        return outFileExt;
    }

    public void setOutFileExt(String outFileExt) {
        this.outFileExt = outFileExt;
    }

    private static void copyToStdErr(Reader r) throws IOException {
        BufferedReader reader = new BufferedReader(r);
        String line;

        while((line = reader.readLine()) != null) {
            System.err.println(line);
        }
        reader.close();
    }

    public void generateRPlot(File workingDir, File inputFile, File outputFile, String plotTitle) throws IOException {
        File tmpScriptFile = new File(workingDir, "dendo_script_tmp.R");
        PrintStream out = new PrintStream(tmpScriptFile);
        out.println(commandTemplate
                .replace(PLOT_FILE_NAME, outputFile.getAbsolutePath() + "_tree" + outFileExt)
                .replace(PLOT_TITLE, plotTitle)
                .replace(PLOT_INPUT, inputFile.getAbsolutePath())
                .replace(HEATMAP_FILE_NAME, outputFile.getAbsolutePath() + outFileExt));
        out.close();

        ProcessBuilder builder = new ProcessBuilder(RPath, "CMD", "BATCH", tmpScriptFile.getAbsolutePath());
        builder.directory(workingDir);
        Process p = builder.start();

        try {
            if(p.waitFor() != 0) {
                System.err.println("RPlotter failed to run");
                try {
                    copyToStdErr(new InputStreamReader(p.getErrorStream()));
                    copyToStdErr(new FileReader(tmpScriptFile.getAbsolutePath() + "out"));
                } catch(Exception ignore) {}
                System.err.println();
            }
        } catch(InterruptedException e) {
            throw new RuntimeException(e);
        }
    }

    public static void main(String [] args) throws Exception {
        RPlotter plotter = new RPlotter();
        plotter.setRPath("/usr/bin/R");
        plotter.setOutFileExt(".png");
        plotter.setCommandTemplate("png(filename=\"%{plot_file_name}\", width=1024, height=768)\n" +
                "distFile <- read.table(\"%{plot_input_matrix}\", header= TRUE, sep=\"\", dec=\".\", fill = TRUE, row.names = 1)\n" +
                "hc <- hclust(as.dist(as.matrix(distFile)), \"average\")\n" +
                "plot(hc, hang=-1, cex=1, xlab = \"\", ylab=\"distance\", main = \"%{plot_title}\", sub=\"UPGMA Clustering\")\n");

        File workingDir = new File("/home/fishjord/tmp");
        File inFile = new File("/home/fishjord/tmp/jaccard.txt");
        File out = new File("/home/fishjord/tmp/out");

        plotter.generateRPlot(workingDir, inFile, out, "Test plot");
    }
}

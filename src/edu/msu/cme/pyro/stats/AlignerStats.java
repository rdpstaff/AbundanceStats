/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.msu.cme.pyro.stats;

import edu.msu.cme.pyro.derep.IdMapping;
import edu.msu.cme.pyro.derep.SampleMapping;
import edu.msu.cme.rdp.readseq.readers.SeqReader;
import edu.msu.cme.rdp.readseq.readers.SequenceReader;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.readseq.utils.ReferenceSeqHelper;
import edu.msu.cme.rdp.readseq.utils.SequenceStats;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.geom.RectangularShape;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.renderer.xy.XYBarPainter;
import org.jfree.chart.renderer.xy.XYBarRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RectangleEdge;

/**
 *
 * @author fishjord
 */
public class AlignerStats {

    private static class SeqEndPoints {

        private String seqid;
        private int start;
        private int end;

        public SeqEndPoints(String seqid, SequenceStats stats) {
            this.seqid = seqid;
            this.start = stats.getStart();
            this.end = stats.getEnd();
        }

        public SeqEndPoints(String seqid, SequenceStats stats, String model) {
            this.seqid = seqid;
            start = ReferenceSeqHelper.getInstance().translate(model, stats.getStart());
            end = ReferenceSeqHelper.getInstance().translate(model, stats.getEnd());
        }

        public int getStart() {
            return start;
        }

        public int getEnd() {
            return end;
        }

        public String getSeqid() {
            return seqid;
        }
    }

    private static void inc(Map<Integer, Integer> map, int key) {
        int v = 0;
        if (map.containsKey(key)) {
            v = map.get(key);
        }

        map.put(key, v + 1);
    }
    private List<SeqEndPoints> seqEndPoints = new ArrayList();
    private List<Double> bitsSaved = new ArrayList();
    private int modelPositionSigma = 0;
    private int lengthSigma = 0;
    private int allGapSeqs = 0;
    private String model;

    public AlignerStats() {
        this.model = null;
    }

    public AlignerStats(String model) {
        this.model = model;
    }

    public void readStats(File seqFile, File scoresFile) throws IOException {
        String line;
        String[] lexemes;
        Map<String, Double> bitsMap = new HashMap();
        
        BufferedReader reader = new BufferedReader(new FileReader(scoresFile));
        try {
            while((line = reader.readLine()) != null) {
                lexemes = line.split("\\s+");
                if(lexemes.length != 12) {
                    continue;
                }

                bitsMap.put(lexemes[1], Double.valueOf(lexemes[6]));
            }
        } finally {
            reader.close();
        }
        
        Sequence seq;
        SeqReader seqReader = new SequenceReader(seqFile);
        try {
            while((seq = seqReader.readNextSequence()) != null) {
                Double bits = bitsMap.get(seq.getSeqName());
                if(bits == null) {
                    bits = -1.0d;
                }
                
                recordSequence(seq, bits);
            }
        } finally {
            seqReader.close();
        }
    }
    
    public void recordSequence(Sequence seq, double bitsSaved) {
        recordSequence(seq.getSeqName(), seq.getSeqString(), bitsSaved);
    }

    public void recordSequence(String seqid, String seqString, double bitsSaved) {
        SequenceStats seqStats = new SequenceStats(seqString);

        if (seqStats.getLength() == 0) {
            allGapSeqs++;
        } else {
            modelPositionSigma += seqStats.getModelPositions();
            lengthSigma += seqStats.getLength();

            if (model == null) {
                seqEndPoints.add(new SeqEndPoints(seqid, seqStats));
            } else {
                seqEndPoints.add(new SeqEndPoints(seqid, seqStats, model));
            }
            this.bitsSaved.add(bitsSaved);
        }
    }

    public void writeStats(File statsFile, File statsImg) throws IOException {
        Map<Integer, Integer> startPosMap = new LinkedHashMap();
        Map<Integer, Integer> endPosMap = new LinkedHashMap();

        int numSeqs = seqEndPoints.size();
        Collections.sort(bitsSaved);

        if (numSeqs == 0) {
            return;
        }

        String modelRefSeq = "(none)";
        if (model != null) {
            modelRefSeq = ReferenceSeqHelper.getInstance().getRefSeq(model);
        }

        PrintStream statsOut = new PrintStream(statsFile);
        statsOut.println("Total Sequences\t" + numSeqs);
        statsOut.println("Alignment Model\t" + ((model == null) ? "Unspecified" : model));
        statsOut.println("Alignment Reference Sequence\t" + modelRefSeq);
        statsOut.println("All gap sequences\t" + allGapSeqs);
        statsOut.println("Average number of comparable positions\t" + modelPositionSigma / numSeqs);
        statsOut.println("Average sequence length (model position only)\t" + lengthSigma / numSeqs);
        /*statsOut.println("Minimum bits saved\t" + bitsSaved.get(0));
        statsOut.println("Maximum bits saved\t" + bitsSaved.get(bitsSaved.size() - 1));
        statsOut.println("Average bits saved\t" + StatsHelper.getAverage(bitsSaved));*/
        statsOut.println();

        statsOut.println("Sequence Id\tStart position\tEnd position");
        for (SeqEndPoints endPoints : seqEndPoints) {
            statsOut.println(endPoints.getSeqid() + "\t" + endPoints.getStart() + "\t" + endPoints.getEnd());
            inc(startPosMap, endPoints.getStart());
            inc(endPosMap, endPoints.getEnd());
        }

        statsOut.close();


        XYSeriesCollection startEndPosColl = new XYSeriesCollection();
        XYSeries startPosSeries = new XYSeries("Start positions");
        XYSeries endPosSeries = new XYSeries("End positions");

        for (int key : startPosMap.keySet()) {
            startPosSeries.add(key, startPosMap.get(key));
        }
        startEndPosColl.addSeries(startPosSeries);

        for (int key : endPosMap.keySet()) {
            endPosSeries.add(key, endPosMap.get(key));
        }
        startEndPosColl.addSeries(endPosSeries);

        JFreeChart chart = ChartFactory.createXYBarChart("Aligned Sequence Start & End Position Histogram", modelRefSeq + " Position", false, "Number of sequences", startEndPosColl, PlotOrientation.VERTICAL, true, true, false);

        XYBarRenderer renderer = (XYBarRenderer) chart.getXYPlot().getRenderer();
        renderer.setShadowVisible(false);
        renderer.setBarPainter(new XYBarPainter() {

            public void paintBar(Graphics2D arg0, XYBarRenderer arg1, int arg2, int arg3, RectangularShape arg4, RectangleEdge arg5) {
                Rectangle r = arg4.getBounds();
                arg0.setPaint(arg1.getItemPaint(arg2, arg3));
                arg0.fillRect((int) r.getX(), (int) r.getY(), (int) r.getWidth(), (int) r.getHeight());
            }

            public void paintBarShadow(Graphics2D arg0, XYBarRenderer arg1, int arg2, int arg3, RectangularShape arg4, RectangleEdge arg5, boolean arg6) {
                throw new UnsupportedOperationException("Not supported yet.");
            }
        });

        ChartUtilities.writeChartAsPNG(new PrintStream(statsImg), chart, 640, 640);
    }

    public static void main(String... args) throws Exception {
        if (args.length < 2) {
            System.out.println("USAGE: AlignerStats [-m <model>] [-i <id-mapping>] [-s <sample mapping>] <outfile_prefix> <fasta_file>...");
            return;
        }

        int argIndex = 0;
        Map<String, AlignerStats> statCalculators = new LinkedHashMap();
        String maskSeq = null;
        IdMapping<Integer> idMapping = null;
        SampleMapping<String> sampleMapping = null;
        Map<String, Integer> reverseMapping = null;

        if (args[argIndex].equals("-m")) {
            argIndex++;
            maskSeq = args[argIndex++];
        }

        if (args[argIndex].equals("-i")) {
            argIndex++;
            idMapping = IdMapping.fromFile(new File(args[argIndex++]));
            reverseMapping = idMapping.getReverseMapping();
        }

        if (args[argIndex].equals("-s")) {
            argIndex++;
            sampleMapping = SampleMapping.fromFile(new File(args[argIndex++]));

            for (String sample : sampleMapping.getSampleList()) {

                statCalculators.put(sample, new AlignerStats(maskSeq));
            }
        }

        String prefix = args[argIndex++];

        for (int index = argIndex; index < args.length; index++) {
            File seqFile = new File(args[index]);
            if (sampleMapping == null) {
                statCalculators.put(seqFile.getName(), new AlignerStats(maskSeq));
            }

            SequenceReader reader = new SequenceReader(seqFile);
            Sequence seq;
            while ((seq = reader.readNextSequence()) != null) {
                if (seq.getSeqName().startsWith("#")) {
                    continue;


                }
                List<String> seqids;
                if (idMapping != null) {
                    Integer exid = reverseMapping.get(seq.getSeqName());
                    if (exid == null) {
                        System.out.println("No id mapping found for " + seq.getSeqName());
                        seqids = Arrays.asList(seq.getSeqName());
                    } else {
                        seqids = idMapping.getIds(exid);
                    }
                } else {
                    seqids = Arrays.asList(seq.getSeqName());
                }

                for (String seqid : seqids) {
                    if (sampleMapping == null) {
                        statCalculators.get(seqFile.getName()).recordSequence(seqid, seq.getSeqString(), 0);
                    } else {
                        String sample = sampleMapping.getSampleById(seqid);
                        if (sample == null) {
                            System.out.println("No sample mapping found for seqid " + seqid);
                        } else {
                            statCalculators.get(sample).recordSequence(seqid, seq.getSeqString(), 0);
                        }
                    }
                }
            }

            reader.close();
        }

        for (String sample : statCalculators.keySet()) {
            statCalculators.get(sample).writeStats(new File(prefix + "_" + sample + "_stats.txt"), new File(prefix + "_" + sample + "_chart.png"));
        }
    }
}

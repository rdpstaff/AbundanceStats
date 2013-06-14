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
package edu.msu.cme.rdp.seqstats;

import edu.msu.cme.rdp.readseq.readers.IndexedSeqReader;
import edu.msu.cme.rdp.readseq.readers.SequenceReader;
import edu.msu.cme.rdp.readseq.readers.SeqReader;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Arrays;

/**
 *
 * @author fishjord
 */
public class MutualInformation {

    private static final int[] baseMap = new int[127];

    static {
        Arrays.fill(baseMap, -1);
        baseMap['a'] = baseMap['A'] = 0;
        baseMap['c'] = baseMap['C'] = 1;
        baseMap['g'] = baseMap['G'] = 2;
        baseMap['t'] = baseMap['T'] = 3;
        baseMap['-'] = 4;
    }
    
    private static double mutualInformation(float totalSeqs, int x, int y, int[][] posFreqs, int[][][][] pairFreqs) {
        double mxy = 0;
        
        for(int b1 = 0;b1 < 5;b1++) {
            for(int b2 = 0;b2 < 5;b2++) {
                if(posFreqs[x][b1] == 0 || posFreqs[y][b2] == 0 || pairFreqs[x][y][b1][b2] == 0) {
                    continue;
                }
                mxy += (pairFreqs[x][y][b1][b2] / totalSeqs) * Math.log((pairFreqs[x][y][b1][b2] / totalSeqs) / ( (posFreqs[x][b1] / totalSeqs) * (posFreqs[y][b2] / totalSeqs)));
            }
        }
        
        return mxy;
    }

    public static void main(String[] args) throws Exception {
        if (args.length != 2 && args.length != 3) {
            System.err.println("USAGE: MutualInformation <input_sequences> <mi_matrix_out> [mask sequence]");
            System.exit(1);
        }

        SeqReader reader;
        File outFile = new File(args[1]);

        if (args.length == 2) {
            reader = new SequenceReader(new File(args[0]));
        } else {
            reader = new IndexedSeqReader(new File(args[0]), args[2]);
        }

        Sequence seq = reader.readNextSequence();

        if (seq == null) {
            throw new IOException("Sequence file is empty");
        }


        int[][] posFreqs = new int[seq.getSeqString().length()][5];
        int[][][][] pairFreqs = new int[seq.getSeqString().length()][seq.getSeqString().length()][5][5];

        int seqCount = 0;
        do {
            seqCount++;
            if (seq.getSeqName().startsWith("#")) {
                continue;
            }
            char[] bases = seq.getSeqString().toCharArray();
            for (int index = 0; index < bases.length; index++) {
                try {
                    posFreqs[index][baseMap[bases[index]]]++;
                } catch (ArrayIndexOutOfBoundsException e) {
                    throw new RuntimeException("Sequence " + seq.getSeqName() + " bases[" + index + "] = " + bases[index]);
                }

                for (int pairWith = 0; pairWith < bases.length; pairWith++) {
                    if (pairWith == index) {
                        continue;
                    }
                    pairFreqs[index][pairWith][baseMap[bases[index]]][baseMap[bases[pairWith]]]++;
                }
            }
        } while ((seq = reader.readNextSequence()) != null);

        PrintStream out = new PrintStream(outFile);
        for (int x = 0; x < posFreqs.length; x++) {
            out.print(x);
            for (int y = 0; y < posFreqs.length; y++) {
                out.print("\t" + mutualInformation(seqCount, x, y, posFreqs, pairFreqs));
            }
            out.println();
        }
        out.close();
    }
}

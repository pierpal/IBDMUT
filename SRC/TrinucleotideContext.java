/*
Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements.  See the NOTICE file
distributed with this work for additional information
regarding copyright ownership.  The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied.  See the License for the
specific language governing permissions and limitations
under the License.
 */

package IBDMUT;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.TreeMap;
import java.util.zip.GZIPInputStream;

/**
 *
 * @author Pier Palamara
 */
public class TrinucleotideContext {

    public static final String[] trinucleotideListDirectional = {"AAT	--	ACT", "GGA	--	GCA",
        "TTC	--	TAC", "CAA	--	CTA", "TTG	--	TCG", "TAG	--	TGG",
        "TTT	--	TAT", "CCG	--	CAG", "AAG	--	AGG", "ATG	--	ACG",
        "ATC	--	AAC", "CAA	--	CCA", "GAA	--	GTA", "ATT	--	AAT",
        "CCG	--	CTG", "TTC	--	TCC", "TAC	--	TGC", "TCA	--	TAA",
        "GAA	--	GCA", "TTT	--	TCT", "CTA	--	CAA", "TAT	--	TGT",
        "GCG	--	GAG", "ATC	--	ACC", "CCC	--	CAC", "AAC	--	AGC",
        "ACA	--	AAA", "TCA	--	TTA", "CCT	--	CAT", "ATT	--	ACT",
        "AAT	--	AGT", "GCG	--	GTG", "CCC	--	CTC", "TTG	--	TGG",
        "ACA	--	ATA", "CCT	--	CTT", "GTA	--	GAA", "CTA	--	CCA",
        "ATG	--	AGG", "GCC	--	GAC", "CAA	--	CGA", "GCT	--	GAT",
        "GCC	--	GTC", "TTC	--	TGC", "GCT	--	GTT", "GTA	--	GCA",
        "GAA	--	GGA", "TTT	--	TGT", "CCG	--	CGG", "ATC	--	AGC",
        "ATT	--	AGT", "TCA	--	TGA", "CTA	--	CGA", "GCG	--	GGG",
        "CCC	--	CGC", "ACA	--	AGA", "CCT	--	CGT", "GTA	--	GGA",
        "GCC	--	GGC", "GCT	--	GGT", "CGG	--	CAG", "CGG	--	CTG",
        "TGA	--	TAA", "GGG	--	GAG", "CGG	--	CCG", "CGC	--	CAC",
        "GGG	--	GTG", "AGA	--	AAA", "TGA	--	TTA", "CGT	--	CAT",
        "CGC	--	CTC", "CAG	--	CTG", "AGA	--	ATA", "GGG	--	GCG",
        "CGT	--	CTT", "TGA	--	TCA", "GGC	--	GAC", "CGC	--	CCC",
        "GGT	--	GAT", "CAG	--	CCG", "GGC	--	GTC", "AGA	--	ACA",
        "CGT	--	CCT", "TAA	--	TTA", "GGT	--	GTT", "CAC	--	CTC",
        "GAG	--	GTG", "GGC	--	GCC", "AAA	--	ATA", "TAA	--	TCA",
        "CAT	--	CTT", "GGT	--	GCT", "TCG	--	TAG", "CTG	--	CAG",
        "CAC	--	CCC", "GAG	--	GCG", "AAA	--	ACA", "ACG	--	AAG",
        "TCG	--	TTG", "CAT	--	CCT", "GAC	--	GTC", "ACG	--	ATG",
        "GAT	--	GTT", "TTA	--	TAA", "TCC	--	TAC", "CAG	--	CGG",
        "CTC	--	CAC", "CTG	--	CCG", "GTG	--	GAG", "GAC	--	GCC",
        "TCT	--	TAT", "ATA	--	AAA", "GAT	--	GCT", "TCC	--	TTC",
        "ACC	--	AAC", "CTT	--	CAT", "ACT	--	AAT", "TCT	--	TTT",
        "ACC	--	ATC", "TTA	--	TCA", "TAA	--	TGA", "GTG	--	GCG",
        "GTC	--	GAC", "CTC	--	CCC", "GTT	--	GAT", "CAC	--	CGC",
        "ACT	--	ATT", "GAG	--	GGG", "CCA	--	CAA", "AAA	--	AGA",
        "ATA	--	ACA", "CTT	--	CCT", "CAT	--	CGT", "CCA	--	CTA",
        "TCG	--	TGG", "CTG	--	CGG", "GAC	--	GGC", "GTT	--	GCT",
        "GTC	--	GCC", "GCA	--	GAA", "ACG	--	AGG", "GAT	--	GGT",
        "GCA	--	GTA", "TTA	--	TGA", "TCC	--	TGC", "CTC	--	CGC",
        "GTG	--	GGG", "ATA	--	AGA", "TCT	--	TGT", "ACC	--	AGC",
        "CTT	--	CGT", "ACT	--	AGT", "GTC	--	GGC", "GTT	--	GGT",
        "TGG	--	TAG", "CCA	--	CGA", "TGG	--	TTG", "AGG	--	AAG",
        "AGG	--	ATG", "TGG	--	TCG", "TGC	--	TAC", "GCA	--	GGA",
        "TGT	--	TAT", "TGC	--	TTC", "AGC	--	AAC", "AGG	--	ACG",
        "TAG	--	TTG", "TGT	--	TTT", "AGT	--	AAT", "AGC	--	ATC",
        "AAG	--	ATG", "TGC	--	TCC", "TAG	--	TCG", "AGT	--	ATT",
        "CGA	--	CAA", "TGT	--	TCT", "AGC	--	ACC", "AAG	--	ACG",
        "TAC	--	TTC", "CGA	--	CTA", "AGT	--	ACT", "TAT	--	TTT",
        "AAC	--	ATC", "GGA	--	GAA", "TTG	--	TAG", "TAC	--	TCC",
        "CGA	--	CCA", "AAT	--	ATT", "TAT	--	TCT", "GGA	--	GTA",
        "AAC	--	ACC", "ATG	--	AAG"};

    public static final String[] trinucleotideListNoDirection = {"AAA	--	ACA", "AAA	--	AGA", "AAA	--	ATA",
        "AAC	--	ACC", "AAC	--	AGC", "AAC	--	ATC", "AAG	--	ACG",
        "AAG	--	AGG", "AAG	--	ATG", "AAT	--	ACT", "AAT	--	AGT",
        "AAT	--	ATT", "ACA	--	AGA", "ACA	--	ATA", "ACC	--	AGC",
        "ACC	--	ATC", "ACG	--	AGG", "ACG	--	ATG", "ACT	--	AGT",
        "ACT	--	ATT", "AGA	--	ATA", "AGC	--	ATC", "AGG	--	ATG",
        "AGT	--	ATT", "CAA	--	CCA", "CAA	--	CGA", "CAA	--	CTA",
        "CAC	--	CCC", "CAC	--	CGC", "CAC	--	CTC", "CAG	--	CCG",
        "CAG	--	CGG", "CAG	--	CTG", "CAT	--	CCT", "CAT	--	CGT",
        "CAT	--	CTT", "CCA	--	CGA", "CCA	--	CTA", "CCC	--	CGC",
        "CCC	--	CTC", "CCG	--	CGG", "CCG	--	CTG", "CCT	--	CGT",
        "CCT	--	CTT", "CGA	--	CTA", "CGC	--	CTC", "CGG	--	CTG",
        "CGT	--	CTT", "GAA	--	GCA", "GAA	--	GGA", "GAA	--	GTA",
        "GAC	--	GCC", "GAC	--	GGC", "GAC	--	GTC", "GAG	--	GCG",
        "GAG	--	GGG", "GAG	--	GTG", "GAT	--	GCT", "GAT	--	GGT",
        "GAT	--	GTT", "GCA	--	GGA", "GCA	--	GTA", "GCC	--	GGC",
        "GCC	--	GTC", "GCG	--	GGG", "GCG	--	GTG", "GCT	--	GGT",
        "GCT	--	GTT", "GGA	--	GTA", "GGC	--	GTC", "GGG	--	GTG",
        "GGT	--	GTT", "TAA	--	TCA", "TAA	--	TGA", "TAA	--	TTA",
        "TAC	--	TCC", "TAC	--	TGC", "TAC	--	TTC", "TAG	--	TCG",
        "TAG	--	TGG", "TAG	--	TTG", "TAT	--	TCT", "TAT	--	TGT",
        "TAT	--	TTT", "TCA	--	TGA", "TCA	--	TTA", "TCC	--	TGC",
        "TCC	--	TTC", "TCG	--	TGG", "TCG	--	TTG", "TCT	--	TGT",
        "TCT	--	TTT", "TGA	--	TTA", "TGC	--	TTC", "TGG	--	TTG",
        "TGT	--	TTT"};

    public static TreeMap<String, Pair<String, String>> context = new TreeMap<String, Pair<String, String>>();

    private TreeMap<String, Integer> trinucleotideContextCounts = new TreeMap<String, Integer>();

    public TreeMap<String, Integer> getTrinucleotideCounts() {
        return trinucleotideContextCounts;
    }

    public double getTrinucleotideRawCount(String trinucleotide) {
        return (!trinucleotideContextCounts.containsKey(trinucleotide)) ? 0. : trinucleotideContextCounts.get(trinucleotide);
    }

    public double getTrinucleotideProbability(String trinucleotide) {
        double norm = 0.;
        for (String t : trinucleotideContextCounts.keySet()) {
            norm += trinucleotideContextCounts.get(t);
        }
        return (!trinucleotideContextCounts.containsKey(trinucleotide)) ? 0. : trinucleotideContextCounts.get(trinucleotide) / norm;
    }

    public void setTrinucleotideCounts(TreeMap<String, Integer> trinucleotideContextCounts) {
        this.trinucleotideContextCounts = trinucleotideContextCounts;
    }

    public void increaseTrinucleotideCountsBy(String from, String to, int increase) {
        String key = from + "\t--\t" + to;
        int counts = (!trinucleotideContextCounts.containsKey(key)) ? 0 : trinucleotideContextCounts.get(key);
        trinucleotideContextCounts.put(key, counts + increase);
    }

    public static TreeMap<String, Integer> mergeTrinucleotideCounts(TrinucleotideContext tn1, TrinucleotideContext tn2) throws IOException {
        TreeMap<String, Integer> tnc1 = tn1.trinucleotideContextCounts;
        TreeMap<String, Integer> tnc2 = tn2.trinucleotideContextCounts;
        TreeMap<String, Integer> merged = new TreeMap<String, Integer>();
        for (String trinucleotide : tnc1.keySet()) {
            merged.put(trinucleotide, tnc1.get(trinucleotide));
        }
        for (String trinucleotide : tnc2.keySet()) {
            int counts = (!merged.containsKey(trinucleotide)) ? 0 : merged.get(trinucleotide);
            merged.put(trinucleotide, counts + tnc2.get(trinucleotide));
        }
        return merged;
    }

    public void printTrinucleotideResults() {
        Tools.printResult("Trinucleotide results");
        for (String trinucleotide : trinucleotideContextCounts.keySet()) {
            Tools.printResult("\t" + trinucleotide + "\t" + trinucleotideContextCounts.get(trinucleotide));
        }
    }

    public static void setTrinucleotideContext(String filename) throws IOException {
        BufferedReader br = null;
        try {
            InputStream fileStream = new FileInputStream(filename);
            InputStream gzipStream = new GZIPInputStream(fileStream);
            Reader decoder = new InputStreamReader(gzipStream, "UTF-8");
            br = new BufferedReader(decoder);
            String line = br.readLine();
            int cnt = 0;
            while (line != null) {
                String[] split = line.split("\\s+");
                context.put(split[0] + ":" + split[1], new Pair(split[2], split[3]));
                line = br.readLine();
                cnt++;
            }
            Tools.printVerboseProgressLevel1("Read " + cnt + " lines from trinucleotide context file.");
        } catch (Exception ex) {
            Tools.exit("Could not read trinucleotide context file " + filename);
        }
    }

    public static Pair<String, String> getTrinucleotides(String mut) {
        return context.get(mut);
    }

}

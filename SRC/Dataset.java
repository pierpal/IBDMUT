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
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashSet;
import java.util.TreeMap;
import java.util.zip.GZIPInputStream;

/**
 *
 * @author Pier Palamara
 */
public class Dataset {

    private ArrayList<Variant> variants = new ArrayList<Variant>();
    private ArrayList<String> SNP_IDs = new ArrayList<String>();
    private TreeMap<Double, Integer> genToIndex = new TreeMap<Double, Integer>();
    private TreeMap<Integer, Integer> physToIndex = new TreeMap<Integer, Integer>();
    private TreeMap<Double, Integer> genToPhys = new TreeMap<Double, Integer>();
    private TreeMap<String, Double> IDToGen = new TreeMap<String, Double>();
    private TreeMap<String, Double> IDToFreq = new TreeMap<String, Double>();
    private TreeMap<String, Integer> IDToVarCounts = new TreeMap<String, Integer>();
    private ArrayList<Double> genPosArray = new ArrayList<Double>();
    private ArrayList<Integer> physPos = new ArrayList<Integer>();
    private TreeMap<String, Integer> IDtoPos = new TreeMap<String, Integer>();
    private TreeMap<Integer, String> physToID = new TreeMap<Integer, String>();
    private TreeMap<String, Individual> individuals = new TreeMap<String, Individual>();
    private Map map = null;
    private Mask includeMask;
    private Mask excludeMask;
    private boolean haveFreq = false;
    private String chr = "-1";
    private static boolean enterTainUser = false;
    
    public static void setEntertainment() {
        enterTainUser = true;
    }
    
    public static void setRoundTo(double RT) {
        roundTo = RT;
    }

    private static double roundTo = 10.0;

    private double totDiff = 0.0;
    private double totInMask = 0.0;
    private int inMask = 0;
    private int notInMask = 0;
    private boolean warnedPostNotContainsInd = false;

    void readPedFile(String root) throws IOException, Exception {
        String mapFile = root + ".map";
        String pedFile = root + ".ped.gz";
        BufferedReader br = null;
        Tools.printProgress("Reading ped file " + pedFile);
        try {
            InputStream fileStream = new FileInputStream(pedFile);
            InputStream gzipStream = new GZIPInputStream(fileStream);
            Reader decoder = new InputStreamReader(gzipStream, "UTF-8");
            br = new BufferedReader(decoder);
        } catch (FileNotFoundException ex) {
            Tools.exit("Could not open ped file " + pedFile);
        }
        try {
            String line = br.readLine();
            int lineCnt = 0;
            double miss = 0;
            double cnt = 0;
            int length = (line.split("\\s+").length - 6) / 2;
            for (int i = 0; i < length; i++) {
                getVariants().add(new Variant());
            }
            while (line != null) {
                String[] strSplit = line.split("\\s+");
                lineCnt++;
                String famId = strSplit[0];
                String indId = strSplit[1];
                String patId = strSplit[2];
                String matId = strSplit[3];
                String sex = strSplit[4];
                String pheno = strSplit[5];

                String ID = famId + "-" + indId;

                if (patId.compareToIgnoreCase("0") != 0 || matId.compareToIgnoreCase("0") != 0) {
                    Tools.warning("Warning: skipping trio/duo child " + ID);
                    line = br.readLine();
                    continue;
                }

                BitSet bitSeq0 = new BitSet(getVariants().size());
                BitSet bitSeq1 = new BitSet(getVariants().size());
                BitSet mask0 = new BitSet(getVariants().size());
                mask0.flip(0, mask0.size());
                BitSet mask1 = new BitSet(getVariants().size());
                mask1.flip(0, mask1.size());
                int varIndex = 0;
                for (int i = 6; i < strSplit.length; i += 2) {
                    char geno0 = strSplit[i].charAt(0);
                    char geno1 = strSplit[i + 1].charAt(0);;
                    Variant thisVar = getVariants().get(varIndex);
                    cnt++;
                    if (geno0 == '0') {
                        // mask this position out
                        mask0.flip(varIndex);
                        miss++;
                    } else {
                        // record sequence by flipping allele
                        if (thisVar.alleleToBit(geno0)) {
                            bitSeq0.flip(varIndex);
                        }
                    }
                    if (geno1 == '0') {
                        // mask this position out
                        mask1.flip(varIndex);
                        miss++;
                    } else {
                        // record sequence by flipping allele
                        if (thisVar.alleleToBit(geno1)) {
                            bitSeq1.flip(varIndex);
                        }
                    }
                    varIndex++;
                }
                Individual ind0 = new Individual(famId, indId + ".0", patId, matId, sex, pheno, bitSeq0, mask0);
                Individual ind1 = new Individual(famId, indId + ".1", patId, matId, sex, pheno, bitSeq1, mask1);
                getIndividuals().put(famId + "\t" + indId + ".0", ind0);
                getIndividuals().put(famId + "\t" + indId + ".1", ind1);
                line = br.readLine();
                if (EstimateMutationRateFromIBD.maxThreads <= 1) {
                    if (!enterTainUser) {
                        Tools.printVerboseProgressSameLine("Read " + lineCnt + " individuals.");
                    } else {
                        Tools.printAnimation();
                    }
                }
            }
            Tools.printVerboseProgressLevel1("Read " + lineCnt + " individuals. Missing fraction: " + miss / (2. * cnt));
        } catch (IOException ex) {
            Tools.exit("Could not read ped file " + pedFile);
        }
        Tools.printVerboseProgressLevel1("Reading file " + mapFile);
        try {
            br = new BufferedReader(new FileReader(mapFile));
        } catch (FileNotFoundException ex) {
            Tools.exit("Could not open map file " + mapFile);
        }
        try {
            String line = br.readLine();
            int cnt = 0, nonIncreasing = 0;
            while (line != null) {
                cnt++;
                String[] strSplit = line.split("\\s+");
                if (cnt > 1 && getChr().compareTo(strSplit[0]) != 0) {
                    Tools.exit("found two chromosomes in map.");
                }
                setChr(strSplit[0]);
                getSNP_IDs().add(strSplit[1]);
                double gen = Double.parseDouble(strSplit[2]);
                if (cnt > 1 && gen < getGenToIndex().lastKey()) {
                    gen = getGenToIndex().lastKey();
                    nonIncreasing++;
//                    Tools.warning("Warning: found two non-increasing genetic positions in map: " + genToIndex.lastKey() + "\t" + gen + ". Substituting with largest value.");
                }
                if (getGenToIndex().get(gen) == null) { // otherwise there's an earlier SNP at this gen dist
                    getGenToIndex().put(gen, cnt - 1);
                }
                getIDtoPos().put(strSplit[1], cnt - 1);
                getIdToGen().put(strSplit[1], gen);
                getGenPosArray().add(gen);
                int phys = Integer.parseInt(strSplit[3]);
                getPhysToIndex().put(phys, cnt - 1);
                getPhysToID().put(phys, strSplit[1]);
                getGenToPhys().put(gen, phys);
                if (cnt > 1 && phys < getPhysPos().get(getPhysPos().size() - 1)) {
                    Tools.exit("found two non-increasing physical positions in map.");
                }
                getPhysPos().add(phys);
                line = br.readLine();
            }
            if (nonIncreasing > 0) {
                Tools.warning("WARNING: found " + nonIncreasing + " non-increasing genetic positions in map. Set them to higher value.");
            }
            setMap(new Map(getChr(), getSNP_IDs(), getGenToIndex(), getPhysPos()));
            getMap().setGenPosArray(getGenPosArray());
            Tools.printVerboseProgressLevel1("Read " + getMap().getSize() + " markers.");
        } catch (IOException ex) {
            Tools.exit("Could not read map file " + mapFile);
        }
    }

    public int genCoordToPhys(double pos) {

        double fromGen = (getGenToPhys().floorKey(pos) != null)
                ? getGenToPhys().floorKey(pos) : getGenToPhys().ceilingKey(pos);
        double fromPhys = (getGenToPhys().floorEntry(pos) != null)
                ? getGenToPhys().floorEntry(pos).getValue() : getGenToPhys().ceilingEntry(pos).getValue();
        if (fromGen == pos) {
            // exact match in map
            return getGenToPhys().floorEntry(pos).getValue();
        }

        double toGen = (getGenToPhys().ceilingKey(pos) != null)
                ? getGenToPhys().ceilingKey(pos) : getGenToPhys().floorKey(pos);
        double toPhys = (getGenToPhys().ceilingEntry(pos) != null)
                ? getGenToPhys().ceilingEntry(pos).getValue() : getGenToPhys().floorEntry(pos).getValue();
        double totPhysDist = toPhys - fromPhys;
        double totGenDist = toGen - fromGen;
        double fraction = (pos - fromGen) / totGenDist;
        return (int) Math.round(fromPhys + fraction * totPhysDist);

    }

    public void addFreqFile(String freqFile) {
        BufferedReader br = null;
        Tools.printVerboseProgressLevel1("Reading frequency file " + freqFile);
        try {
            br = new BufferedReader(new FileReader(freqFile));
        } catch (FileNotFoundException ex) {
            Tools.exit("Could not open frequency file " + freqFile);
        }
        try {
            String line = br.readLine();
            double cnt = 0;
            while (line != null) {
                String[] strSplit = line.trim().split("\\s+");
                if (strSplit[0].compareToIgnoreCase("CHR") == 0) {
                    line = br.readLine();
                    continue;
                }
                String ID = strSplit[1];
                double freq = Double.parseDouble(strSplit[4]);
                double variantCount = Double.parseDouble(strSplit[5]);
                int copies = (int) Math.round(freq * variantCount);
                getIDToFreq().put(ID, freq);
                getIDToVarCounts().put(ID, copies);
                line = br.readLine();
                cnt++;
            }
            Tools.printVerboseProgressLevel1("Read " + cnt + " frequency lines.");
        } catch (IOException ex) {
            Tools.exit("Could not read frequency file " + freqFile);
        }
    }

    public void readPosteriorsFile(String filename, double posteriorFrom, double posteriorTo, boolean storePosterior) throws IOException {
        BufferedReader br = null;
        try {
            InputStream fileStream = new FileInputStream(filename);
            InputStream gzipStream = new GZIPInputStream(fileStream);
            Reader decoder = new InputStreamReader(gzipStream, "UTF-8");
            br = new BufferedReader(decoder);
        } catch (FileNotFoundException ex) {
            Tools.exit("Could not open posteriors file " + filename);
        }
        try {
            HashSet<String> availableIndividuals = new HashSet<String>();
            for (String ind : getIndividuals().keySet()) {
                availableIndividuals.add(ind);
            }
            Tools.printVerboseProgressLevel1("Reading posterior file " + filename);
            if (storePosterior) {
                Tools.warning("Warning: Storing posteriors in memory (requires lots of memory)");
            }
            String line = br.readLine();
            int lineCnt = 0;
            double miss = 0;
            double cnt = 0;
            int length = (line.split("\\s+").length - 6) / 2;
            for (int i = 0; i < length; i++) {
                getVariants().add(new Variant());
            }
            ArrayList<String> IDorder = new ArrayList<String>();
            int currIndex = 0;
            int lastPos = -1;
            while (line != null) {
                String[] strSplit = line.split("\\s+");
                if (strSplit[0].startsWith("#")) {
                } else if (lineCnt == 0) {
                    // header will contain IDs
                    for (int i = 0; i < strSplit.length; i += 2) {
                        if (availableIndividuals.contains(strSplit[i] + "\t" + strSplit[i + 1])) {
                            availableIndividuals.remove(strSplit[i] + "\t" + strSplit[i + 1]);
                        }
                        IDorder.add(strSplit[i] + "\t" + strSplit[i + 1]);
                    }
                    if (availableIndividuals.size() != 0) {
                        Tools.exit("Coulnd not find all individuals in posterior file. Missing " + availableIndividuals.size() + ".");
                    }
                    lineCnt++;
                } else {
                    // body line
                    int count = 0;
                    int parsedPos = Integer.parseInt(strSplit[0]);
                    if (parsedPos <= lastPos) {
                        Tools.warning("Warning: skipping duplicate line or position " + parsedPos + ", which is less than previous in posterior file.");
                        line = br.readLine();
                        lastPos = parsedPos;
                        lineCnt++;
                        continue;
                    }
                    if (lastPos > -1 && parsedPos > this.getPhysPos().get(currIndex)) {
                        Tools.warning("Warning: Missing position " + this.getPhysPos().get(currIndex) + " in posterior file. Setting mask to 0 for everyone.");
                        for (Individual ind : getIndividuals().values()) {
                            ind.getMask().set(currIndex, false);
                        }
                        currIndex++;
                        lineCnt++;
                        continue;
                    }
                    if (!physToID.containsKey(parsedPos)) {
//                        System.out.println("Nope " + parsedPos);
                        line = br.readLine();
                        lineCnt++;
                        continue;
                    }
                    if (parsedPos != this.getPhysPos().get(currIndex) || currIndex != getPhysToIndex().get(parsedPos)) {
                        Tools.exit("Mismatch of index/position in posteriors. " + parsedPos + "\t"
                                + this.getPhysPos().get(currIndex) + "\t" + currIndex + "\t" + getPhysToIndex().get(parsedPos));
                    }
                    for (int i = 3; i < strSplit.length; i += 4) {
                        float probAA = Float.parseFloat(strSplit[i]);
                        float probRR = Float.parseFloat(strSplit[i + 1]);
                        float probAR = Float.parseFloat(strSplit[i + 2]);
                        float probRA = Float.parseFloat(strSplit[i + 3]);
                        float probHet = probAR + probRA;
                        float posterior;
                        if (probHet > probAA && probHet > probRR) {
                            posterior = Math.max(probAR, probRA);
                        } else {
                            posterior = Math.max(probAA, probRR);
                        }
                        if (!individuals.containsKey(IDorder.get(count))) {
                            if (!isWarnedPostNotContainsInd()) {
                                setWarnedPostNotContainsInd(true);
                                Tools.warning("Warning: Some individuals in the posteriors file are not present in the loaded data, e.g.: " + IDorder.get(count) + ". There could be more. Will not print this warning any longer.");
                            }
                        } else {
                            Individual ind = getIndividuals().get(IDorder.get(count));
                            if (!(posterior >= posteriorFrom && posterior <= posteriorTo)) {
                                ind.getMask().set(currIndex, false);
                            }
                            if (storePosterior) {
                                ind.getPosteriors().add(posterior);
                            }
                        }
                        count++;
                    }
                    currIndex++;
                    lineCnt++;
                    lastPos = parsedPos;
                }
                line = br.readLine();
            }
        } catch (IOException ex) {
            Tools.exit("Could not read posteriors file " + filename + ". This may be due to a formatting issue. Check that the first line of the file contains a space divided list of family ID and ID.");
        }
    }

    /**
     * @return the variants
     */
    public ArrayList<Variant> getVariants() {
        return variants;
    }

    /**
     * @param variants the variants to set
     */
    public void setVariants(ArrayList<Variant> variants) {
        this.variants = variants;
    }

    /**
     * @return the SNP_IDs
     */
    public ArrayList<String> getSNP_IDs() {
        return SNP_IDs;
    }

    /**
     * @param SNP_IDs the SNP_IDs to set
     */
    public void setSNP_IDs(ArrayList<String> SNP_IDs) {
        this.SNP_IDs = SNP_IDs;
    }

    /**
     * @return the genToIndex
     */
    public TreeMap<Double, Integer> getGenToIndex() {
        return genToIndex;
    }

    /**
     * @param genToIndex the genToIndex to set
     */
    public void setGenToIndex(TreeMap<Double, Integer> genToIndex) {
        this.genToIndex = genToIndex;
    }

    /**
     * @return the physToIndex
     */
    public TreeMap<Integer, Integer> getPhysToIndex() {
        return physToIndex;
    }

    /**
     * @param physToIndex the physToIndex to set
     */
    public void setPhysToIndex(TreeMap<Integer, Integer> physToIndex) {
        this.physToIndex = physToIndex;
    }

    /**
     * @return the genToPhys
     */
    public TreeMap<Double, Integer> getGenToPhys() {
        return genToPhys;
    }

    /**
     * @param genToPhys the genToPhys to set
     */
    public void setGenToPhys(TreeMap<Double, Integer> genToPhys) {
        this.genToPhys = genToPhys;
    }

    /**
     * @return the idToGen
     */
    public TreeMap<String, Double> getIdToGen() {
        return IDToGen;
    }

    /**
     * @param idToGen the idToGen to set
     */
    public void setIdToGen(TreeMap<String, Double> idToGen) {
        this.IDToGen = idToGen;
    }

    /**
     * @return the IDToFreq
     */
    public TreeMap<String, Double> getIDToFreq() {
        return IDToFreq;
    }

    /**
     * @param IDToFreq the IDToFreq to set
     */
    public void setIDToFreq(TreeMap<String, Double> IDToFreq) {
        this.IDToFreq = IDToFreq;
    }

    /**
     * @return the IDToVarCounts
     */
    public TreeMap<String, Integer> getIDToVarCounts() {
        return IDToVarCounts;
    }

    /**
     * @param IDToVarCounts the IDToVarCounts to set
     */
    public void setIDToVarCounts(TreeMap<String, Integer> IDToVarCounts) {
        this.IDToVarCounts = IDToVarCounts;
    }

    /**
     * @return the genPosArray
     */
    public ArrayList<Double> getGenPosArray() {
        return genPosArray;
    }

    /**
     * @param genPosArray the genPosArray to set
     */
    public void setGenPosArray(ArrayList<Double> genPosArray) {
        this.genPosArray = genPosArray;
    }

    /**
     * @return the physPos
     */
    public ArrayList<Integer> getPhysPos() {
        return physPos;
    }

    /**
     * @param physPos the physPos to set
     */
    public void setPhysPos(ArrayList<Integer> physPos) {
        this.physPos = physPos;
    }

    /**
     * @return the IDtoPos
     */
    public TreeMap<String, Integer> getIDtoPos() {
        return IDtoPos;
    }

    /**
     * @param IDtoPos the IDtoPos to set
     */
    public void setIDtoPos(TreeMap<String, Integer> IDtoPos) {
        this.IDtoPos = IDtoPos;
    }

    /**
     * @return the PhysToID
     */
    public TreeMap<Integer, String> getPhysToID() {
        return physToID;
    }

    /**
     * @param PhysToID the PhysToID to set
     */
    public void setPhysToID(TreeMap<Integer, String> PhysToID) {
        this.physToID = PhysToID;
    }

    /**
     * @return the individuals
     */
    public TreeMap<String, Individual> getIndividuals() {
        return individuals;
    }

    /**
     * @param individuals the individuals to set
     */
    public void setIndividuals(TreeMap<String, Individual> individuals) {
        this.individuals = individuals;
    }

    /**
     * @return the map
     */
    public Map getMap() {
        return map;
    }

    /**
     * @param map the map to set
     */
    public void setMap(Map map) {
        this.map = map;
    }

    /**
     * @return the includeMask
     */
    public Mask getIncludeMask() {
        return includeMask;
    }

    /**
     * @param includeMask the includeMask to set
     */
    public void setIncludeMask(Mask includeMask) {
        this.includeMask = includeMask;
    }

    /**
     * @return the excludeMask
     */
    public Mask getExcludeMask() {
        return excludeMask;
    }

    /**
     * @param excludeMask the excludeMask to set
     */
    public void setExcludeMask(Mask excludeMask) {
        this.excludeMask = excludeMask;
    }

    /**
     * @return the haveFreq
     */
    public boolean isHaveFreq() {
        return haveFreq;
    }

    /**
     * @param haveFreq the haveFreq to set
     */
    public void setHaveFreq(boolean haveFreq) {
        this.haveFreq = haveFreq;
    }

    /**
     * @return the roundTo
     */
    public double getRoundTo() {
        return roundTo;
    }

    /**
     * @return the totDiff
     */
    public double getTotDiff() {
        return totDiff;
    }

    /**
     * @param totDiff the totDiff to set
     */
    public void setTotDiff(double totDiff) {
        this.totDiff = totDiff;
    }

    /**
     * @return the totInMask
     */
    public double getTotInMask() {
        return totInMask;
    }

    /**
     * @param totInMask the totInMask to set
     */
    public void setTotInMask(double totInMask) {
        this.totInMask = totInMask;
    }

    /**
     * @return the inMask
     */
    public int getInMask() {
        return inMask;
    }

    /**
     * @param inMask the inMask to set
     */
    public void setInMask(int inMask) {
        this.inMask = inMask;
    }

    /**
     * @return the notInMask
     */
    public int getNotInMask() {
        return notInMask;
    }

    /**
     * @param notInMask the notInMask to set
     */
    public void setNotInMask(int notInMask) {
        this.notInMask = notInMask;
    }

    /**
     * @return the warnedPostNotContainsInd
     */
    public boolean isWarnedPostNotContainsInd() {
        return warnedPostNotContainsInd;
    }

    /**
     * @param warnedPostNotContainsInd the warnedPostNotContainsInd to set
     */
    public void setWarnedPostNotContainsInd(boolean warnedPostNotContainsInd) {
        this.warnedPostNotContainsInd = warnedPostNotContainsInd;
    }

    /**
     * @return the chr
     */
    public String getChr() {
        return chr;
    }

    /**
     * @param chr the chr to set
     */
    public void setChr(String chr) {
        this.chr = chr;
    }

}

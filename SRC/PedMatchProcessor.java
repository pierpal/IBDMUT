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
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.Reader;
import java.io.Writer;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.concurrent.Callable;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

/**
 *
 * @author Pier Palamara
 */
class PedMatchProcessor implements Callable<TreeMap<String, Results>> {

    private final String pedFile;
    private final boolean useExcludeMask;
    private final String maskFile;
    private final String excludeFile;
    private final String matchFile;
    private final boolean usePosteriors;
    private final double offsetCM;
    private final ArrayList<Integer> maxVarCountList;
    private final double minLen;
    private String posteriorsFile;
    private final double posteriorFrom;
    private final double posteriorTo;
    private final boolean storePosterior;
    private final boolean printMutMatchOut;
    private final boolean saveBin;
    private final String saveBinSuffix;
    private final boolean loadBin;
    private final boolean writeMutMatchOut;
    private final String loadBinSuffix;
    private final ArrayList<String> maskFiles;

    private static TreeSet<String> onlyIncludeSNPs;
    private static boolean useOnlyIncludeSNPs = false;

    private static boolean cumulativeMaAFRegression = true;

    public static void setCumulativeMaAFRegression(boolean value) throws IOException {
        cumulativeMaAFRegression = value;
    }
//    private final Dataset data;

    public static void setOnlyIncludeSNPs(String file) throws IOException {
        Tools.printVerboseProgressLevel1("Reading SNP inclusion file " + file);
        BufferedReader br = null;
        onlyIncludeSNPs = new TreeSet<String>();
        useOnlyIncludeSNPs = true;
        try {
            br = new BufferedReader(new FileReader(file));
        } catch (FileNotFoundException ex) {
            Tools.exit("Could not open SNP excusion file " + file);
        }
        String line = br.readLine();
        while (line != null) {
            String[] strSplit = line.trim().split("\\s+");
            String ID = strSplit[0];
            onlyIncludeSNPs.add(ID);
            line = br.readLine();
        }
        Tools.printVerboseProgressLevel1("Read " + onlyIncludeSNPs.size() + " SNPs");
    }

    public PedMatchProcessor(String pedFile, boolean useMask, boolean useExcludeMask, String maskFile, String excludeFile, String matchFile,
            double offsetCM, ArrayList<Integer> maxVarCountList, double minLen, boolean usePosteriors, String posteriorsFile, double posteriorFrom,
            double posteriorTo, boolean storePosterior, boolean printMutMatchOut,
            boolean saveBin, String saveBinSuffix, boolean loadBin, String loadBinSuffix,
            boolean writeMutMatchOut, ArrayList<String> maskFiles) {
        this.pedFile = pedFile;
        this.useExcludeMask = useExcludeMask;
        this.maskFile = maskFile;
        this.excludeFile = excludeFile;
        this.matchFile = matchFile;
        this.offsetCM = offsetCM;
        this.maxVarCountList = maxVarCountList;
        this.minLen = minLen;
        this.usePosteriors = usePosteriors;
        this.posteriorsFile = posteriorsFile;
        this.posteriorFrom = posteriorFrom;
        this.posteriorTo = posteriorTo;
        this.storePosterior = storePosterior;
        this.printMutMatchOut = printMutMatchOut;
        this.saveBin = saveBin;
        this.saveBinSuffix = saveBinSuffix;
        this.loadBin = loadBin;
        this.loadBinSuffix = loadBinSuffix;
        this.writeMutMatchOut = writeMutMatchOut;
        this.maskFiles = maskFiles;
    }

    public TreeMap<String, Results> call() throws FileNotFoundException, IOException {

        Dataset data = new DataLoader(loadBinSuffix, saveBinSuffix, pedFile, saveBin, loadBin, usePosteriors, posteriorFrom, posteriorTo, storePosterior).call();
        TreeMap<String, Results> resultsList = new TreeMap<String, Results>();
        for (String maskFile : maskFiles) {
            boolean useTrinucleotideContext = (TrinucleotideContext.context != null);
            boolean didAlreadyWarnMatchFileHasNoSNPs = false;

            for (Integer count : maxVarCountList) {
                resultsList.put(maskFile + "\t" + count, new Results(maskFile, count));
            }

            FileOutputStream MUToutput = null;
            Writer MUTwriter = null;

            boolean useMask = maskFile.compareToIgnoreCase("") != 0;

            if (useMask) {
                data.setIncludeMask(new Mask(maskFile, data.getChr()));
            }
            if (useExcludeMask) {
                data.setExcludeMask(new Mask(excludeFile, data.getChr()));
            }

            Tools.printVerboseProgressLevel1("Reading file " + matchFile);
            if (this.writeMutMatchOut) {
                File f = new File(maskFile);
                String mask = f.getName();
                MUToutput = new FileOutputStream(matchFile + "." + mask + ".mut.gz");
                MUTwriter = new OutputStreamWriter(new GZIPOutputStream(MUToutput), "UTF-8");
            }
            BufferedReader br = null;
            try {
                InputStream fileStream = new FileInputStream(matchFile);
                InputStream gzipStream = new GZIPInputStream(fileStream);
                Reader decoder = new InputStreamReader(gzipStream, "UTF-8");
                br = new BufferedReader(decoder);
            } catch (FileNotFoundException ex) {
                Tools.exit("Could not open match file " + matchFile);

            } catch (IOException ex) {
                Logger.getLogger(PedMatchProcessor.class
                        .getName()).log(Level.SEVERE, null, ex);
                System.exit(
                        1);
            }

            // PROCESS MATCH
            try {
                String line = br.readLine();
                while (line != null) {
                    ArrayList<Double> freqOfMismatchingSNPsOnSegment = new ArrayList<Double>();
                    ArrayList<String> IDsOfMismatchingSNPsOnSegment = new ArrayList<String>();
                    String[] strSplit = line.split("\\s+");
                    String FamID1 = strSplit[0];
                    String IndID1 = strSplit[1];
                    String FamID2 = strSplit[2];
                    String IndID2 = strSplit[3];
                    String chromosome = strSplit[4];
                    int fromPos = Integer.parseInt(strSplit[5]);
                    int toPos = Integer.parseInt(strSplit[6]);
                    String fromSNP = strSplit[7];
                    String toSNP = strSplit[8];
                    double length = Double.parseDouble(strSplit[10]);
                    String unit = strSplit[11];
                    if (unit.compareToIgnoreCase("Mb") == 0) {
                        Tools.warning("Warning: parsed match which is not in cM unit.");
                    }

                    if (length < minLen || length - 2 * offsetCM <= 0) {
                        line = br.readLine();
                        continue;
                    }

                    int mapFrom, mapTo;
                    try {
                        mapFrom = data.getIDtoPos().get(fromSNP);
                        mapTo = data.getIDtoPos().get(toSNP);
                    } catch (Exception ex) {
                        if (!didAlreadyWarnMatchFileHasNoSNPs) {
                            Tools.warning("Warning: Could not find SNPs indicated in match file. Will always find them using physical start/end of segments.");
                            didAlreadyWarnMatchFileHasNoSNPs = true;
                        }
                        fromSNP = (data.getPhysToID().ceilingEntry(fromPos) != null) ? data.getPhysToID().ceilingEntry(fromPos).getValue() : data.getPhysToID().floorEntry(fromPos).getValue();
                        mapFrom = data.getIDtoPos().get(fromSNP);
                        toSNP = (data.getPhysToID().floorEntry(toPos) != null) ? data.getPhysToID().floorEntry(toPos).getValue() : data.getPhysToID().ceilingEntry(toPos).getValue();
                        mapTo = data.getIDtoPos().get(toSNP);
                    }
                    Individual ID1 = data.getIndividuals().get(FamID1 + "\t" + IndID1);
                    Individual ID2 = data.getIndividuals().get(FamID2 + "\t" + IndID2);
                    TreeMap<Integer, Integer> diffList = new TreeMap<Integer, Integer>();
                    for (Integer maxVarCount : maxVarCountList) {
                        diffList.put(maxVarCount, 0);
                    }
                    BitSet bitSeq1 = ID1.getBitSeq();
                    BitSet mask1 = ID1.getMask();
                    BitSet bitSeq2 = ID2.getBitSeq();
                    BitSet mask2 = ID2.getMask();
                    double fromGen = data.getGenPosArray().get(mapFrom);
                    double toGen = data.getGenPosArray().get(mapTo);
                    int firstPos = -1;
                    int inMaskSize = 0;
                    for (int ind = mapFrom; ind <= mapTo; ind++) {
                        double distanceFromLeftEdge = data.getGenPosArray().get(ind) - fromGen;
                        double distanceFromRightEdge = toGen - data.getGenPosArray().get(ind);
                        if (distanceFromLeftEdge < offsetCM || distanceFromRightEdge < offsetCM) {
                            continue;
                        }
                        if (firstPos == -1) {
                            firstPos = data.getPhysPos().get(ind);
                        }
                        int SNPpos = data.getPhysPos().get(ind);

                        if ((useExcludeMask && data.getExcludeMask().contains(SNPpos)) || (useMask && !data.getIncludeMask().contains(SNPpos))) {
                            data.setNotInMask(data.getNotInMask() + 1);
                            continue;
                        }

                        if (useOnlyIncludeSNPs && !onlyIncludeSNPs.contains(data.getPhysToID().get(SNPpos))) {
                            data.setNotInMask(data.getNotInMask() + 1);
                            continue;
                        }

                        data.setInMask(data.getInMask() + 1);
                        if (!mask1.get(ind) || !mask2.get(ind)) {
                            //HANDLE MISSING
                        } else {
                            if (bitSeq1.get(ind) != bitSeq2.get(ind)) {
                                if (data.isHaveFreq() && (printMutMatchOut || writeMutMatchOut)) {
                                    // if they are the same
                                }
                                for (int cntInd = 0; cntInd < maxVarCountList.size(); cntInd++) {
                                    int maxVarCount = maxVarCountList.get(cntInd);
                                    Results results = resultsList.get(maskFile + "\t" + maxVarCount);
                                    boolean condition = (!cumulativeMaAFRegression)
                                            ? (cntInd != 0)
                                            && (data.getIDToVarCounts().get(data.getPhysToID().get(SNPpos)) < maxVarCount
                                            && data.getIDToVarCounts().get(data.getPhysToID().get(SNPpos)) >= maxVarCountList.get(cntInd - 1))
                                            : (data.getIDToVarCounts().get(data.getPhysToID().get(SNPpos)) <= maxVarCount);
                                    if (condition) {
                                        int diff = diffList.get(maxVarCount);
                                        diff++;
                                        if (useTrinucleotideContext) {
                                            Pair<String, String> trinucleotides = TrinucleotideContext.getTrinucleotides(data.getPhysToID().get(SNPpos));
                                            if (trinucleotides != null) {
                                                results.increaseTrinucleotideCounts(trinucleotides.getKey(), trinucleotides.getValue());
                                            }
                                        }
                                        diffList.put(maxVarCount, diff);
                                        double lenRounded = Math.floor(length * 2) / 2.0;
                                        if (data.getVariants().get(ind).isTransition()) {
                                            results.setTransitions(results.getTransitions() + 1);
                                            double val = (results.getTransitionPerCM().containsKey(lenRounded)) ? results.getTransitionPerCM().get(lenRounded) + 1 : 1;
                                            results.getTransitionPerCM().put(lenRounded, val);
                                        } else {
                                            results.setTransversions(results.getTransversions() + 1);
                                            double val = (results.getTransversionPerCM().containsKey(lenRounded)) ? results.getTransversionPerCM().get(lenRounded) + 1 : 1;
                                            results.getTransversionPerCM().put(lenRounded, val);
                                        }
                                        if (data.isHaveFreq()) {
                                            double frq = data.getIDToFreq().get(data.getPhysToID().get(SNPpos));
                                            int frqCounts = 1;
                                            if (results.getFreqCounts().containsKey(frq)) {
                                                frqCounts += results.getFreqCounts().get(frq);
                                            }
                                            results.getFreqCounts().put(frq, frqCounts);
                                        }
                                    }
                                }
                            } else {
                                // if they are different
                                double frq = data.getIDToFreq().get(data.getPhysToID().get(SNPpos));
                                IDsOfMismatchingSNPsOnSegment.add(data.getPhysToID().get(SNPpos));
                                freqOfMismatchingSNPsOnSegment.add(frq);
                            }
                        }
                    }
                    int exactFromPhys = data.genCoordToPhys(data.getIdToGen().get(fromSNP) + offsetCM);
                    int exactToPhys = data.genCoordToPhys(data.getIdToGen().get(toSNP) - offsetCM);
                    inMaskSize = (exactToPhys <= exactFromPhys) ? 0 : (exactToPhys - exactFromPhys);
                    if (useMask && inMaskSize > 0) {
                        inMaskSize = data.getIncludeMask().getMaskOverlapWithRegion(exactFromPhys, exactToPhys);
                    }
                    if (inMaskSize > 0) {
                        for (Integer maxVarCount : maxVarCountList) {
                            int diff = diffList.get(maxVarCount);
                            Results results = resultsList.get(maskFile + "\t" + maxVarCount);
                            results.addTotInMask(inMaskSize);
                            results.setTotalSegmentCount(results.getTotalSegmentCount() + 1);
                            results.setTotalSegmentLength(results.getTotalSegmentLength() + length);
//                        TODO here is where it's added. May fix'
                            results.addHetCount(length, diff, inMaskSize);
                            StringBuilder additionalInfoString = null;
                            if (printMutMatchOut || writeMutMatchOut) {
                                additionalInfoString = new StringBuilder();
                                for (String k : IDsOfMismatchingSNPsOnSegment) {
                                    additionalInfoString.append(k);
                                    additionalInfoString.append(" ");
                                }
                            }
                            if (printMutMatchOut) {
                                System.out.println(line + "\t" + diff + "\t" + inMaskSize + "\t" + additionalInfoString.toString());
                            }
                            if (writeMutMatchOut) {
                                MUTwriter.write(line + "\t" + diff + "\t" + inMaskSize + "\t" + maxVarCount + "\t" + additionalInfoString.toString() + "\n");
                            }
                            double histBin = (double) Math.round(length * data.getRoundTo()) / data.getRoundTo();
                            double counts = diff;
                            data.setTotDiff(data.getTotDiff() + diff);
                            if (results.getCountsHistogram().containsKey(histBin)) {
                                counts += results.getCountsHistogram().get(histBin);
                            }
                            results.getCountsHistogram().put(histBin, counts);

                            double inMaskThisSegment = inMaskSize;
                            data.setTotInMask(data.getTotInMask() + inMaskSize);
                            if (results.getInMaskHistogram().containsKey(histBin)) {
                                inMaskThisSegment += results.getInMaskHistogram().get(histBin);
                            }
                            results.getInMaskHistogram().put(histBin, inMaskThisSegment);
                        }
                    }
                    line = br.readLine();
                }
            } catch (IOException ex) {
                Tools.exit("Could not read match file " + matchFile);
            }

            if (this.writeMutMatchOut) {
                MUTwriter.flush();
                MUTwriter.close();
            }

            if (useMask) {
                Tools.printVerboseProgressLevel2("Finished analyzing file " + pedFile + ". In mask:\t" + data.getInMask() + "\tnot in mask:\t" + data.getNotInMask() + "\tratio:\t" + ((double) data.getInMask()) / data.getNotInMask());
            }
        }
        return resultsList;
    }

}

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
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.TreeMap;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

/**
 *
 * @author Pier Palamara
 */
public class EstimateMutationRateFromIBD {

    static int maxThreads = 1;
    public static final int largeValue = 100000000;

    public static TreeMap<String, Pair<Double, Double>> inverseHeterozygosity;
    private static TreeMap<String, Pair<Integer, Integer>> regionBoundaries = new TreeMap<String, Pair<Integer, Integer>>();

    public static void unitTestRegression() {
        double[][] x = {{-0.5}, {0.0}, {0.5}, {1.0}, {1.5}, {2.0}, {2.5}, {3.0}, {3.5}, {4.0}, {4.5}};
        double[] y = {1.566927e-08, 1.662871e-08, 1.469284e-08, 1.705565e-08, 1.909188e-08, 1.298154e-08, 2.251976e-08, 1.055619e-08, 1.245717e-08, 1.775048e-08, 1.805694e-08};
        double[] w = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};

        double reg[] = Tools.simpleRegression(x, y, w);
        for (int i = 0; i < reg.length; i++) {
            System.out.println(reg[i]);
        }

//        LinearRegression linReg = new LinearRegression();
//        double[] w = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};
////        double[][] xx = {{-0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5}};
//        double[][] xx = {{1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.}, {-0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5}};
////        double[][] xx = {{-0.5}, {0.0}, {0.5}, {1.0}, {1.5}, {2.0}, {2.5}, {3.0}, {3.5}, {4.0}, {4.5}};
//        linReg.regress(y, xx, w);
//        for (int i = 0; i < linReg.getCoefficients().length; i++) {
//            System.out.println(linReg.getCoefficients()[i]);
//        }
    }

    public static void unitTestDemography() throws IOException {
//        Demography demo = new Demography("GoNL.demo");
        Demography demo = new Demography("2000.simple");
        double from = 1.6, interval = .1, to = 5.;
        for (double u = from; u <= to; u += interval) {
            double age = demo.getAgeOfSegment(u, 1000000., true);
            System.out.println(u + "\t" + age + "\t" + 2 * age);
        }
        for (double u = from; u <= to; u += interval) {
            double age = demo.getAgeOfSegment(u, u + interval, false);
            System.out.println(u + "\t" + age + "\t" + 2 * age);
        }
    }

    public static void unitTestJackknife() {
        Double[] val = {-0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5};
        Double[] weights = {1.566927e-08, 1.662871e-08, 1.469284e-08, 1.705565e-08, 1.909188e-08, 1.298154e-08, 2.251976e-08, 1.055619e-08, 1.245717e-08, 1.775048e-08, 1.805694e-08};
        Double global = 1.0;
        Pair<Double, Double> res = Tools.weightedBlockJackknife(val, weights, global);
        System.out.println(res.getKey() + "\t" + res.getValue());
    }

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws Exception {

//        Tools.testColors();
//        unitTestDemography();
//        unitTestRegression();
//        unitTestJackknife();
//        if (true) {
//            return;
//        }

        String pedFile = "", listFile = "", matchFile = "", maskFile = "", excludeFile = "",
                posteriorsFile = "", outMatchFile = "", demographyFile = "", logFile = "",
                trinucleotideContextFile = "", maskListFile = "", regionBoundariesFile = "",
                onlyIncludeSNPsFile = "", weightsFile = "";
        String saveBinSuffix = "", loadBinSuffix = "";
        ArrayList<String> pedList = new ArrayList<String>();
        ArrayList<String> matchList = new ArrayList<String>();
        double posteriorFrom = 0.0, posteriorTo = 0.0;
        boolean useMask = false;
        boolean useMaskList = false;
        boolean useExcludeMask = false;
        boolean usePosteriors = false;
        boolean useMaxCount = false;
        boolean useMaxCountRangeRegression = false;
        boolean printMutMatchOut = false;
        boolean storePosterior = false;
        boolean haveDemography = false;
        boolean doJackKnife = false;
        boolean saveBin = false;
        boolean loadBin = false;
        boolean onlyIncludeSNPs = false;
        boolean trinucleotideNoDirection = true;
        boolean computeHeterozygosity = false;
        boolean writeMutMatchOut = false;
        double minimumLength = 2.0;
        double maximumLength = 5.0;
        Demography demography = null;
        double offsetCM = 0.0;
        double minLen = 0.0;
        int maxVarCountMin = largeValue;
        int maxVarCountInterval = largeValue;
        int maxVarCountMax = largeValue;

        int argIndex = 0;
        while (argIndex < args.length) {
            String arg = args[argIndex++];
            if (arg.equals("--plink")) {
                pedFile = args[argIndex];
                argIndex++;
            } else if (arg.equals("--printMutMatch")) {
                printMutMatchOut = true;
            } else if (arg.equals("--printRegression")) {
                LinearRegression.setPrintReg(true);
            } else if (arg.equals("--weightedTMRCARegression")) {
                Results.setUseJackknifeWeigthsIntMRCAReg(true);
            } else if (arg.equals("--printJackknife")) {
                Tools.setPrintJackknife(true);
            } else if (arg.equals("--writeMutMatch")) {
                writeMutMatchOut = true;
            } else if (arg.equals("--noCumulativeMaAF")) {
                PedMatchProcessor.setCumulativeMaAFRegression(false);
                Results.setCumulativeMaAFRegression(false);
            } else if (arg.equals("--directionalTrinucleotides")) {
                trinucleotideNoDirection = false;
            } else if (arg.equals("--trinucleotideContext")) {
                trinucleotideContextFile = args[argIndex];
                argIndex++;
            } else if (arg.equals("--saveBin")) {
                saveBin = true;
                saveBinSuffix = args[argIndex];
                argIndex++;
            } else if (arg.equals("--loadBin")) {
                loadBin = true;
                loadBinSuffix = args[argIndex];
                argIndex++;
            } else if (arg.equals("--plinkList")) {
                listFile = args[argIndex];
                argIndex++;
            } else if (arg.equals("--quiet")) {
                Tools.makeQuiet();
            } else if (arg.equals("--threads")) {
                maxThreads = Integer.parseInt(args[argIndex]);;
                argIndex++;
            } else if (arg.equals("--weightMinObs")) {
                Results.setMinimumNumberOfPointsToComputeAWeight(Integer.parseInt(args[argIndex]));
                argIndex++;
            } else if (arg.equals("--MaAF")) {
                if (useMaxCountRangeRegression) {
                    Tools.printHelpAndExit("Cannot use both --maxCountRangeRegression and --maxCount");
                }
                useMaxCount = true;
                maxVarCountMin = Integer.parseInt(args[argIndex]);
                maxVarCountInterval = maxVarCountMin;
                maxVarCountMax = maxVarCountMin;
                argIndex++;
            } else if (arg.equals("--MaAFRegression")) {
                if (useMaxCount) {
                    Tools.printHelpAndExit("Cannot use both --maxCountRangeRegression and --maxCount");
                }
                useMaxCountRangeRegression = true;
                maxVarCountMin = Integer.parseInt(args[argIndex]);
                argIndex++;
                maxVarCountInterval = Integer.parseInt(args[argIndex]);
                argIndex++;
                maxVarCountMax = Integer.parseInt(args[argIndex]);
                argIndex++;
            } else if (arg.equals("--match")) {
                matchFile = args[argIndex];
                argIndex++;
            } else if (arg.equals("--offsetCM")) {
                offsetCM = Double.parseDouble(args[argIndex]);
                argIndex++;
            } else if (arg.equals("--roundTo")) {
                Dataset.setRoundTo(Double.parseDouble(args[argIndex]));
                argIndex++;
            } else if (arg.equals("--mask")) {
                useMask = true;
                maskFile = args[argIndex];
                argIndex++;
            } else if (arg.equals("--demography")) {
                haveDemography = true;
                demographyFile = args[argIndex];
                argIndex++;
            } else if (arg.equals("--exclude")) {
                useExcludeMask = true;
                excludeFile = args[argIndex];
                argIndex++;
            } else if (arg.equals("--setTMRCAWeights")) {
                Results.setUseJackknifeWeigthsIntMRCAReg(true);
                weightsFile = args[argIndex];
                argIndex++;
            } else if (arg.equals("--onlyIncludeSNPs")) {
                onlyIncludeSNPs = true;
                onlyIncludeSNPsFile = args[argIndex];
                argIndex++;
            } else if (arg.equals("--log")) {
                Tools.setLoggerIsOn(true);
                logFile = args[argIndex];
                Tools.initLog(logFile);
                argIndex++;
            } else if (arg.equals("--maskList")) {
                useMaskList = true;
                maskListFile = args[argIndex];
                argIndex++;
            } else if (arg.equals("--storePosteriors")) {
                storePosterior = true;
            } else if (arg.equals("--nonOverlappingIntervals")) {
                storePosterior = true;
            } else if (arg.equals("--entertainUser")) {
                Dataset.setEntertainment();
            } else if (arg.equals("--computeHeterozygosity")) {
                computeHeterozygosity = true;
                regionBoundariesFile = args[argIndex];
                argIndex++;
            } else if (arg.equals("--minLength")) {
                minimumLength = Double.parseDouble(args[argIndex]);
                argIndex++;
            } else if (arg.equals("--lenRange")) {
                minimumLength = Double.parseDouble(args[argIndex]);
                minLen = minimumLength;
                argIndex++;
                maximumLength = Double.parseDouble(args[argIndex]);
                argIndex++;
            } else if (arg.equals("--posteriorRange")) {
                usePosteriors = true;
                posteriorFrom = Double.parseDouble(args[argIndex]);
                argIndex++;
                posteriorTo = Double.parseDouble(args[argIndex]);
                argIndex++;
            } else if (arg.equals("--jackknife")) {
                doJackKnife = true;
            } else if (arg.equals("--noCumulative")) {
                Results.setUseCumulative(false);
            } else {
                Tools.printHelpAndExit("Unsupported argument " + arg + ".");
            }
        }

        if (haveDemography) {
            demography = new Demography(demographyFile);
        }

        boolean useTrinucleotide = !"".equals(trinucleotideContextFile);
        if (useTrinucleotide) {
            TrinucleotideContext.setTrinucleotideContext(trinucleotideContextFile);
        }

        if (pedFile.compareToIgnoreCase("") == 0 && matchFile.compareToIgnoreCase("") == 0 && listFile.compareToIgnoreCase("") == 0) {
            Tools.printHelpAndExit("Both --plink/--match and --plinkList (or neither) were used.");
        }
        if (!haveDemography) {
            Tools.printHelpAndExit("Did not specify a demography.");
        }

        if (onlyIncludeSNPs) {
            PedMatchProcessor.setOnlyIncludeSNPs(onlyIncludeSNPsFile);
        }

        if (!"".equals(weightsFile)) {
            BufferedReader br = null;
            try {
                br = new BufferedReader(new FileReader(weightsFile));
            } catch (FileNotFoundException ex) {
                Tools.exit("Could not open tMRCA regression weights file " + weightsFile);
            }
            String line = br.readLine();
            TreeMap<Double, Double> weightsMap = new TreeMap<Double, Double>();
            int cnt = 0;
            while (line != null) {
                String[] tokens = line.split("\\s+");
                weightsMap.put(Double.parseDouble(tokens[0]), Double.parseDouble(tokens[1]));
                line = br.readLine();
                cnt++;
            }
            Results.setTMRCARegressionWeightsMap(weightsMap);
            Tools.printVerboseProgressLevel1("Read " + cnt + " tMRCA weights.");
        }

        if (computeHeterozygosity) {
            BufferedReader br = null;
            try {
                br = new BufferedReader(new FileReader(regionBoundariesFile));
            } catch (FileNotFoundException ex) {
                Tools.exit("Could not open region boundaries file " + regionBoundariesFile);
            }
            String line = br.readLine();
            while (line != null) {
                String[] tokens = line.split("\\s+");
                String thisFile = tokens[0];
                int from = Integer.parseInt(tokens[1]);
                int to = Integer.parseInt(tokens[2]);
                regionBoundaries.put(thisFile, new Pair(from, to));
                line = br.readLine();
            }
        }

        if (pedFile.compareToIgnoreCase("") != 0) {
            pedList.add(pedFile);
            matchList.add(matchFile);
        } else if (listFile.compareToIgnoreCase("") != 0) {
            BufferedReader br = null;
            try {
                br = new BufferedReader(new FileReader(listFile));
            } catch (FileNotFoundException ex) {
                Tools.exit("Could not open list file " + listFile);
            }
            String line = br.readLine();
            while (line != null) {
                String[] files = line.split("\\s+");
                pedFile = files[0];
                matchFile = files[1];
                pedList.add(pedFile);
                matchList.add(matchFile);
                line = br.readLine();
            }
        } else {
            Tools.printHelpAndExit("Did not set ped file or list.");
        }

        if (pedList.size() < 2 && doJackKnife) {
            Tools.exit("Trying to run Jackknife with one region.");
        }

        // RUN THREADS THAT READ DATA AND COMPUTE HISTOGRAMS FOR THE MASK, FOR EACH MAXCOUNT VALUE
        ArrayList<Integer> maxVarCountList = new ArrayList<Integer>();

        for (int count = maxVarCountMin; count <= maxVarCountMax; count += maxVarCountInterval) {
            maxVarCountList.add(count);
        }

        ArrayList<String> maskFiles = new ArrayList<String>();
        if (!useMask && !useMaskList) {
            maskFiles.add("");
        } else {
            if (useMaskList && useMask) {
                Tools.printHelpAndExit("--mask and --maskList should not be used together");
            }
            if (useMask) {
                maskFiles.add(maskFile);
            }
            if (useMaskList) {
                Tools.printProgress("Reading mask list file " + maskListFile);
                BufferedReader br = null;
                Tools.printVerboseProgressLevel1("Reading mask list file " + maskListFile);
                try {
                    br = new BufferedReader(new FileReader(maskListFile));
                } catch (FileNotFoundException ex) {
                    Tools.exit("Could not open mask list file " + maskListFile);
                }
                try {
                    String line = br.readLine();
                    int cnt = 0;
                    while (line != null) {
                        cnt++;
                        maskFiles.add(line);
                        line = br.readLine();
                    }
                    Tools.printVerboseProgressLevel1("Read " + cnt + " mask files.");
                } catch (IOException ex) {
                    Tools.exit("Could not read mask list file " + maskListFile);
                }
            }
        }

        if (computeHeterozygosity) {
            // NOTE: this is loading frequencies again, should be done in data loader. Fine for now.
            inverseHeterozygosity = getInverseHeterozygosityForAllMasks(maskFiles, pedList);
        }

        ArrayList<Future<TreeMap<String, Results>>> resultsList = new ArrayList<Future<TreeMap<String, Results>>>();
        ExecutorService executor = Executors.newFixedThreadPool(maxThreads);

        if (maxVarCountList.size() > 1) {
            // add results with no constraint, used to get difference between GC correction and uncorrected.
            maxVarCountList.add(EstimateMutationRateFromIBD.largeValue);
        }
        for (int i = 0; i < pedList.size(); i++) {
            pedFile = pedList.get(i);
            matchFile = matchList.get(i);
            if (Tools.loggerIsOn()) {
                Tools.writeLog("Analyzing files in " + pedFile + "\tmatch:\t" + matchFile + "\tmask:\t" + maskFile
                        + "\tminVarCount:\t" + maxVarCountMin + "\tmaxVarCountMin:\t" + maxVarCountMin + "\tmaxVarCountMin:\t" + maxVarCountMin);
            }
            Callable thread = new PedMatchProcessor(pedFile, useMask, useExcludeMask, maskFile, excludeFile, matchFile,
                    offsetCM, maxVarCountList, minLen, usePosteriors, posteriorsFile, posteriorFrom,
                    posteriorTo, storePosterior, printMutMatchOut,
                    saveBin, saveBinSuffix, loadBin, loadBinSuffix, writeMutMatchOut, maskFiles);
            resultsList.add(executor.submit(thread));

        }
        executor.shutdown();
        executor.awaitTermination(Integer.MAX_VALUE, TimeUnit.DAYS);
        Tools.finalizeLogIfOpen();

        for (String mask : maskFiles) {
            maskFile = mask;
            Results.processResults(doJackKnife, maxVarCountMin, maxVarCountInterval, maxVarCountMax, maxVarCountList,
                    maskFile, resultsList, haveDemography, useMaxCountRangeRegression, useTrinucleotide,
                    demography, minimumLength, maximumLength, trinucleotideNoDirection);
        }
    }

    private static TreeMap<String, Pair<Double, Double>> getInverseHeterozygosityForAllMasks(ArrayList<String> maskFiles, ArrayList<String> pedFiles) {
        TreeMap<String, Pair<Double, Double>> results = new TreeMap<String, Pair<Double, Double>>();
        TreeMap<String, TreeMap<String, Pair<Double, Double>>> heterozygosityForEachFileAndMask = new TreeMap<String, TreeMap<String, Pair<Double, Double>>>();
        for (String pedFile : pedFiles) {
            TreeMap<String, Pair<Double, Double>> heterozygosityForEachMask = new TreeMap<String, Pair<Double, Double>>();
            String freqFile = pedFile + ".frq";
            BufferedReader br = null;
            Tools.printVerboseProgressLevel1("Reading frequency file " + freqFile);
            try {
                br = new BufferedReader(new FileReader(freqFile));
            } catch (FileNotFoundException ex) {
                Tools.exit("Could not open frequency file " + freqFile);
            }
            try {
                String chr = null;
                String line = br.readLine();
                ArrayList<Double> freq = new ArrayList<Double>();
                ArrayList<Integer> pos = new ArrayList<Integer>();
                while (line != null) {
                    String[] strSplit = line.trim().split("\\s+");
                    if (strSplit[0].compareToIgnoreCase("CHR") == 0) {
                        line = br.readLine();
                        continue;
                    }
                    String ID = strSplit[1];
                    freq.add(Double.parseDouble(strSplit[4]));
                    String[] posString = strSplit[1].split(":"); // NOTE: SNP IDs must be coded as chr:physpos
                    pos.add(Integer.parseInt(posString[1]));
                    chr = strSplit[0];
                    line = br.readLine();
                }
                for (String maskFile : maskFiles) {
                    Mask mask = new Mask(maskFile, chr);
                    Pair<Integer, Integer> region = regionBoundaries.get(pedFile);
                    double numerator = 0, denominator = mask.getSize(region.getKey(), region.getValue());
                    for (int i = 0; i < pos.size(); i++) {
                        if (mask.contains(pos.get(i))) {
                            numerator += 2 * freq.get(i) * (1. - freq.get(i));
                        }
                    }
                    heterozygosityForEachMask.put(maskFile, new Pair(numerator, denominator));
                }
                heterozygosityForEachFileAndMask.put(pedFile, heterozygosityForEachMask);
            } catch (IOException ex) {
                Tools.exit("Could not read frequency file " + freqFile);
            }
        }
        // compute jackknife heterozygosity for each mask
        for (String maskFile : maskFiles) {
            Double mergedAll = 0.;
            Double[] weights = new Double[pedFiles.size()];
            Double[] allBut = new Double[pedFiles.size()];
            double mergedNum = 0, mergedDen = 0;
            int i = 0;
            double[] numerators = new double[pedFiles.size()];
            for (String pedFile : pedFiles) {
                // merge all
                numerators[i] = heterozygosityForEachFileAndMask.get(pedFile).get(maskFile).getKey();
                mergedNum += numerators[i];
                weights[i] = heterozygosityForEachFileAndMask.get(pedFile).get(maskFile).getValue();
                mergedDen += weights[i];
                i++;
            }
            mergedAll = 1. / (mergedNum / mergedDen);
            for (int j = 0; j < i; j++) {
                double numerator = 0, denominator = 0;
                // exclude j-th element
                for (int k = 0; k < i; k++) {
                    if (k != j) {
                        numerator += numerators[k];
                        denominator += weights[k];
                    }
                }
                allBut[j] = 1. / (numerator / denominator);
            }
            Pair<Double, Double> heterozygosityForThisMask = Tools.weightedBlockJackknife(allBut, weights, mergedAll);
            results.put(maskFile, heterozygosityForThisMask);
        }
        return results;
    }

}

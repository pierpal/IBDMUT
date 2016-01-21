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

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.TreeMap;
import java.util.concurrent.Future;

/**
 *
 * @author Pier Palamara
 */
public class Results {

    private static boolean useCumulative = true;
    private static int significantDigitsInLengthForRegression = 1;
    private static double binLengthInrease = 1 / Math.pow(10., significantDigitsInLengthForRegression);
    private static int significantDigitsInLengthForSpectrum = 1;
    private static int significantDigitsInLengthForTsTv = 1;
    private static boolean cumulativeMaAFRegression = true;

    private static boolean useWeigthsIntMRCAReg = false;
    private static int minimumNumberOfPointsToComputeAWeight = 5;
    private static final boolean printWeights = false;

    private static double[] tMRCARegressionWeights;

    private static TreeMap<Double, Double> tMRCARegressionWeightsMap;

    public static void setMinimumNumberOfPointsToComputeAWeight(int value) {
        minimumNumberOfPointsToComputeAWeight = value;
    }

    public static void setTMRCARegressionWeightsMap(TreeMap<Double, Double> map) {
        tMRCARegressionWeightsMap = map;
    }

    public static void setUseJackknifeWeigthsIntMRCAReg(boolean value) {
        useWeigthsIntMRCAReg = value;
    }

    public static boolean getUseWeigthsIntMRCAReg() {
        return useWeigthsIntMRCAReg;
    }

    public static double[] getTMRCARegressionWeights() {
        return tMRCARegressionWeights;
    }

    public static void setCumulativeMaAFRegression(boolean value) throws IOException {
        cumulativeMaAFRegression = value;
    }

    /**
     * @return the useCumulative
     */
    public static boolean isUseCumulative() {
        return useCumulative;
    }

    /**
     * @param aUseCumulative the useCumulative to set
     */
    public static void setUseCumulative(boolean aUseCumulative) {
        useCumulative = aUseCumulative;
    }

    /**
     * @return the significantDigitsInLengthForRegression
     */
    public static int getSignificantDigitsInLengthForRegression() {
        return significantDigitsInLengthForRegression;
    }

    public static double getLengthIncrement() {
        return binLengthInrease;
    }

    /**
     * @return the significantDigitsInLengthForSpectrum
     */
    public static int getSignificantDigitsInLengthForSpectrum() {
        return significantDigitsInLengthForSpectrum;
    }

    /**
     * @return the significantDigitsInLengthForTsTv
     */
    public static int getSignificantDigitsInLengthForTsTv() {
        return significantDigitsInLengthForTsTv;
    }

    /**
     * @param aSignificantDigitsInLengthForRegression the
     * significantDigitsInLengthForRegression to set
     */
    public static void setSignificantDigitsInLengthForRegression(int aSignificantDigitsInLengthForRegression) {
        significantDigitsInLengthForRegression = aSignificantDigitsInLengthForRegression;
    }

    /**
     * @param aSignificantDigitsInLengthForSpectrum the
     * significantDigitsInLengthForSpectrum to set
     */
    public static void setSignificantDigitsInLengthForSpectrum(int aSignificantDigitsInLengthForSpectrum) {
        significantDigitsInLengthForSpectrum = aSignificantDigitsInLengthForSpectrum;
    }

    /**
     * @param aSignificantDigitsInLengthForTsTv the
     * significantDigitsInLengthForTsTv to set
     */
    public static void setSignificantDigitsInLengthForTsTv(int aSignificantDigitsInLengthForTsTv) {
        significantDigitsInLengthForTsTv = aSignificantDigitsInLengthForTsTv;
    }

    private static void setWeightsFromMap(TreeMap<Double, BasePairsAndHets> hist) {
        Tools.printProgress("Using user-specified weights in tMRCA regression.");
        int indexWeights = 0;
        for (Double u : hist.keySet()) {
            tMRCARegressionWeights[indexWeights] = (tMRCARegressionWeightsMap.containsKey(u))
                    ? tMRCARegressionWeightsMap.get(u)
                    : 0.;
            if (printWeights) {
                System.out.println(u + "\tSet to " + tMRCARegressionWeights[indexWeights]);
            }
            indexWeights++;
        }
    }

    static void processResults(boolean doJackKnife, int maxVarCountMin, int maxVarCountInterval, int maxVarCountMax, ArrayList<Integer> maxVarCountList, String maskFile, ArrayList<Future<TreeMap<String, Results>>> resultsList, boolean haveDemography, boolean useMaxCountRangeRegression, boolean useTrinucleotide, Demography demography, double minimumLength, double maximumLength, boolean trinucleotideNoDirection) throws Exception {
        if (doJackKnife) {
            Tools.printProgress("Computing Jackknife results.");
            boolean applyingGeneConversionCorrection = false;
            TreeMap<Integer, Results> allMergedMap = new TreeMap<Integer, Results>();
            for (int count : maxVarCountList) {
                Results mergedAllThisCount = new Results(maskFile, count);
                for (int i = 0; i < resultsList.size(); i++) {
                    String key = maskFile + "\t" + count;
                    Results toBeMerged = resultsList.get(i).get().get(key);
                    mergedAllThisCount.mergeResults(toBeMerged);
                }
                allMergedMap.put(count, mergedAllThisCount);
            }

            if (allMergedMap.size() > 1) {
                applyingGeneConversionCorrection = true;
            } else {
                for (Results r : allMergedMap.values()) {
                    Tools.printResult("Ts/Tv:\t" + r.getTransitions() / r.getTransversions());
                }
            }

            if (applyingGeneConversionCorrection) {
                Tools.printProgress("Performing regression on maximum allele count frequency.");
            }
            ArrayList<TreeMap<Integer, Results>> allButOneList = new ArrayList<TreeMap<Integer, Results>>();
            for (int exclude = 0; exclude < resultsList.size(); exclude++) {
                TreeMap<Integer, Results> allButOneMap = new TreeMap<Integer, Results>();
                for (int count : maxVarCountList) {
                    Results mergedAllButOneThisCount = new Results(maskFile, count);
                    for (int i = 0; i < resultsList.size(); i++) {
                        if (i == exclude) {
                            continue;
                        }
                        String key = maskFile + "\t" + count;
                        Results toBeMerged = resultsList.get(i).get().get(key);
                        mergedAllButOneThisCount.mergeResults(toBeMerged);
                    }
                    allButOneMap.put(count, mergedAllButOneThisCount);
                }
                allButOneList.add(allButOneMap);
            }

            // compute weights for regression
            if (useWeigthsIntMRCAReg) {
                String key = maskFile + "\t" + EstimateMutationRateFromIBD.largeValue;
                TreeMap<Double, BasePairsAndHets> allHist = (isUseCumulative())
                        ? allMergedMap.get(EstimateMutationRateFromIBD.largeValue).getCumulativeHistogram()
                        : allMergedMap.get(EstimateMutationRateFromIBD.largeValue).getHistogram();
                tMRCARegressionWeights = new double[allHist.keySet().size()];

                if (tMRCARegressionWeightsMap != null) {
                    setWeightsFromMap(allHist);
                } else {
                    Tools.printProgress("Computing jackknife weights for inverse-variance weighting of tMRCA regression.");
                    int indexWeights = 0;
                    for (Double u : allHist.keySet()) {
                        // for each length bin, cumulative or not
                        int numerOfRegionsWhereThisLengthWasObserved = 0;
                        // this code is a bit slow, since getCumulativeHistogram() costs something.
                        for (int i = 0; i < resultsList.size(); i++) {
                            TreeMap<Double, BasePairsAndHets> thisHist = (isUseCumulative())
                                    ? allButOneList.get(i).get(EstimateMutationRateFromIBD.largeValue).getCumulativeHistogram()
                                    : allButOneList.get(i).get(EstimateMutationRateFromIBD.largeValue).getHistogram();
                            if (thisHist.containsKey(u)) {
                                numerOfRegionsWhereThisLengthWasObserved++;
                                if (printWeights) {
                                    System.out.println(i + "\tcontains\t" + u + "\t" + thisHist.get(u).count);
                                }
                            } else {
                                if (isUseCumulative() && thisHist.ceilingEntry(u) != null) {
                                    numerOfRegionsWhereThisLengthWasObserved++;
                                    if (printWeights) {
                                        System.out.println(i + "\tcontains\t" + u + "\t" + thisHist.ceilingEntry(u).getValue().count);
                                    }
                                } else {
                                    if (printWeights) {
                                        System.out.println(i + "\tdoes not contain\t" + u);
                                    }

                                }
                            }
                        }
                        if (printWeights) {
                            System.out.println(u + "\t" + numerOfRegionsWhereThisLengthWasObserved);
                        }
                        if (numerOfRegionsWhereThisLengthWasObserved < minimumNumberOfPointsToComputeAWeight) {
                            // skip this length bin if not enough observations
                            tMRCARegressionWeights[indexWeights] = Double.MAX_VALUE;
                            tMRCARegressionWeights[indexWeights] = 0.;
                            if (printWeights) {
                                System.out.println(u + "\tSet to " + tMRCARegressionWeights[indexWeights]);
                            }
                            indexWeights++;
                            continue;
                        }
                        // being lazy, looping twice
                        Double[] hetAllBut = new Double[numerOfRegionsWhereThisLengthWasObserved];
                        Double[] hetAllButWeights = new Double[numerOfRegionsWhereThisLengthWasObserved];
                        Double hetAll = allHist.get(u).mutRate / (double) allHist.get(u).count;
                        int cnt = 0;
                        for (int i = 0; i < resultsList.size(); i++) {
                            TreeMap<Double, BasePairsAndHets> thisHist = (isUseCumulative())
                                    ? allButOneList.get(i).get(EstimateMutationRateFromIBD.largeValue).getCumulativeHistogram()
                                    : allButOneList.get(i).get(EstimateMutationRateFromIBD.largeValue).getHistogram();
                            if (thisHist.containsKey(u)) {
                                hetAllBut[cnt] = thisHist.get(u).mutRate / (double) thisHist.get(u).count;
                                System.out.println("SEEN\t" + u + "\t" + hetAllBut[cnt]);
                                hetAllButWeights[cnt] = resultsList.get(i).get().get(key).getTotalBasePairsObserved();
                                System.out.println(hetAll + "\t" + hetAllBut[cnt] + "\t" + hetAllButWeights[cnt]);
                                cnt++;
                            } else {
                                if (isUseCumulative() && thisHist.ceilingEntry(u) != null) {
                                    hetAllBut[cnt] = thisHist.ceilingEntry(u).getValue().mutRate / (double) thisHist.ceilingEntry(u).getValue().count;
                                    hetAllButWeights[cnt] = resultsList.get(i).get().get(key).getTotalBasePairsObserved();
                                    cnt++;
                                }
                            }
                        }
                        Pair<Double, Double> jackRes = Tools.weightedBlockJackknife(hetAllBut, hetAllButWeights, hetAll);
                        tMRCARegressionWeights[indexWeights] = jackRes.getValue();
                        tMRCARegressionWeights[indexWeights] = 1 / (jackRes.getValue() * jackRes.getValue());
                        if (printWeights) {
                            System.out.println(u + "\tSet to " + tMRCARegressionWeights[indexWeights] + "\testimate\t" + jackRes.getKey() + "\ts.e.\t" + jackRes.getValue());
                        }
                        indexWeights++;
                    }
                }
            }

            // get results
            Double[] muMergedAllBut = new Double[resultsList.size()];
            Double[] secondParamMergedAllBut = new Double[resultsList.size()];
            Double[] noGCcorrectionMuMergedAllBut = new Double[resultsList.size()];
            Double[] noGCcorrectionMuDiffMergedAllBut = new Double[resultsList.size()];
            for (int exclude = 0; exclude < resultsList.size(); exclude++) {
                TreeMap<Integer, Results> allButOneMap = allButOneList.get(exclude);
                double[] muResAllButOne = (applyingGeneConversionCorrection) ? Results.regressMuForListOfMaxCounts(exclude, allButOneMap, maskFile, maxVarCountMin, maxVarCountInterval, maxVarCountMax, demography, minimumLength, maximumLength) : Results.regressMuNoMaxCounts(allButOneMap, maskFile, demography, minimumLength, maximumLength);
                muMergedAllBut[exclude] = muResAllButOne[0];
                secondParamMergedAllBut[exclude] = muResAllButOne[1];
                if (applyingGeneConversionCorrection) {
                    if (LinearRegression.getPrintReg()) {
                        System.out.println("Regression\tExcluding region " + exclude + "\tMaAF\t" + EstimateMutationRateFromIBD.largeValue);
                    }
                    noGCcorrectionMuMergedAllBut[exclude] = allButOneMap.get(EstimateMutationRateFromIBD.largeValue).getMu(demography, minimumLength, maximumLength)[1];
                    noGCcorrectionMuDiffMergedAllBut[exclude] = noGCcorrectionMuMergedAllBut[exclude] - muMergedAllBut[exclude];
                }
            }

            double[] muRes = (applyingGeneConversionCorrection)
                    ? Results.regressMuForListOfMaxCounts(-1, allMergedMap, maskFile,
                            maxVarCountMin, maxVarCountInterval, maxVarCountMax,
                            demography, minimumLength, maximumLength)
                    : Results.regressMuNoMaxCounts(allMergedMap, maskFile,
                            demography, minimumLength, maximumLength);
            Double muMergedAll = muRes[0];
            Double secondParamMergedAll = muRes[1];
            Double noGCcorrectionMuAll = null;
            Double noGCcorrectionMuDiffAll = null;
            if (applyingGeneConversionCorrection) {
                if (LinearRegression.getPrintReg()) {
                    System.out.println("Regression\tExcluding region " + (-1) + "\tMaAF\t" + EstimateMutationRateFromIBD.largeValue);
                }
                noGCcorrectionMuAll = allMergedMap.get(EstimateMutationRateFromIBD.largeValue).getMu(demography, minimumLength, maximumLength)[1];
                noGCcorrectionMuDiffAll = noGCcorrectionMuAll - muMergedAll;
            }
            Double[] mergedAllButWeights = new Double[resultsList.size()];
            for (int exclude = 0; exclude < resultsList.size(); exclude++) {
                double weight = 0.0;
                for (int i = 0; i < resultsList.size(); i++) {
                    if (i == exclude) {
                        continue;
                    }
                    String key = maskFile + "\t" + maxVarCountMin;
                    Results toBeMerged = resultsList.get(i).get().get(key);
                    weight += toBeMerged.getTotalBasePairsObserved();
                }
                mergedAllButWeights[exclude] = weight;
            }
            Pair<Double, Double> finalMuAndStd = Tools.weightedBlockJackknife(muMergedAllBut, mergedAllButWeights, muMergedAll);
            Pair<Double, Double> finalSecondParamAndStd = Tools.weightedBlockJackknife(secondParamMergedAllBut, mergedAllButWeights, secondParamMergedAll);
            Pair<Double, Double> finalMuAndStdUncorrected = null;
            Pair<Double, Double> finalMuDiffAndStdUncorrected = null;
            if (applyingGeneConversionCorrection) {
                finalMuAndStdUncorrected = Tools.weightedBlockJackknife(noGCcorrectionMuMergedAllBut, mergedAllButWeights, noGCcorrectionMuAll);
                finalMuDiffAndStdUncorrected = Tools.weightedBlockJackknife(noGCcorrectionMuDiffMergedAllBut, mergedAllButWeights, noGCcorrectionMuDiffAll);
            }
            String resultString = "Estimated mutation rate for mask " + maskFile + "\t" + finalMuAndStd.getKey() + "\ts.e.\t" + finalMuAndStd.getValue() + "\tC.I.\t" + (finalMuAndStd.getKey() - 1.959964 * finalMuAndStd.getValue()) + "\t" + (finalMuAndStd.getKey() + 1.959964 * finalMuAndStd.getValue());
            Tools.printResult(resultString);
            if (Tools.loggerIsOn()) {
                Tools.writeLog(resultString);
            }
            resultString = "Estimated second parameter for mask " + maskFile + "\t" + finalSecondParamAndStd.getKey() + "\ts.e.\t" + finalSecondParamAndStd.getValue() + "\tC.I.\t" + (finalSecondParamAndStd.getKey() - 1.959964 * finalSecondParamAndStd.getValue()) + "\t" + (finalSecondParamAndStd.getKey() + 1.959964 * finalSecondParamAndStd.getValue());
            Tools.printResult(resultString);
            if (Tools.loggerIsOn()) {
                Tools.writeLog(resultString);
            }
            if (applyingGeneConversionCorrection) {
                resultString = "Estimated uncorrected mu for mask " + maskFile + "\t" + finalMuAndStdUncorrected.getKey() + "\ts.e.\t" + finalMuAndStdUncorrected.getValue() + "\tC.I.\t" + (finalMuAndStdUncorrected.getKey() - 1.959964 * finalMuAndStdUncorrected.getValue()) + "\t" + (finalMuAndStdUncorrected.getKey() + 1.959964 * finalMuAndStdUncorrected.getValue());
                Tools.printResult(resultString);
                if (Tools.loggerIsOn()) {
                    Tools.writeLog(resultString);
                }
                resultString = "Estimated uncorrected mu difference for mask " + maskFile + "\t" + finalMuDiffAndStdUncorrected.getKey() + "\ts.e.\t" + finalMuDiffAndStdUncorrected.getValue() + "\tC.I.\t" + (finalMuDiffAndStdUncorrected.getKey() - 1.959964 * finalMuDiffAndStdUncorrected.getValue()) + "\t" + (finalMuDiffAndStdUncorrected.getKey() + 1.959964 * finalMuDiffAndStdUncorrected.getValue());
                Tools.printResult(resultString);
                if (Tools.loggerIsOn()) {
                    Tools.writeLog(resultString);
                }
                if (EstimateMutationRateFromIBD.inverseHeterozygosity != null) {
                    double meanDiff = finalMuDiffAndStdUncorrected.getKey();
                    double varDiff = finalMuDiffAndStdUncorrected.getValue() * finalMuDiffAndStdUncorrected.getValue();
                    double meanInverseHet = EstimateMutationRateFromIBD.inverseHeterozygosity.get(maskFile).getKey();
                    double varInverseHet = EstimateMutationRateFromIBD.inverseHeterozygosity.get(maskFile).getValue();
                    double productExpectation = meanDiff * meanInverseHet;
                    double productVar = meanDiff * meanDiff * varInverseHet + meanInverseHet * meanInverseHet * varDiff + varDiff * varInverseHet;
                    double productStdev = Math.sqrt(productVar);
                    resultString = "Estimated gene conversion rate for mask " + maskFile + "\t" + productExpectation + "\ts.e.\t" + productStdev + "\tC.I.\t" + (productExpectation - 1.959964 * productStdev) + "\t" + (productExpectation + 1.959964 * productStdev);
                    Tools.printResult(resultString);
                    if (Tools.loggerIsOn()) {
                        Tools.writeLog(resultString);
                    }
                }
            }
            // report average segment length spanning non-zero regions
            double mergedAverageLen = 0.0;
            double cnt = 0.0;
            for (int region = 0; region < resultsList.size(); region++) {
                String key = maskFile + "\t" + maxVarCountMin;
                mergedAverageLen += resultsList.get(region).get().get(key).getTotalSegmentLength();
                cnt += resultsList.get(region).get().get(key).getTotalSegmentCount();
            }
            mergedAverageLen /= cnt;
            Double[] avgLenMergedAllButValue = new Double[resultsList.size()];
            Double[] avgLenMergedAllButWeights = new Double[resultsList.size()];
            for (int exclude = 0; exclude < resultsList.size(); exclude++) {
                double AllButAverage = 0.0;
                double AllButCnt = 0.0;
                for (int region = 0; region < resultsList.size(); region++) {
                    if (region == exclude) {
                        continue;
                    }
                    String key = maskFile + "\t" + maxVarCountMin;
                    AllButAverage += resultsList.get(region).get().get(key).getTotalSegmentLength();
                    AllButCnt += resultsList.get(region).get().get(key).getTotalSegmentCount();
                }
                avgLenMergedAllButValue[exclude] = AllButAverage / AllButCnt;
                avgLenMergedAllButWeights[exclude] = AllButCnt;
            }
            Pair<Double, Double> segmentLenAverageAndStd = Tools.weightedBlockJackknife(avgLenMergedAllButValue, avgLenMergedAllButWeights, mergedAverageLen);
            String avgLenResultString = "Average IBD segment length spanning mask " + maskFile + "\t" + segmentLenAverageAndStd.getKey() + "\ts.e.\t" + segmentLenAverageAndStd.getValue() + "\tC.I.\t" + (segmentLenAverageAndStd.getKey() - 1.959964 * segmentLenAverageAndStd.getValue()) + "\t" + (segmentLenAverageAndStd.getKey() + 1.959964 * segmentLenAverageAndStd.getValue());
            Tools.printResult(avgLenResultString);
            if (Tools.loggerIsOn()) {
                Tools.writeLog(avgLenResultString);
            }
            if (useTrinucleotide) {
                if (TrinucleotideContext.context != null) {
                    Tools.warning("Processing trinucleotide results for maxVarCountMin only.");
                    Tools.printResult("Trinucleotide results");
                    Results[] muMergedAllButTrinucleotide = new Results[resultsList.size()];
                    for (int exclude = 0; exclude < resultsList.size(); exclude++) {
                        TreeMap<Integer, Results> allButOneMap = new TreeMap<Integer, Results>();
                        Results mergedAllButOne = new Results(maskFile, maxVarCountMin);
                        for (int i = 0; i < resultsList.size(); i++) {
                            if (i == exclude) {
                                continue;
                            }
                            String key = maskFile + "\t" + maxVarCountMin;
                            Results toBeMerged = resultsList.get(i).get().get(key);
                            mergedAllButOne.mergeResults(toBeMerged);
                        }
                        muMergedAllButTrinucleotide[exclude] = mergedAllButOne;
                    }
                    String[] referencelist = (trinucleotideNoDirection) ? TrinucleotideContext.trinucleotideListNoDirection : TrinucleotideContext.trinucleotideListDirectional;
                    for (String trinucleotide : referencelist) {
                        Results mergedAll = allMergedMap.get(maxVarCountMin);
                        double mergedAllResultsTrinucleotide = mergedAll.getTrinucleotideContextResults().getTrinucleotideProbability(trinucleotide);
                        Double[] mergedAllButResultsTrinucleotide = new Double[resultsList.size()];
                        Double[] mergedAllButWeightsTrinucleotide = new Double[resultsList.size()];
                        int cntNonZeroWeigts = 0;
                        for (int region = 0; region < resultsList.size(); region++) {
                            String key = maskFile + "\t" + maxVarCountMin;
                            Results allButThisRegion = muMergedAllButTrinucleotide[region];
                            mergedAllButResultsTrinucleotide[region] = allButThisRegion.getTrinucleotideContextResults().getTrinucleotideProbability(trinucleotide);
                            mergedAllButWeightsTrinucleotide[region] = allButThisRegion.getTrinucleotideContextResults().getTrinucleotideRawCount(trinucleotide);
                            if (mergedAllButWeightsTrinucleotide[region] > 0) {
                                cntNonZeroWeigts++;
                            }
                        }
                        if (cntNonZeroWeigts > 0) {
                            Pair<Double, Double> results = Tools.weightedBlockJackknife(mergedAllButResultsTrinucleotide, mergedAllButWeightsTrinucleotide, mergedAllResultsTrinucleotide);
                            Tools.printResult("\t" + trinucleotide + "\t" + mergedAllResultsTrinucleotide + "\t" + results.getKey() + "\t" + results.getValue() + "\t" + cntNonZeroWeigts);
                        } else {
                            Tools.printResult("\t" + trinucleotide + "\t" + mergedAllResultsTrinucleotide + "\tNA\tNA\t0");
                        }
                    }
                }
            }
        } else {
            Tools.warning("The code for non-jackknife computation needs testing and expansion.");
            if (useWeigthsIntMRCAReg) {
                if (tMRCARegressionWeightsMap == null) {
                    Tools.exit("Inverse-variance regression weights require jackknife.");
                } else {
                    TreeMap<Integer, Results> allMergedMap = new TreeMap<Integer, Results>();
                    for (int count : maxVarCountList) {
                        Results mergedAllThisCount = new Results(maskFile, count);
                        for (int i = 0; i < resultsList.size(); i++) {
                            String key = maskFile + "\t" + count;
                            Results toBeMerged = resultsList.get(i).get().get(key);
                            mergedAllThisCount.mergeResults(toBeMerged);
                        }
                        allMergedMap.put(count, mergedAllThisCount);
                    }
                    TreeMap<Double, BasePairsAndHets> allHist = (isUseCumulative())
                            ? allMergedMap.get(EstimateMutationRateFromIBD.largeValue).getCumulativeHistogram()
                            : allMergedMap.get(EstimateMutationRateFromIBD.largeValue).getHistogram();
                    tMRCARegressionWeights = new double[allHist.keySet().size()];
                    setWeightsFromMap(allHist);
                }
            }
            if (resultsList.get(resultsList.size() - 1).get().get(maskFile + "\t" + EstimateMutationRateFromIBD.largeValue).inMaskTot == 0) {
                Tools.printResult("Mask " + maskFile + " did not span any mutations in the mask.");
                return;
            }
            if (useMaxCountRangeRegression) {
                Tools.printProgress("Performing regression on maximum allele count frequency.");
                double avgIBDLen = 0.0;
                double totIBDseg = 0.0;
                TreeMap<Integer, Results> maxCountResList = new TreeMap<Integer, Results>();
                for (int count : maxVarCountList) {
                    for (int i = 0; i < resultsList.size(); i++) {
                        String key = maskFile + "\t" + count;
                        Results result = resultsList.get(i).get().get(key);
                        maxCountResList.put(count, result);
                        avgIBDLen = result.getTotalSegmentLength() / result.getTotalSegmentCount();
                        totIBDseg = result.getTotalSegmentCount();
                    }
                }
                double[] maxCountResultInterceptSlopeRsquare = Results.regressMuForListOfMaxCounts(-1, maxCountResList, maskFile, maxVarCountMin, maxVarCountInterval, maxVarCountMax, demography, minimumLength, maximumLength);
                Tools.printResult("Estimated mutation rate for mask " + maskFile + ". Intercept (mu)\t" + maxCountResultInterceptSlopeRsquare[0] + "\tslope (~GC)\t" + maxCountResultInterceptSlopeRsquare[1]);
                Tools.printResult("Average IBD segment length spanning mask " + maskFile + "\t" + avgIBDLen);
                Tools.printResult("Total number of IBD segments spanning mask " + maskFile + "\t" + totIBDseg);
            } else {
                for (Future<TreeMap<String, Results>> futureResult : resultsList) {
                    Results result = futureResult.get().firstEntry().getValue();
                    if (haveDemography) {
                        double[] interceptSlopeRsquare = result.getMu(demography, minimumLength, maximumLength);
                        Tools.printResult("Estimated mutation rate for mask " + maskFile + ". Intercept\t" + interceptSlopeRsquare[0] + "\tslope\t" + interceptSlopeRsquare[1]);
                        String avgLenResultString = "Average IBD segment length spanning mask " + maskFile + "\t" + result.getTotalSegmentLength() / result.getTotalSegmentCount();
                        Tools.printResult(avgLenResultString);
                    }
                }
            }
        }
    }

    public Results(String maskName) {
        this.maskName = maskName;
    }

    public Results(String maskName, int maxCount) {
        this.maskName = maskName;
        this.maxCount = maxCount;
    }

    // DATA STRUCTURES
    private double totalSegmentLength;
    private double totalSegmentCount;
    private TreeMap<Double, BasePairsAndHets> mutCountHistogram = new TreeMap<Double, BasePairsAndHets>();
    private TrinucleotideContext trinucleotideContext = new TrinucleotideContext();
    private int maxCount = EstimateMutationRateFromIBD.largeValue;
    private int inMaskTot = 0;
    private String maskName;
    private TreeMap<Double, Double> spectrumOnSegments = new TreeMap<Double, Double>();
    private TreeMap<Double, Double> transitionPerCM = new TreeMap<Double, Double>();
    private TreeMap<Double, Double> transversionPerCM = new TreeMap<Double, Double>();
    private double transitions = 0.0;
    private double transversions = 0.0;
    private TreeMap<Double, Double> countsHistogram = new TreeMap<Double, Double>();
    private TreeMap<Double, Double> inMaskHistogram = new TreeMap<Double, Double>();
    private TreeMap<Double, Integer> freqCounts = new TreeMap<Double, Integer>();
    private HashSet<String> variantsInFiles = new HashSet<String>();

    public void addTotInMask(int add) {
        inMaskTot += add;
    }

    public int getTotInMask(int add) {
        return inMaskTot;
    }

    public void increaseTrinucleotideCounts(String from, String to) {
        trinucleotideContext.increaseTrinucleotideCountsBy(from, to, 1);
    }

    public TrinucleotideContext getTrinucleotideContextResults() {
        return trinucleotideContext;
    }

    public void mergeResults(Results otherResults) throws IOException {
        if (otherResults.getMaskName().compareToIgnoreCase(getMaskName()) != 0) {
            Tools.exit("Trying to merge results from different bed files.");
        }
        if (otherResults.getMaxCount() != this.getMaxCount()) {
            Tools.exit("Trying to merge results from different max variant count values.");
        }
        transitions += otherResults.transitions;
        transversions += otherResults.transversions;
        // merge het historgram
        for (Double len : otherResults.getMutCountHistogram().keySet()) {
            BasePairsAndHets bpHets = (getMutCountHistogram().containsKey(len))
                    ? getMutCountHistogram().get(len)
                    : new BasePairsAndHets();
            bpHets.basePairs += otherResults.getMutCountHistogram().get(len).basePairs;
            bpHets.count += otherResults.getMutCountHistogram().get(len).count;
            bpHets.hets += otherResults.getMutCountHistogram().get(len).hets;
            bpHets.mutRate += otherResults.getMutCountHistogram().get(len).mutRate;
            getMutCountHistogram().put(len, bpHets);
        }
        this.trinucleotideContext.setTrinucleotideCounts(TrinucleotideContext.mergeTrinucleotideCounts(this.trinucleotideContext, otherResults.trinucleotideContext));
    }

    public void addHetCount(double length, int hets, int basePairs) {
        double roundedLength = Tools.floorToSignificantDigits(length, getSignificantDigitsInLengthForRegression());
        BasePairsAndHets basePairsAndHets = null;
        if (getMutCountHistogram().containsKey(roundedLength)) {
            basePairsAndHets = getMutCountHistogram().get(roundedLength);
        } else {
            basePairsAndHets = new BasePairsAndHets();
        }
        basePairsAndHets.basePairs += basePairs;
        basePairsAndHets.hets += hets;
        basePairsAndHets.count++;
        basePairsAndHets.mutRate += hets / (double) basePairs;
        getMutCountHistogram().put(roundedLength, basePairsAndHets);
    }

    public TreeMap<Double, BasePairsAndHets> getHistogram() {
        return getMutCountHistogram();
    }

    public TreeMap<Double, BasePairsAndHets> getCumulativeHistogram() {
        TreeMap<Double, BasePairsAndHets> cumHist = new TreeMap<Double, BasePairsAndHets>();
        int previousHets = 0;
        double previousBasePairs = 0;
        int previousCounts = 0;
        double previousMutRate = 0;
        for (Double len : getMutCountHistogram().descendingKeySet()) {
            previousHets += getMutCountHistogram().get(len).hets;
            previousBasePairs += getMutCountHistogram().get(len).basePairs;
            previousCounts += getMutCountHistogram().get(len).count;
            previousMutRate += getMutCountHistogram().get(len).mutRate;
            BasePairsAndHets thisBin = new BasePairsAndHets(previousHets, previousBasePairs, previousCounts, previousMutRate);
            cumHist.put(len, thisBin);
        }
        return cumHist;
    }

    public double getTotalBasePairsObserved() {
        double res = 0.0;
        for (BasePairsAndHets bpHet : getMutCountHistogram().values()) {
            res += bpHet.basePairs;
        }
        return res;
    }

    public String getMaskMinCountString() {
        return this.getMaskName() + "\t" + this.getMaxCount();
    }

    public double[] getMu(Demography demography, double minimumLength, double maximumLength) {
        TreeMap<Double, BasePairsAndHets> histogram = (isUseCumulative())
                ? getCumulativeHistogram()
                : getHistogram();

        if (demography == null) {
            Tools.exit("Did not initialize demography");
        }
        ArrayList<Double> obsMutVector = new ArrayList<Double>();
        ArrayList<Double> ageVector = new ArrayList<Double>();
        int cnt = 0;
        for (double len : histogram.keySet()) {
            BasePairsAndHets thisRes = histogram.get(len);
            if (thisRes.basePairs > 0 && len >= minimumLength && len <= maximumLength) { //  should not be 0, but anyway...
                double thisObsMutRate = thisRes.mutRate / (double) thisRes.count;
                if (thisRes.count > 0) {
                    obsMutVector.add(thisObsMutRate);
                    ageVector.add(2 * demography.getAgeOfSegment(len, len + Results.getLengthIncrement(), isUseCumulative()));
                    cnt++;
                } else {
                    Tools.warning("Found zero segments for " + len);
                }
            }
        }

        double[] obsMutVectorArray = new double[cnt];
        double[][] ageVectorArray = new double[cnt][1];

        for (int i = 0; i < obsMutVector.size(); i++) {
            obsMutVectorArray[i] = obsMutVector.get(i);
            ageVectorArray[i][0] = ageVector.get(i);
        }

        double[] w = new double[ageVectorArray.length];
        if (Results.getUseWeigthsIntMRCAReg()) {
            w = Results.getTMRCARegressionWeights();
        } else {
            for (int k = 0; k < w.length; k++) {
                w[k] = 1.;
            }
        }
        return Tools.simpleRegression(ageVectorArray, obsMutVectorArray, w);
    }

    public static double[] regressMuNoMaxCounts(TreeMap<Integer, Results> resultsMap, String maskFile,
            Demography demography, double minimumLength, double maximumLength) {
        // will not do regression on maxCount, only one value available
        double intercept = resultsMap.firstEntry().getValue().getMu(demography, minimumLength, maximumLength)[0];
        double mu = resultsMap.firstEntry().getValue().getMu(demography, minimumLength, maximumLength)[1];
        double rSquare = resultsMap.firstEntry().getValue().getMu(demography, minimumLength, maximumLength)[2];
        if (Tools.loggerIsOn()) {
            Tools.writeLog("\nNo GC correction for mask " + maskFile + "\t" + mu + "\t" + intercept);
        }
        return new double[]{mu, intercept, rSquare};
    }

    public static double[] regressMuForListOfMaxCounts(int region, TreeMap<Integer, Results> resultsMap, String maskFile, int countFrom, int countInterval, int countTo,
            Demography demography, double minimumLength, double maximumLength) {
        double mean = 0.;
        // will do regression on maxCount
        // size-1 because of extra bin for uncorrected value
        double[] muOfCount = new double[resultsMap.size() - 1];
        double[][] counts = new double[resultsMap.size() - 1][1];
        int i = 0;
        for (int count = countFrom; count <= countTo; count += countInterval) {
            if (LinearRegression.getPrintReg()) {
                System.out.println("Regression\tExcluding region " + region + "\tMaAF\t" + count);
            }
            Results thisCountResults = resultsMap.get(count);
            if (Tools.loggerIsOn()) {
                Tools.writeLog("\nCalling getMu for mask " + maskFile + " maxCount " + count);
            }
            muOfCount[i] = thisCountResults.getMu(demography, minimumLength, maximumLength)[1];
            mean += muOfCount[i];
            if (Tools.loggerIsOn()) {
                Tools.writeLog("\nResult: " + muOfCount[i]);
            }
            counts[i][0] = count;
            i++;
        }
        // do rgression on maxCounts
        if (Tools.loggerIsOn()) {
            Tools.writeLog("\nRegression on max counts for " + maskFile);
        }

        if (!cumulativeMaAFRegression) {
            return new double[]{mean / i / countInterval, 0, 0};
        } else {
            double[] w = new double[counts.length];
            for (int k = 0; k < w.length; k++) {
                w[k] = 1.;
            }
            if (LinearRegression.getPrintReg()) {
                System.out.println("Regression\tMaAF-Regression for all but " + region);
            }
            return Tools.simpleRegression(counts, muOfCount, w);
        }
    }

    /**
     * @return the mutCountHistogram
     */
    public TreeMap<Double, BasePairsAndHets> getMutCountHistogram() {
        return mutCountHistogram;
    }

    /**
     * @param mutCountHistogram the mutCountHistogram to set
     */
    public void setMutCountHistogram(TreeMap<Double, BasePairsAndHets> mutCountHistogram) {
        this.mutCountHistogram = mutCountHistogram;
    }

    /**
     * @return the maxCount
     */
    public int getMaxCount() {
        return maxCount;
    }

    /**
     * @param maxCount the maxCount to set
     */
    public void setMaxCount(int maxCount) {
        this.maxCount = maxCount;
    }

    /**
     * @return the maskName
     */
    public String getMaskName() {
        return maskName;
    }

    /**
     * @param maskName the maskName to set
     */
    public void setMaskName(String maskName) {
        this.maskName = maskName;
    }

    /**
     * @return the spectrumOnSegments
     */
    public TreeMap<Double, Double> getSpectrumOnSegments() {
        return spectrumOnSegments;
    }

    /**
     * @param spectrumOnSegments the spectrumOnSegments to set
     */
    public void setSpectrumOnSegments(TreeMap<Double, Double> spectrumOnSegments) {
        this.spectrumOnSegments = spectrumOnSegments;
    }

    /**
     * @return the transitionPerCM
     */
    public TreeMap<Double, Double> getTransitionPerCM() {
        return transitionPerCM;
    }

    /**
     * @param transitionPerCM the transitionPerCM to set
     */
    public void setTransitionPerCM(TreeMap<Double, Double> transitionPerCM) {
        this.transitionPerCM = transitionPerCM;
    }

    /**
     * @return the transversionPerCM
     */
    public TreeMap<Double, Double> getTransversionPerCM() {
        return transversionPerCM;
    }

    /**
     * @param transversionPerCM the transversionPerCM to set
     */
    public void setTransversionPerCM(TreeMap<Double, Double> transversionPerCM) {
        this.transversionPerCM = transversionPerCM;
    }

    /**
     * @return the transitions
     */
    public double getTransitions() {
        return transitions;
    }

    /**
     * @param transitions the transitions to set
     */
    public void setTransitions(double transitions) {
        this.transitions = transitions;
    }

    /**
     * @return the transversions
     */
    public double getTransversions() {
        return transversions;
    }

    /**
     * @param transversions the transversions to set
     */
    public void setTransversions(double transversions) {
        this.transversions = transversions;
    }

    /**
     * @return the countsHistogram
     */
    public TreeMap<Double, Double> getCountsHistogram() {
        return countsHistogram;
    }

    /**
     * @param countsHistogram the countsHistogram to set
     */
    public void setCountsHistogram(TreeMap<Double, Double> countsHistogram) {
        this.countsHistogram = countsHistogram;
    }

    /**
     * @return the inMaskHistogram
     */
    public TreeMap<Double, Double> getInMaskHistogram() {
        return inMaskHistogram;
    }

    /**
     * @param inMaskHistogram the inMaskHistogram to set
     */
    public void setInMaskHistogram(TreeMap<Double, Double> inMaskHistogram) {
        this.inMaskHistogram = inMaskHistogram;
    }

    /**
     * @return the freqCounts
     */
    public TreeMap<Double, Integer> getFreqCounts() {
        return freqCounts;
    }

    /**
     * @param freqCounts the freqCounts to set
     */
    public void setFreqCounts(TreeMap<Double, Integer> freqCounts) {
        this.freqCounts = freqCounts;
    }

    /**
     * @return the variantsInFiles
     */
    public HashSet<String> getVariantsInFiles() {
        return variantsInFiles;
    }

    /**
     * @param variantsInFiles the variantsInFiles to set
     */
    public void setVariantsInFiles(HashSet<String> variantsInFiles) {
        this.variantsInFiles = variantsInFiles;
    }

    /**
     * @return the totalSegmentLength
     */
    public double getTotalSegmentLength() {
        return totalSegmentLength;
    }

    /**
     * @param totalSegmentLength the totalSegmentLength to set
     */
    public void setTotalSegmentLength(double totalSegmentLength) {
        this.totalSegmentLength = totalSegmentLength;
    }

    /**
     * @return the totalSegmentCount
     */
    public double getTotalSegmentCount() {
        return totalSegmentCount;
    }

    /**
     * @param totalSegmentCount the totalSegmentCount to set
     */
    public void setTotalSegmentCount(double totalSegmentCount) {
        this.totalSegmentCount = totalSegmentCount;
    }

}

class BasePairsAndHets {

    double basePairs = 0;
    int hets = 0;

    int count = 0;
    double mutRate = 0.0;

    public BasePairsAndHets() {
        this.hets = 0;
        this.basePairs = 0;
        this.mutRate = 0;
        this.count = 0;
    }

    public BasePairsAndHets(int hets, double basePairs, int count, double mutRate) {
        this.hets = hets;
        this.basePairs = basePairs;
        this.mutRate = mutRate;
        this.count = count;
    }

}

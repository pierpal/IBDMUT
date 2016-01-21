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

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import org.apache.commons.math3.stat.regression.RegressionResults;
import org.apache.commons.math3.stat.regression.SimpleRegression;

/**
 *
 * @author Pier Palamara
 */
public class Tools {

    private static boolean quiet = false;
    private static boolean loggerIsOn = false;
    private static PrintWriter logFileStream;
    private static boolean useApacheRegression = false;
    private static boolean printJackknife = false;

    private static int jackknifeWarningCount = 10;

    public static void setPrintJackknife(boolean value) {
        printJackknife = value;
    }

    public static void setUseApacheReg() {
        useApacheRegression = true;
    }

    public static void makeQuiet() {
        quiet = true;
    }

    public static void initLog(String logFile) throws FileNotFoundException, UnsupportedEncodingException {
        setLogFileStream(new PrintWriter(logFile, "UTF-8"));
    }

    public static void finalizeLogIfOpen() throws FileNotFoundException, UnsupportedEncodingException {
        if (getLogFileStream() != null) {
            getLogFileStream().close();
        }
    }

    public static void writeLog(String s) {
        getLogFileStream().println(s);
    }

    public static double[] simpleRegression(double[][] x, double[] y, double[] w) {
        if (useApacheRegression) {
            Tools.exit("Apache Regression should be removed.");
            return Tools.ApacheRegression(x, y);
        } else {
            double[][] transp = new double[x[0].length + 1][x.length];
            for (int i = 0; i < x.length; i++) {
                transp[0][i] = 1.;
            }
            for (int i = 0; i < x.length; i++) {
//                w[i] = 1.   ;
                for (int j = 0; j < x[0].length; j++) {
                    transp[j + 1][i] = x[i][j];
                }
            }
            LinearRegression linReg = new LinearRegression();
            if (linReg.regress(y, transp, w)) {
                double res[] = new double[3];
                res[0] = linReg.getCoefficients()[0];
                res[1] = linReg.getCoefficients()[1];
                res[2] = linReg.getCorrelationCoefficient();
                return res;
            } else {
                Tools.warning("Linear regression failed due to too few data points or non-invertible matrix");
                return new double[]{Double.NaN, Double.NaN, Double.NaN};
            }
        }
    }

    public static double[] ApacheRegression(double[][] x, double[] y) {
        if (x.length < 2) {
//            Tools.warning("******************************************************");
//            Tools.warning("******************************************************");
//            Tools.warning("******************************************************");
//            Tools.warning("Trying to run regression with " + x.length + " points.");
//            Tools.warning("******************************************************");
//            Tools.warning("******************************************************");
//            Tools.warning("******************************************************");
            exit("Trying to run regression with " + x.length + " points.");
        }
        if (Tools.loggerIsOn()) {
            Tools.writeLog("Regression");
            Tools.writeLog("\tx\ty");
            for (int i = 0; i < x.length; i++) {
                Tools.writeLog("\t" + x[i][0] + "\t" + y[i]);
            }
        }
        SimpleRegression reg = new SimpleRegression(true);
        reg.addObservations(x, y);
        RegressionResults regRes = reg.regress();
        double[] regResValues = regRes.getParameterEstimates();
        double intercept = regResValues[0];
        double slope = regResValues[1];
        return new double[]{intercept, slope, regRes.getRSquared()};
    }

    public static Pair<Double, Double> weightedBlockJackknife(Double[] originalEstimates, Double[] originalWeights, Double estTot) {

        if (originalEstimates.length != originalWeights.length) {
            exit("Jacknife weights are not same length as estimates vector.");
        }

        if (printJackknife) {
            System.out.println("Jackknife estimateAll:\t" + estTot);
            for (int k = 0; k < originalWeights.length; k++) {
                System.out.println("Jackknife data:\t" + originalEstimates[k] + "\t" + originalWeights[k]);
            }
        }

        // get rid of zero weights
        int cnt = 0;
        for (int i = 0; i < originalWeights.length; i++) {
            if (originalWeights[i] < 0) {
                exit("Jacknife weights must be >= 0.");
            }
            if (originalWeights[i] > 0) {
                cnt++;
            }
        }

        if (cnt < jackknifeWarningCount) {
            warning("Running Jackknife with " + cnt + " data points.");
        }

        Double[] estimates = new Double[cnt];
        Double[] weights = new Double[cnt];
        int curr = 0;
        for (int i = 0; i < originalWeights.length; i++) {
            if (originalWeights[i] > 0) {
                estimates[curr] = originalEstimates[i];
                weights[curr] = originalWeights[i];
                curr++;
            }
        }

        if (estimates.length < 2) {
            exit("Trying to run jackknife with " + estimates.length + " replicate.");
        }

        int numEstimates = weights.length;
        double totWeight = 0.;
        for (int i = 0; i < weights.length; i++) {
            totWeight += weights[i];
        }

        double[] inverseContribution = new double[weights.length];
        for (int i = 0; i < inverseContribution.length; i++) {
            inverseContribution[i] = totWeight / weights[i];
        }
        double jackKnifeEstimate = numEstimates * estTot;
        for (int i = 0; i < numEstimates; i++) {
            jackKnifeEstimate -= (1 - weights[i] / totWeight) * estimates[i];
        }
        double jackKnifeVariance = 0;
        for (int i = 0; i < numEstimates; i++) {
            double currPseudo = inverseContribution[i] * estTot - (inverseContribution[i] - 1) * estimates[i] - jackKnifeEstimate;
            jackKnifeVariance += 1 / (inverseContribution[i] - 1) * currPseudo * currPseudo;
        }
        jackKnifeVariance /= (double) numEstimates;
        double se = Math.sqrt(jackKnifeVariance);
        if (Tools.loggerIsOn()) {
            Tools.writeLog("JackKnife");
            Tools.writeLog("\tglobal estimate:\t" + estTot);
            Tools.writeLog("\tweight\testimate");
            for (int i = 0; i < weights.length; i++) {
                Tools.writeLog("\t" + weights[i] + "\t" + estimates[i]);
            }
            Tools.writeLog("\tJackKnife estimate:\t" + jackKnifeEstimate + "\ts.e.\t" + se);
            Tools.writeLog(Red);
        }
        if (printJackknife) {
            System.out.println("Jackknife result:\t" + jackKnifeEstimate + "\t" + se);
        }
        return new Pair<Double, Double>(jackKnifeEstimate, se);
    }

    public static double floorToSignificantDigits(double value, int places) {
//        System.out.println(value + "\t" + places + "\t");
        return Math.floor(value * Math.pow(10., places)) / Math.pow(10., places);
    }

    public static void warning(String warning) {
        if (!quiet) {
            System.err.println(getLightRed() + "Warning:\t" + warning + getResetColor());
        }
    }

    public static void printProgress(String message) {
        if (!quiet) {
            System.err.println(getLightGreen() + message + getResetColor());
        }
    }

    public static void printVerboseProgressLevel1(String message) {
        if (!quiet) {
            System.err.println("\t" + getGreen() + message + getResetColor());
        }
    }

    public static void printVerboseProgressLevel2(String message) {
        if (!quiet) {
            System.err.println("\t" + getLightCyan() + message + getResetColor());
        }
    }

    public static void printResult(String message) {
        System.err.println(getYellow() + message + getResetColor());
    }

    public static void printVerboseProgressSameLine(String message) {
        if (!quiet) {
            System.err.print("\r\t" + getGreen() + message + getResetColor());
        }
    }

    public static void printVerboseProgressSameLineNoTab(String message) {
        if (!quiet) {
            System.err.print("\r" + getGreen() + message + getResetColor());
        }
    }

    public static void exit(String error) {
        System.err.println(getLightRed() + "ERROR:\t" + error + getResetColor());
        System.exit(1);
    }

    public static void testColors() {
        System.err.println(getBlack() + "Black" + getResetColor());
        System.err.println(getDarkGray() + "DarkGray" + getResetColor());
        System.err.println(getBlue() + "Blue" + getResetColor());
        System.err.println(getLightBlue() + "LightBlue" + getResetColor());
        System.err.println(getGreen() + "Green" + getResetColor());
        System.err.println(getLightGreen() + "LightGreen" + getResetColor());
        System.err.println(getCyan() + "Cyan" + getResetColor());
        System.err.println(getLightCyan() + "LightCyan" + getResetColor());
        System.err.println(getRed() + "Red" + getResetColor());
        System.err.println(getLightRed() + "LightRed" + getResetColor());
        System.err.println(getPurple() + "Purple" + getResetColor());
        System.err.println(getLightPurple() + "LightPurple" + getResetColor());
        System.err.println(getBrown() + "Brown" + getResetColor());
        System.err.println(getYellow() + "Yellow" + getResetColor());
        System.err.println(getLightGray() + "LightGray" + getResetColor());
        System.err.println(getWhite() + "White" + getResetColor());
    }

    private static String Black = "\033[0;30m";
    private static String DarkGray = "\033[1;30m";
    private static String Blue = "\033[0;34m";
    private static String LightBlue = "\033[1;34m";
    private static String Green = "\033[0;32m";
    private static String LightGreen = "\033[1;32m";
    private static String Cyan = "\033[0;36m";
    private static String LightCyan = "\033[1;36m";
    private static String Red = "\033[0;31m";
    private static String LightRed = "\033[1;31m";
    private static String Purple = "\033[0;35m";
    private static String LightPurple = "\033[1;35m";
    private static String Brown = "\033[0;33m";
    private static String Yellow = "\033[1;33m";
    private static String LightGray = "\033[0;37m";
    private static String White = "\033[1;37m";
    private static String ResetColor = "\u001B[0m";

    private static int waveCounter = 0;
    private static final String[] waves = new String[]{
        "▁ ▂ ▃ ▄ ▅ ▆ ▇ █ ▇ ▆ ▅ ▄ ▃ ▂ ",
        "▂ ▃ ▄ ▅ ▆ ▇ █ ▇ ▆ ▅ ▄ ▃ ▂ ▁ ",
        "▃ ▄ ▅ ▆ ▇ █ ▇ ▆ ▅ ▄ ▃ ▂ ▁ ▂ ",
        "▄ ▅ ▆ ▇ █ ▇ ▆ ▅ ▄ ▃ ▂ ▁ ▂ ▃ ",
        "▅ ▆ ▇ █ ▇ ▆ ▅ ▄ ▃ ▂ ▁ ▂ ▃ ▄ ",
        "▆ ▇ █ ▇ ▆ ▅ ▄ ▃ ▂ ▁ ▂ ▃ ▄ ▅ ",
        "▇ █ ▇ ▆ ▅ ▄ ▃ ▂ ▁ ▂ ▃ ▄ ▅ ▆ ",
        "█ ▇ ▆ ▅ ▄ ▃ ▂ ▁ ▂ ▃ ▄ ▅ ▆ ▇ ",
        "▇ ▆ ▅ ▄ ▃ ▂ ▁ ▂ ▃ ▄ ▅ ▆ ▇ █ ",
        "▆ ▅ ▄ ▃ ▂ ▁ ▂ ▃ ▄ ▅ ▆ ▇ █ ▇ ",
        "▅ ▄ ▃ ▂ ▁ ▂ ▃ ▄ ▅ ▆ ▇ █ ▇ ▆ ",
        "▄ ▃ ▂ ▁ ▂ ▃ ▄ ▅ ▆ ▇ █ ▇ ▆ ▅ ",
        "▃ ▂ ▁ ▂ ▃ ▄ ▅ ▆ ▇ █ ▇ ▆ ▅ ▄ ",
        "▂ ▁ ▂ ▃ ▄ ▅ ▆ ▇ █ ▇ ▆ ▅ ▄ ▃ "
    };

    private static final String[] wavesNoSpace = new String[]{
        "▁▂▃▄▅▆▇█▇▆▅▄▃▂",
        "▂▃▄▅▆▇█▇▆▅▄▃▂▁",
        "▃▄▅▆▇█▇▆▅▄▃▂▁▂",
        "▄▅▆▇█▇▆▅▄▃▂▁▂▃",
        "▅▆▇█▇▆▅▄▃▂▁▂▃▄",
        "▆▇█▇▆▅▄▃▂▁▂▃▄▅",
        "▇█▇▆▅▄▃▂▁▂▃▄▅▆",
        "█▇▆▅▄▃▂▁▂▃▄▅▆▇",
        "▇▆▅▄▃▂▁▂▃▄▅▆▇█",
        "▆▅▄▃▂▁▂▃▄▅▆▇█▇",
        "▅▄▃▂▁▂▃▄▅▆▇█▇▆",
        "▄▃▂▁▂▃▄▅▆▇█▇▆▅",
        "▃▂▁▂▃▄▅▆▇█▇▆▅▄",
        "▂▁▂▃▄▅▆▇█▇▆▅▄▃"
    };

    //from http://www.ascii-art.de/ascii/s/stickman.txt
    private static final String[] rotatingMan = new String[]{
        "Watch ...\n\n\n o  \n/|\\ \n/ \\ ",
        "Watch ...\n\n\n     \\ o /\n       |  \n      / \\ ",
        "Watch ...\n\n\n            _ o\n             /\\\n            | \\ ",
        "Watch ...\n\n\n \n                  ___\\o\n                 /)  | ",
        "Watch ...\n\n\n                        __|  \n                          \\o \n                          ( \\ ",
        "Watch ...\n\n\n                               \\ / \n                                |  \n                               /o\\ ",
        "Watch ...\n\n\n                                      |__\n                                     o/   \n                                    / )   ",
        "Watch ...\n\n\n                                           \n                                          o/__ \n                                           |  (\\n",
        "Watch ...\n\n                                                  o _\n                                                  /\\\n                                                  / |\n",
        "Watch ...\n\n\n                                                       \\ o /\n                                                         | \n                                                        / \\ ",
        "Watch ...\n\n\n                                                               o\n                                                              /|\\\n                                                              / \\"
    };

    public static void printAnimation() {
        String[] animation;
        animation = rotatingMan;
//        animation = wavesNoSpace;
//        animation = waves;
        System.out.print("\u001b[2J");
        System.out.flush();
        Tools.printVerboseProgressSameLineNoTab(animation[waveCounter]);
        waveCounter = (waveCounter + 1) % animation.length;
    }

    /**
     * @return the logFileStream
     */
    public static PrintWriter getLogFileStream() {
        return logFileStream;
    }

    /**
     * @param aLogFileStream the logFileStream to set
     */
    public static void setLogFileStream(PrintWriter aLogFileStream) {
        logFileStream = aLogFileStream;
    }

    /**
     * @return the Black
     */
    public static String getBlack() {
        return Black;
    }

    /**
     * @return the DarkGray
     */
    public static String getDarkGray() {
        return DarkGray;
    }

    /**
     * @return the Blue
     */
    public static String getBlue() {
        return Blue;
    }

    /**
     * @return the LightBlue
     */
    public static String getLightBlue() {
        return LightBlue;
    }

    /**
     * @return the Green
     */
    public static String getGreen() {
        return Green;
    }

    /**
     * @param aGreen the Green to set
     */
    public static void setGreen(String aGreen) {
        Green = aGreen;
    }

    /**
     * @return the LightGreen
     */
    public static String getLightGreen() {
        return LightGreen;
    }

    /**
     * @param aLightGreen the LightGreen to set
     */
    public static void setLightGreen(String aLightGreen) {
        LightGreen = aLightGreen;
    }

    /**
     * @return the Cyan
     */
    public static String getCyan() {
        return Cyan;
    }

    /**
     * @param aCyan the Cyan to set
     */
    public static void setCyan(String aCyan) {
        Cyan = aCyan;
    }

    /**
     * @return the LightCyan
     */
    public static String getLightCyan() {
        return LightCyan;
    }

    /**
     * @param aLightCyan the LightCyan to set
     */
    public static void setLightCyan(String aLightCyan) {
        LightCyan = aLightCyan;
    }

    /**
     * @return the Red
     */
    public static String getRed() {
        return Red;
    }

    /**
     * @param aRed the Red to set
     */
    public static void setRed(String aRed) {
        Red = aRed;
    }

    /**
     * @return the LightRed
     */
    public static String getLightRed() {
        return LightRed;
    }

    /**
     * @param aLightRed the LightRed to set
     */
    public static void setLightRed(String aLightRed) {
        LightRed = aLightRed;
    }

    /**
     * @return the Purple
     */
    public static String getPurple() {
        return Purple;
    }

    /**
     * @param aPurple the Purple to set
     */
    public static void setPurple(String aPurple) {
        Purple = aPurple;
    }

    /**
     * @return the LightPurple
     */
    public static String getLightPurple() {
        return LightPurple;
    }

    /**
     * @param aLightPurple the LightPurple to set
     */
    public static void setLightPurple(String aLightPurple) {
        LightPurple = aLightPurple;
    }

    /**
     * @return the Brown
     */
    public static String getBrown() {
        return Brown;
    }

    /**
     * @param aBrown the Brown to set
     */
    public static void setBrown(String aBrown) {
        Brown = aBrown;
    }

    /**
     * @return the Yellow
     */
    public static String getYellow() {
        return Yellow;
    }

    /**
     * @param aYellow the Yellow to set
     */
    public static void setYellow(String aYellow) {
        Yellow = aYellow;
    }

    /**
     * @return the LightGray
     */
    public static String getLightGray() {
        return LightGray;
    }

    /**
     * @param aLightGray the LightGray to set
     */
    public static void setLightGray(String aLightGray) {
        LightGray = aLightGray;
    }

    /**
     * @return the White
     */
    public static String getWhite() {
        return White;
    }

    /**
     * @param aWhite the White to set
     */
    public static void setWhite(String aWhite) {
        White = aWhite;
    }

    /**
     * @return the ResetColor
     */
    public static String getResetColor() {
        return ResetColor;
    }

    /**
     * @param aResetColor the ResetColor to set
     */
    public static void setResetColor(String aResetColor) {
        ResetColor = aResetColor;
    }

    /**
     * @return the loggerIsOn
     */
    public static boolean loggerIsOn() {
        return loggerIsOn;
    }

    /**
     * @param aLoggerIsOn the loggerIsOn to set
     */
    public static void setLoggerIsOn(boolean aLoggerIsOn) {
        loggerIsOn = aLoggerIsOn;
    }

    /**
     * @param aBlack the Black to set
     */
    public static void setBlack(String aBlack) {
        Black = aBlack;
    }

    /**
     * @param aDarkGray the DarkGray to set
     */
    public static void setDarkGray(String aDarkGray) {
        DarkGray = aDarkGray;
    }

    /**
     * @param aBlue the Blue to set
     */
    public static void setBlue(String aBlue) {
        Blue = aBlue;
    }

    /**
     * @param aLightBlue the LightBlue to set
     */
    public static void setLightBlue(String aLightBlue) {
        LightBlue = aLightBlue;
    }

    static void printHelp() {
        String s = "IBDMUT 0.11.12.15\thttps://github.com/pierpal/IBDMUT\n" +
                "Contact: Pier Palamara, ppalama AT hsph DOT harvard DOTAGAIN edu\n\n" +
                "Options: \n"+
                "\t--demography [file] (demographic model with format \"generation\tsize\" for each line) \n" +
                "\t--plink [fileRoot] (plink file root; expects .ped.gz, .map, and .frq) \n" +
                "\t--match [file] (match.gz file in GERMLINE format) \n" +
                "\nSome of the optional or alternative flags:\n" +
                "\t--plinkList [file] (substitutes --plink and --match for many files) \n" +
                "\t--maskList [file] (substitutes --mask for many files) \n" +
                "\t--lenRange [fromLen] [toLen] (default: 2.0 5.0) \n" +
                "\t--jackknife (only if --plinkList is used with several independent regions) \n" +
                "\t--saveBin [suffix] (saves a binary file, which will load much faster than the ped.gz file) \n" +
                "\t--loadBin [suffix] (load a binary file) \n" +
                "\t--offsetCM [value] (distance to be excluded from edges; default: 0.0) \n" +
                "\t--mask [file] (bed file with regions to be included in analysis) \n" +
                "\t--threads [value] (default: 1) \n" +
                "\t--MaAFRegression [intValueFromCount] [intValueStep] [intValueToCount] (from MAF counts, interval MAF counts, to MAF counts) \n";
//                "\t--posteriorRange [doubleValue] [doubleValue] \n" +
//                "\t--computeHeterozygosity [file]";
        Tools.printResult(s);
    }

    static void printHelpAndExit(String s) {
        printHelp();
        Tools.exit(s);
    }

}

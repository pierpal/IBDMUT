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

/**
 *
 * @author Pier Palamara
 */
public class Demography {

    ArrayList<Double> timesInput = new ArrayList<Double>();
    ArrayList<Double> sizesInput = new ArrayList<Double>();
    ArrayList<Double> sizePerGenStartingAtZero = new ArrayList<Double>();

    static final int maxGen = 10000;

    public Demography(String fileName) throws IOException {
        Tools.printProgress("Reading demography from file " + fileName);
        BufferedReader br = null;
        try {
            br = new BufferedReader(new FileReader(fileName));
        } catch (FileNotFoundException ex) {
            Tools.exit("Could not open demography file " + fileName);
        }
        String line = br.readLine();
        int cnt = 0;
        while (line != null) {
            String[] strSplit = line.split("\\s+");
            if (strSplit.length != 2) {
                Tools.exit("Parsed line with " + strSplit.length + " entries in demographic file " + fileName + ". Use format\ntime\tsize");
            }
            timesInput.add(Double.parseDouble(strSplit[0]));
            sizesInput.add(Double.parseDouble(strSplit[1]));
            line = br.readLine();
        }
        if (timesInput.get(0) != 0) {
            Tools.exit("First line does not contain time = 0 in demographic file " + fileName);
        }
        setSize();
        Tools.printVerboseProgressLevel1("Demography ends at generation " + (sizePerGenStartingAtZero.size() - 1));
    }

    // TODO: this is slow. Could use closed forms for intervals
    private void setSize() {
        for (int i = 0; i < timesInput.size() - 1; i++) {
            double timeFrom = timesInput.get(i);
            double timeTo = timesInput.get(i + 1);
            double sizeCurrent = sizesInput.get(i);
            double sizeAncestral = sizesInput.get(i + 1);
            for (int gen = 0; gen < timeTo - timeFrom; gen++) {
                double timeInterval = timeTo - timeFrom;
                double thisSize = sizeCurrent * Math.pow((sizeAncestral / sizeCurrent), gen / timeInterval);
                sizePerGenStartingAtZero.add(thisSize);
            }
        }
        sizePerGenStartingAtZero.add(sizesInput.get(sizesInput.size() - 1));
    }

    public double getAgeOfSegment(double fromLen, double toLen, boolean upperInfinite) {
        fromLen /= 100.; //in cM
        toLen /= 100.; //in cM
        double pNotCoal = 1;
        ArrayList<Double> segCoalPerGenStartingAtZero = new ArrayList<Double>();
        segCoalPerGenStartingAtZero.add(0.); // can't coal at gen 0;
        double tot = 0;
        double ancestral = sizePerGenStartingAtZero.get(sizePerGenStartingAtZero.size() - 1);
        for (int gen = 1; gen < maxGen; gen++) {
            double popSize = gen < sizePerGenStartingAtZero.size() ? sizePerGenStartingAtZero.get(gen) : ancestral;
            double pCoal = 1 / popSize;
            double thisSegCoal;
            if (fromLen == toLen) {
                thisSegCoal = (4 * gen * gen) * Math.exp(-2 * gen * toLen); // exact length
            } else {
                thisSegCoal = (upperInfinite)
                        ? pCoal * pNotCoal * 2 * Math.exp(-2 * gen * fromLen) * gen
                        : pCoal * pNotCoal * ((2 * Math.exp(-2 * gen * toLen) * gen) - (2 * Math.exp(-2 * gen * fromLen) * gen));
            }
            segCoalPerGenStartingAtZero.add(thisSegCoal);
            pNotCoal = pNotCoal * (1 - pCoal);
            tot += thisSegCoal;
        }
        double expAge = 0.0;
        for (int gen = 1; gen < maxGen; gen++) {
            expAge += gen * segCoalPerGenStartingAtZero.get(gen) / tot;
        }
        return expAge;
    }
}

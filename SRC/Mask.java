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
import java.util.TreeMap;

/**
 *
 * @author Pier Palamara
 */
public class Mask {

    private TreeMap<Integer, Integer> mask = new TreeMap<Integer, Integer>();

    private TreeMap<String, Pair<Double, Double>> heterozygosityNumeratorDenominatorForEachRegion = new TreeMap<String, Pair<Double, Double>>();

    private void setHeterozygosityNumeratorDenominator(String region, double num, double den) {
        heterozygosityNumeratorDenominatorForEachRegion.put(region, new Pair<Double, Double>(num, den));
    }

    public Mask(String maskFile, String chr) {
        TreeMap<Integer, Integer> myMask = new TreeMap<Integer, Integer>();
        BufferedReader br = null;
        try {
            br = new BufferedReader(new FileReader(maskFile));
        } catch (FileNotFoundException ex) {
            Tools.exit("Could not open mask file " + maskFile);
        }
        try {
            String line = br.readLine();
            double cnt = 0;
            while (line != null) {
                String[] strSplit = line.split("\\s+");
                if (strSplit.length == 0) {
                    continue;
                }
                String chrom = strSplit[0];
                if (chrom.compareToIgnoreCase(chr) != 0) {
                    line = br.readLine();
                    continue;
                }
                int from = Integer.parseInt(strSplit[1]);
                int to = Integer.parseInt(strSplit[2]);
                myMask.put(from, to);
                line = br.readLine();
                cnt++;
            }
            Tools.printVerboseProgressLevel1("Read " + cnt + " mask lines from " + maskFile);
        } catch (IOException ex) {
            Tools.exit("Could not read mask file " + maskFile);
        }
        mask = myMask;
    }

    public double getSize(int from, int to) {
        double size = 0;
        for (int i : mask.keySet()) {
            if (mask.get(i) < from || i > to) {
                continue;
            }
            size += Math.min(to, mask.get(i)) - Math.max(from, i);
        }
        return size;
    }

    public boolean contains(int pos) {
        java.util.Map.Entry<Integer, Integer> floor = getMask().floorEntry(pos);
        if (floor == null) {
            return false;
        }
        return (pos >= floor.getKey() && pos <= floor.getValue());
    }

    public int getMaskOverlapWithRegion(int fromRegion, int toRegion) {
        int overlap = 0;
        java.util.Map.Entry<Integer, Integer> currentRegion = getMask().floorEntry(fromRegion);
        if (currentRegion == null) {
            // segment starts before first region in map
            currentRegion = getMask().ceilingEntry(fromRegion);
            if (currentRegion == null) {
                // map is empty, should throw an error maybe
                return 0;
            }
        }
        while (currentRegion != null && currentRegion.getKey() <= toRegion) {
            overlap += getOverlapBetweenTwoSegments(fromRegion, toRegion, currentRegion.getKey(), currentRegion.getValue());
            currentRegion = getMask().higherEntry(currentRegion.getKey());
            if (currentRegion != null) {
            }
        }
        return overlap;
    }

    public int getOverlapBetweenTwoSegments(int from1, int to1, int from2, int to2) {
        int overlap = Math.min(to1, to2) - Math.max(from1, from2);
        return (overlap > 0) ? overlap : 0;
    }

    /**
     * @return the mask
     */
    public TreeMap<Integer, Integer> getMask() {
        return mask;
    }

    /**
     * @param mask the mask to set
     */
    public void setMask(TreeMap<Integer, Integer> mask) {
        this.mask = mask;
    }

}

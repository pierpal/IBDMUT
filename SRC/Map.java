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

import java.util.ArrayList;
import java.util.TreeMap;

/**
 *
 * @author Pier Palamara
 */
class Map {

    private String chr;
    private ArrayList<String> ID;
    private TreeMap<Double, Integer> genPos;
    private ArrayList<Double> genPosArray;
    private ArrayList<Integer> physPos;
    private int size;

    public Map(String chr, ArrayList<String> ID, TreeMap<Double, Integer> genPos, ArrayList<Integer> physPos) {
        this.chr = chr;
        this.ID = ID;
        this.genPos = genPos;
        this.physPos = physPos;
        this.size = physPos.size();
        genPosArray = new ArrayList<Double>(genPos.keySet());
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

    /**
     * @return the ID
     */
    public ArrayList<String> getID() {
        return ID;
    }

    /**
     * @param ID the ID to set
     */
    public void setID(ArrayList<String> ID) {
        this.ID = ID;
    }

    /**
     * @return the genPos
     */
    public TreeMap<Double, Integer> getGenPos() {
        return genPos;
    }

    /**
     * @param genPos the genPos to set
     */
    public void setGenPos(TreeMap<Double, Integer> genPos) {
        this.genPos = genPos;
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
     * @return the size
     */
    public int getSize() {
        return size;
    }

    /**
     * @param size the size to set
     */
    public void setSize(int size) {
        this.size = size;
    }
}

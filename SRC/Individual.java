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
import java.util.BitSet;

/**
 *
 * @author Pier Palamara
 */
public class Individual {

    private BitSet bitSeq;
    private BitSet mask;
    
    private ArrayList<Float> posteriors = new ArrayList<Float>();
    
    private String famId;
    private String indId;
    private String patId;
    private String matId;
    private String sex;
    private String pheno;

    public Individual(String famId, String indId, String patId, String matId, String sex, String pheno, BitSet bitSeq, BitSet mask) {
        this.famId = famId;
        this.indId = indId;
        this.matId = matId;
        this.patId = patId;
        this.sex = sex;
        this.pheno = pheno;
        this.bitSeq = bitSeq;
        this.mask = mask;
    }

    public void printFromBitSet(ArrayList<Variant> variants) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < getBitSeq().size(); i++) {
            sb.append(' ');
            Variant currVar = variants.get(i);
            char thisVar = (getBitSeq().get(i)) ? currVar.getVarOne() : currVar.getVarZero();
            sb.append(thisVar);
        }
        System.out.println(getFam() + sb);
    }

    String getFam() {
        return this.getFamId() + "\t" + this.getIndId() + "\t" + this.getMatId() + "\t" + this.getPatId() + "\t" + this.getSex() + "\t" + this.getPheno();
    }

    /**
     * @return the bitSeq
     */
    public BitSet getBitSeq() {
        return bitSeq;
    }

    /**
     * @param bitSeq the bitSeq to set
     */
    public void setBitSeq(BitSet bitSeq) {
        this.bitSeq = bitSeq;
    }

    /**
     * @return the mask
     */
    public BitSet getMask() {
        return mask;
    }

    /**
     * @param mask the mask to set
     */
    public void setMask(BitSet mask) {
        this.mask = mask;
    }

    /**
     * @return the posteriors
     */
    public ArrayList<Float> getPosteriors() {
        return posteriors;
    }

    /**
     * @param posteriors the posteriors to set
     */
    public void setPosteriors(ArrayList<Float> posteriors) {
        this.posteriors = posteriors;
    }

    /**
     * @return the famId
     */
    public String getFamId() {
        return famId;
    }

    /**
     * @param famId the famId to set
     */
    public void setFamId(String famId) {
        this.famId = famId;
    }

    /**
     * @return the indId
     */
    public String getIndId() {
        return indId;
    }

    /**
     * @param indId the indId to set
     */
    public void setIndId(String indId) {
        this.indId = indId;
    }

    /**
     * @return the patId
     */
    public String getPatId() {
        return patId;
    }

    /**
     * @param patId the patId to set
     */
    public void setPatId(String patId) {
        this.patId = patId;
    }

    /**
     * @return the matId
     */
    public String getMatId() {
        return matId;
    }

    /**
     * @param matId the matId to set
     */
    public void setMatId(String matId) {
        this.matId = matId;
    }

    /**
     * @return the sex
     */
    public String getSex() {
        return sex;
    }

    /**
     * @param sex the sex to set
     */
    public void setSex(String sex) {
        this.sex = sex;
    }

    /**
     * @return the pheno
     */
    public String getPheno() {
        return pheno;
    }

    /**
     * @param pheno the pheno to set
     */
    public void setPheno(String pheno) {
        this.pheno = pheno;
    }
}

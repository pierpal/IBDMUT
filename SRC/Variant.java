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

/**
 *
 * @author Pier Palamara
 */

class Variant {

    private char[] vars = new char[2];

    public char getVarZero() {
        return getVars()[0];
    }

    public char getVarOne() {
        return getVars()[1];
    }
    
    public Variant() {
        this.vars[0] = '0';
        this.vars[1] = '0';
    }

    public Variant(char var1) {
        this.vars[0] = var1;
        this.vars[1] = '0';
    }

    public Variant(char var1, char var2) {
        this.vars[0] = var1;
        this.vars[1] = var2;
    }

    public boolean contains(char var) {
        return getVars()[0] == var || getVars()[1] == var;
    }

    public boolean alleleToBit(char var) throws Exception {
        if (this.contains(var)) {
            if (getVars()[0] == var) {
                return false;
            } else if (getVars()[1] == var) {
                return true;
            } else {
                throw new Exception("Requested a site that was not previously recorded: " + var);
            }
        } else {
            add(var); // add it first
            return alleleToBit(var); // call again. Will be there now
        }
    }

    public void add(char var) throws Exception {
        if (getVars()[0] == '0') {
            this.getVars()[0] = var;
        } else if (getVars()[1] == '0') {
            this.getVars()[1] = var;
        } else {
            throw new Exception("triallelic site? Adding " + var + ", have already seen " + getVars()[0] + " and " + getVars()[1]);
        }
    }

    public boolean isTransition() {
        return getVars()[0] == 'A' && getVars()[1] == 'G' || getVars()[0] == 'G' && getVars()[1] == 'A' || getVars()[0] == 'C' && getVars()[1] == 'T' || getVars()[0] == 'T' && getVars()[1] == 'C';
    }

    /**
     * @return the vars
     */
    public char[] getVars() {
        return vars;
    }

    /**
     * @param vars the vars to set
     */
    public void setVars(char[] vars) {
        this.vars = vars;
    }

}

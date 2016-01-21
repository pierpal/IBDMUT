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
 * @author Pier Palamara based on code by Sindharta Tanuwijaya downloaded from
 * http://sin-logs.blogspot.com/2009/04/weighted-linear-regression-in-java-and.html
 * which is a java implementation of C# version described by Walt Fair, Jr. in
 * http://www.codeproject.com/Articles/25335/An-Algorithm-for-Weighted-Linear-Regression
 * which cites Draper, N. R. and H. Smith, Applied Regression Analysis, New
 * York: Wiley (1966)
 */

/*
 * LinearRegression.java
 *

 * USAGE:
 *
 * OPTIONS:
 *
 * CONFIGURATION AND ENVIRONMENT:
 *
 * DEPENDENCIES:
 *
 * INCOMPATIBILITIES:
 * BUGS AND LIMITATIONS:
 *
 * LICENCE AND COPYRIGHT :
 * Walt Fair, Jr.
 * http://www.codeproject.com/KB/recipes/LinReg.aspx
 * Ported from LinearRegression.cs to Java by Sindharta Tanuwijaya, 2009
 * http://sin-memories.blogspot.com/
 *
 * EXAMPLES:
 * FAQ:
 *
 * COMMON USAGE MISTAKES:
 * SEE ALSO:
 *
 * ACKNOWLEDGEMENTS:
 */
public class LinearRegression {

    double[][] V;            // Least squares and var/covar matrix
    public double[] C;      // Coefficients
    public double[] SEC;    // Std Error of coefficients
    double RYSQ;            // Multiple correlation coefficient
    double SDV;             // Standard deviation of errors
    double FReg;            // Fisher F statistic for regression
    double[] Ycalc;         // Calculated values of Y
    double[] DY;            // Residual values of Y

    private static boolean printReg = false;

    public static void setPrintReg(boolean value) {
        printReg = value;
    }

    public static boolean getPrintReg() {
        return printReg;
    }

    /*
     * @return Fisher F statistic for regression
     */
    public double getFisherF() {
        return FReg;
    }

    /*@return Multiple Correlation Coefficient*/
    public double getCorrelationCoefficient() {
        return RYSQ;
    }


    /*@return standard deviation of errors*/
    public double getStandardDeviation() {
        return SDV;
    }

    /*@return calculated values of Y*/
    public double[] getCalculatedValues() {
        return Ycalc;
    }

    /*@return residual values of Y*/
    public double[] getResiduals() {
        return DY;
    }

    /*@return coefficients*/
    public double[] getCoefficients() {
        return C;
    }

    /*@return standard error of coefficients*/
    public double[] getCoefficientsStandardError() {
        return SEC;
    }

    /*@return variance matrix*/
    public double[][] VarianceMatrix() {
        return V;
    }

    /*Perform the regression
     * @param array of Y (dependent) values
     * @param array of array of X (independent) values, e.g Y = a + bx + cx^2 + ..
     * @param array of weights of each observation/input
     * @return whether regression works or not
     */
    public boolean regress(double[] Y, double[][] X, double[] W) {
        if (X.length < 0) {
            return false;
        }

        if (printReg) {
//            System.out.println(X[0].length + "\t" + Y.length + "\t" + W.length);
            for (int j = 0; j < X[1].length; j++) {
                System.out.println("Regression\t" + j + "\t" + X[1][j] + "\t" + Y[j] + "\t" + W[j]);
            }
        }
        int M = Y.length;                   // M = Number of data points
        int N = X.length * X[0].length / M; // N = Number of linear terms
        int NDF = M - N;                    // Degrees of freedom
        Ycalc = new double[M];
        DY = new double[M];
        // If not enough data, don't attempt regression
        if (NDF < 1) {
            System.out.println("Regression\tnot enough degrees of freedom");
            return false;
        }
        V = new double[N][N];
        C = new double[N];
        SEC = new double[N];
        double[] B = new double[N];   // Vector for LSQ

        // Clear the matrices to start out
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                V[i][j] = 0;
            }
        }

        // Form Least Squares Matrix
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                V[i][j] = 0;
                for (int k = 0; k < M; k++) {
                    V[i][j] = V[i][j] + W[k] * X[i][k] * X[j][k];
                }
            }
            B[i] = 0;
            for (int k = 0; k < M; k++) {
                B[i] = B[i] + W[k] * X[i][k] * Y[k];
            }
        }
        // V now contains the raw least squares matrix
        if (!invertSymmetricMatrix(V)) {
            System.out.println("Regression\tcannot invert matrix");
            return false;
        }

        // V now contains the inverted least square matrix
        // Matrix multpily to get coefficients C = VB
        for (int i = 0; i < N; i++) {
            C[i] = 0;
            for (int j = 0; j < N; j++) {
                C[i] = C[i] + V[i][j] * B[j];
            }
        }

        // Calculate statistics
        double TSS = 0;
        double RSS = 0;
        double YBAR = 0;
        double WSUM = 0;
        for (int k = 0; k < M; k++) {
            YBAR = YBAR + W[k] * Y[k];
            WSUM = WSUM + W[k];
        }
        YBAR = YBAR / WSUM;
        for (int k = 0; k < M; k++) {
            Ycalc[k] = 0;
            for (int i = 0; i < N; i++) {
                Ycalc[k] = Ycalc[k] + C[i] * X[i][k];
            }
            DY[k] = Ycalc[k] - Y[k];
            TSS = TSS + W[k] * (Y[k] - YBAR) * (Y[k] - YBAR);
            RSS = RSS + W[k] * DY[k] * DY[k];
        }
        double SSQ = RSS / NDF;
        RYSQ = 1 - RSS / TSS;
        FReg = 9999999;
        if (RYSQ < 0.9999999) {
            FReg = RYSQ / (1 - RYSQ) * NDF / (N - 1);
        }
        SDV = Math.sqrt(SSQ);

        // Calculate var-covar matrix and std error of coefficients
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                V[i][j] = V[i][j] * SSQ;
            }
            SEC[i] = Math.sqrt(V[i][ i]);
        }
        if (printReg) {
            System.out.println("Regression\tresults\t" + getCoefficients()[0] + "\t" + getCoefficients()[1] + "\t" + getCorrelationCoefficient());
        }
        return true;
    }


    /*Invert symmetric matrix*/
    public boolean invertSymmetricMatrix(double[][] V) {
        if (V.length < 1) {
            return false;
        }

        int N = (int) Math.sqrt(V.length * V[0].length);
        double[] t = new double[N];
        double[] Q = new double[N];
        double[] R = new double[N];
        double AB;
        int K, L, M;

        // Invert a symetric matrix in V
        for (M = 0; M < N; M++) {
            R[M] = 1;
        }
        K = 0;
        for (M = 0; M < N; M++) {
            double Big = 0;
            for (L = 0; L < N; L++) {
                AB = Math.abs(V[L][L]);
                if ((AB > Big) && (R[L] != 0)) {
                    Big = AB;
                    K = L;
                }
            }
            if (Big == 0) {
                return false;
            }

            R[K] = 0;
            Q[K] = 1 / V[K][K];
            t[K] = 1;
            V[K][ K] = 0;
            if (K != 0) {
                for (L = 0; L < K; L++) {
                    t[L] = V[L][ K];
                    if (R[L] == 0) {
                        Q[L] = V[L][K] * Q[K];
                    } else {
                        Q[L] = -V[L][K] * Q[K];
                    }
                    V[L][ K] = 0;
                }
            }

            if ((K + 1) < N) {
                for (L = K + 1; L < N; L++) {
                    if (R[L] != 0) {
                        t[L] = V[K][L];
                    } else {
                        t[L] = -V[K][L];
                    }
                    Q[L] = -V[K][L] * Q[K];
                    V[K][L] = 0;
                }
            }

            for (L = 0; L < N; L++) {
                for (K = L; K < N; K++) {
                    V[L][ K] = V[L][K] + t[L] * Q[K];
                }
            }

        }
        M = N;
        L = N - 1;
        for (K = 1; K < N; K++) {
            M = M - 1;
            L = L - 1;
            for (int J = 0; J <= L; J++) {
                V[M][J] = V[J][ M];
            }
        }
        return true;
    }

}

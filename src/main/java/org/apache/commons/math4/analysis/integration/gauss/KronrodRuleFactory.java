package org.apache.commons.math4.analysis.integration.gauss;

import org.apache.commons.math4.exception.DimensionMismatchException;
import org.apache.commons.math4.linear.EigenDecomposition;
import org.apache.commons.math4.linear.RealMatrix;
import org.apache.commons.math4.util.Pair;

import java.util.Arrays;

/**
 * Generates the (2n+1)-point Gauss-Kronrod quadrature rule.
 * TODO: Finish docs
 */
public class KronrodRuleFactory extends BaseRuleFactory<Double> {

    @Override
    protected Pair<Double[], Double[]> computeRule(int numberOfPoints) throws DimensionMismatchException {
        double[] v = Jacobi(numberOfPoints*2);
        double[][] v2 = JacobiKronrod(numberOfPoints, v);
        return Kronrod(2, v2);
    }

    /**
     * Generates n recurrence coefficients for monic Jacobi polynomials. The n alpha-coefficients
     * are stored in the first column, the n beta-coefficients in
     * the second column, of the nx2 array ab.
     * Orthogonal on [0,1], weight function w(t)=(1-t)^a t^b.
     * @param N Number of Points
     * @return
     */
    private double[] Jacobi(int N){
        if (N == 1){
            final Double[] w = new Double[1];
            final Double[] x = new Double[1];
            w[0] = 0d;
            x[0] = 2d;
            ///return new Pair<>(w, x);
            //TODO
        }
        N -= 1;
        int[] n = new int[N];

        // TODO: Merge
        for (int i = 0, j = 1; i < N; i++, j++){
            n[i] = j;
        }
        double[] nab = new double[N];
        for (int i = 0; i < N; i++) {
            nab[i] = (n[i] * 2);
        }
        n = new int[N - 1];
        for (int i = 0, j = 2; i < N-1; i++, j++){
            n[i] = j;
        }
        double[] nab2 = new double[n.length];
        int counter = 0;
        for (int i : n) {
            nab2[counter] = nab[i - 1];
            counter++;
        }
        double b1 = 1d/3d;
        double[] v1 = new double[n.length];
        for (int i = 0; i < n.length; i++){
            v1[i] = (Math.pow(n[i], 4) * 4)/(Math.pow(nab2[i], 2) * (Math.pow(nab2[i], 2) - 1));
        }
        N += 1;
        double[] out1 = new double[N];
        out1[0] = 2d;
        out1[1] = b1;
        System.arraycopy(v1, 0, out1, 2, out1.length - 2);


        double[] n2 = new double[N];
        n2[0] = 1d;
        for (int i = 1; i < N; i++){
            n2[i] = out1[i] / 4d;
        }
        return n2;
    }

    /**
     *
     * @param N
     * @param ab2
     * @return
     */
    private double[][] JacobiKronrod(int N, double[] ab2){
        if (ab2.length < ((3 * N + 2 - 1)/2) + 1){
            //throw new Exception();
            return null;
        }
        double[] a = new double[2 * N + 1];
        double[] b = new double[2 * N + 1];
        Arrays.fill(a, 0.5d);

        // b(k+1)=ab0(k+1,2);
        System.arraycopy(ab2, 0, b, 0, (((3 * N) + 1) / 2) + 1);
        double[] s = new double[N / 2 + 2];
        double[] t = new double[N / 2 + 2];
        t[1] = b[N + 1];

        /////////
        for (int i = 0; i < N - 1; i++) {
            int [] l = new int[(i + 1) / 2 + 1];
            for (int j = (i + 1) / 2, k = 0; j > -1; j--, k++) {
                l[k] = i - j;
            }

            double [] sum = new double[l.length];
            for (int j = 0, k = (i + 1) / 2; j < l.length; j++, k--) {
                sum[j] = (b[k + N + 1] * s[k]) - (b[l[j]] * s[k + 1]);
            }

            sum = CumulativeSum(sum);
            for (int j = 0, k = (i + 1) / 2; j < sum.length; j++, k--) {
                s[k + 1] = sum[j];
            }

            double[] swap = s;
            s = t;
            t = swap;
        }

        System.arraycopy(s, -1 + 1, s, 1, N / 2 + 1);
        for (int i = N - 1; i < ((2 * N) - 3) + 1; i++) {
            int [] l = new int[((i - 1) / 2) - (i + 1 - N) + 1];

            for (int j = (i + 1) - N, k = 0; k < l.length; j++, k++) {
                l[k] = i - j;
            }
            double [] sum = new double[l.length];
            for (int j = 0, k = (i + 1) - N; j < l.length; j++, k++) {
                sum[j] = -b[k + N + 1] * s[N - l[j]] + b[l[j]] * s[N - l[j] + 1];
            }
            sum = CumulativeSum(sum);
            for (int j = 0; j < sum.length; j++) {
                s[N - l[j]] = sum[j];
            }

            if (i % 2 != 0){
                int k = (i + 1) / 2;
                b[k + N + 1] = s[N - l[l.length - 1]] / s[N + 1 - l[l.length - 1]];
            }
            double[] swap = s;
            s = t;
            t = swap;
        }
        return new double[][]{a, b};
    }

    /**
     *
     * @param N
     * @param ab0
     * @return
     */
    private Pair<Double[], Double[]> Kronrod(int N, double[][] ab0){
        double[][] j = new double[2 * N + 1][2 * N + 1];
        double[] diag1 = new double[2 * N];
        double[] diag2 = new double[2 * N + 1];
        Arrays.fill(diag2, 0.5d);

        for (int i = 0; i < j.length - 1; i++) {
            j[i][i] = ab0[0][i];
            j[i][i + 1] = Math.sqrt(ab0[1][i + 1]);
            diag1[i] = j[i][i + 1];
            j[i + 1][i] = j[i][i + 1];
        }

        j[2 * N][2 * N] = ab0[0][2 * N];
        EigenDecomposition eig = new EigenDecomposition(diag2, diag1);
        RealMatrix V = eig.getV();
        RealMatrix D = eig.getD();
        double[] d2 = Diag(D.getData());
        double[] e = new double[ab0[0].length];
        for (int i = 0; i < e.length; i++) {
            e[i] = ab0[1][0] * (Math.pow(V.getEntry(0, i), 2));
        }

        Arrays.sort(d2);
        Double[] abscissae = new Double[d2.length];
        Double[] weights = new Double[e.length];
        for (int i = 0; i < d2.length; i++) {
            abscissae[i] = 2 * d2[i] - 1;
            weights[i] = 2 * e[i];
        }

        return new Pair<>(abscissae, weights);
    }

    /**
     *
     * @param arr
     * @return
     */
    private static double[] CumulativeSum(double[] arr){
        double[] out = new double[arr.length];
        double total = 0;
        for (int i = 0; i < arr.length; i++) {
            total += arr[i];
            out[i] = total;
        }
        return out;
    }

    /**
     *
     * @param arr
     * @return
     */
    private static double[] Diag(double[][] arr){
        double[] out = new double[arr.length];
        for (int i = 0; i < arr.length; i++) {
            out[i] = arr[i][i];
        }
        return out;
    }
}

package lab;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

public class ApproximatingFuncFinder {
    private final RealMatrix xMatrix;
    private RealMatrix yTeorMatrix;
    private final RealMatrix yActMatrix;
    private RealMatrix xAllPowersMatrix;
    private final int maxPower;

    private RealMatrix bMatrix;

    private double lSK = Double.MAX_VALUE;

    public ApproximatingFuncFinder(RealMatrix xMatrix, RealMatrix yActMatrix, int maxPower) {
        this.xMatrix = xMatrix;
        this.yActMatrix = yActMatrix;
        this.maxPower = maxPower;

        createXAllPowersMatrix();
    }

    private void createXAllPowersMatrix() {
        double[][] xAllPM = new double[xMatrix.getRowDimension()][maxPower+1];
        for(int i=0; i<xAllPM.length; i++) {
            for(int j=0; j<xAllPM[i].length; j++) {
                if(j == 0) {
                    xAllPM[i][j] = 1;
                }
                else {
                    xAllPM[i][j] = Math.pow(xMatrix.getEntry(i, 0), j);
                }
            }
        }
        xAllPowersMatrix = new Array2DRowRealMatrix(xAllPM);
    }

    public double calcAndGetLeastSquaresKr() {
        calcBMatrix();
        calcYTeorMatrix();
        double lSK = calcLeastSquaresKr();
        this.lSK = lSK;
        System.out.println("Max power: "+ maxPower +  ", LSK:" + lSK);
        return lSK;
    }

    public void makeEvaluationAndShowResult() {
        double corIndex = calcCorIndex();
        double[] dArr = calcDArr();
        double disp = calcDispersion();
        double[] tArr = calcTArr(dArr, disp);
        double tCritical = new TCritical().getT(xMatrix.getRowDimension()-maxPower-1);
        double[] paramAc = calcParametersAccuracy(tCritical, disp, dArr);

        System.out.println();
        System.out.println("Min least square criterion: " + lSK);
        System.out.println("Max power: " + maxPower);
        for(int i=0; i<bMatrix.getRowDimension(); i++) {
            if(i==0)
                System.out.print("y = " + bMatrix.getEntry(i, 0) + " ");
            else if(i<bMatrix.getRowDimension()-1)
                System.out.print(bMatrix.getEntry(i, 0) + "x^" + i + " + ");
            else
                System.out.println(bMatrix.getEntry(i, 0) + "x^" + i);
        }
        System.out.println("corIndex: " + String.format("%.2f", corIndex));
        System.out.print("D: ");
        for (double v : dArr) {
            System.out.print(v + " ");
        }
        System.out.println();
        System.out.println("Dispersion: " + disp);

        System.out.println("T critical: " + tCritical);
        System.out.print("T: ");
        for(double v: tArr) {
            System.out.print(v + " ");
        }
        System.out.println();
        System.out.println("Parameters accuracy: ");
        for(double v: paramAc) {
            System.out.print(v + " ");
        }
        System.out.println();
    }

    public void calcBMatrix() {
        RealMatrix transposeM = xAllPowersMatrix.transpose();

        RealMatrix tempM = transposeM.multiply(xAllPowersMatrix);
        RealMatrix inverseM = new LUDecomposition(tempM).getSolver().getInverse();
        bMatrix = (inverseM.multiply(transposeM)).multiply(yActMatrix);
    }

    public void calcYTeorMatrix() {
        if(bMatrix == null)
            calcBMatrix();

        double[] yTM = new double[yActMatrix.getRowDimension()];
        for(int i=0; i<yTM.length; i++) {
            double sum = 0d;
            for(int j=0; j<bMatrix.getRowDimension(); j++) {
                sum += bMatrix.getEntry(j, 0)*xAllPowersMatrix.getEntry(i, j);
            }
            yTM[i] = sum;
        }

        yTeorMatrix = new Array2DRowRealMatrix(yTM);
    }

    private double calcLeastSquaresKr() {
        double sum = 0d;
        for(int i=0; i< yActMatrix.getRowDimension(); i++) {
            sum += Math.pow(yTeorMatrix.getEntry(i, 0) - yActMatrix.getEntry(i, 0), 2);
        }
        return sum;
    }

    public double calcCorIndex() {
        double avY = 0d;
        for(int i=0; i< yActMatrix.getRowDimension(); i++) {
            avY += yActMatrix.getEntry(i, 0);
        }
        avY /= yActMatrix.getRowDimension();

        double sigmaFact2 = 0d, sigmaZag2 = 0d;
        for(int i=0; i<yActMatrix.getRowDimension(); i++) {
            sigmaFact2 += Math.pow(yTeorMatrix.getEntry(i, 0) - avY, 2);
            sigmaZag2 += Math.pow(yActMatrix.getEntry(i, 0) - avY, 2);
        }
        int sm1 = yTeorMatrix.getRowDimension() - 1;
        sigmaFact2 /= sm1;
        sigmaZag2 /= sm1;

        return Math.sqrt(sigmaFact2/sigmaZag2);
    }

    public double[] calcDArr() {
        double[] dArr = new double[xAllPowersMatrix.getColumnDimension()];
        RealMatrix transposeM = xAllPowersMatrix.transpose();
        RealMatrix tempM = transposeM.multiply(xAllPowersMatrix);
        RealMatrix inverseM = new LUDecomposition(tempM).getSolver().getInverse();
        for(int i=0; i<dArr.length; i++) {
            dArr[i] = inverseM.getEntry(i, i);
        }
        return dArr;
    }

    public double calcDispersion() {
        double sum = 0d;
        for(int i=0; i<yTeorMatrix.getRowDimension(); i++) {
            sum += Math.pow(yTeorMatrix.getEntry(i, 0) - yActMatrix.getEntry(i, 0), 2);
        }
        sum = sum/(yTeorMatrix.getRowDimension()-maxPower-1);
        return sum;
    }

    public double[] calcTArr(double[] dArr, double dispersion) {
        double[] tArr = new double[dArr.length];
        for(int i=0; i<tArr.length; i++) {
            tArr[i] = Math.abs(bMatrix.getEntry(i, 0))/Math.sqrt(dArr[i]*dispersion);
        }
        return tArr;
    }

    public double[] calcParametersAccuracy(double tCr, double disp, double[] d) {
        double[] paramAc = new double[d.length];
        for(int i=0; i<paramAc.length; i++) {
            paramAc[i] = tCr*Math.sqrt(disp*d[i]);
        }
        return paramAc;
    }
}

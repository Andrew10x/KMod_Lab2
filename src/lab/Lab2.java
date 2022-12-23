package lab;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;

public class Lab2 {
    public static void main(String[] args) {
        doTask();
    }

    public static void doTask() {
        double[] xM = {1,5,9,13,17,21,25,29,33,37,41,45,49,53,57,61};
        double[] yTeorM = { 312.89, 1612, 4225, 8043, 12900, 18560, 24740,
                31070, 37160, 42510, 46600, 48820, 48510, 44960, 37370, 24910};

        ApproximatingFuncFinder[] wArr = new ApproximatingFuncFinder[10];
        int minPos = 0;
        double minLSK = Double.MAX_VALUE;
        for(int i=0; i< wArr.length; i++) {
            wArr[i] = new ApproximatingFuncFinder(new Array2DRowRealMatrix(xM), new Array2DRowRealMatrix(yTeorM), i+1);
            double curLSK = wArr[i].calcAndGetLeastSquaresKr();
            //wArr[i].makeEvaluationAndShowResult();
            if(minLSK > curLSK) {
                minLSK = curLSK;
                minPos = i;
            }
        }

        wArr[minPos].makeEvaluationAndShowResult();
    }


}

package lab;

public class TCritical {
    private final double[] tArr = {12.7, 4.303, 3.182, 2.776, 2.571, 2.447, 2.365,
            2.306, 2.262, 2.228, 2.201, 2.179, 2.16, 2.145, 2.131, 2.12, 2.11, 2.101,
            2.093, 2.086};

    public double getT(int pos) {
        pos -= 1;
        if(pos>tArr.length-1)
            pos=tArr.length-1;

        return tArr[pos];
    }
}

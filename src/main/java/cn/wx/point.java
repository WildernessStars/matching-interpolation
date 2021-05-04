package cn.wx;

public class point {
    private double x;
    private double y;

    public point(double a,double b) throws Exception {
        this.x=a;
        this.y=b;
    }


    public void setX(double x) {
        this.x = x;
    }

    public void setY(double y) {
        this.y = y;
    }

    public double getX(){
        return x;
    }

    public double getY() {
        return y;
    }


}

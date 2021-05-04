package cn.wx;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.analysis.FunctionUtils;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.analysis.integration.UnivariateIntegrator;

public class simpson {
    private static double a;
    private static double b;
    private static double c;
    private static double d;
    private static double x0;
    private static double x1;

    public void setA(double a) {
        this.a = a;
    }
    public void setB(double b) {
        this.b = b;
    }
    public void setC(double c) {
        this.c = c;
    }
    public void setD(double d) {
        this.d = d;
    }
    public void setX0(double x0) {
        this.x0 = x0;
    }
    public void setX1(double x1){
        this.x1 = x1;
    }

    private static double lower;
    private static double upper;    // 1. 2. 4

    public double running() throws Exception {
        func = new MyFunction();
        integrator = new SimpsonIntegrator();

        lower=x0;upper=x1;
        double result = ((SimpsonIntegrator) (integrator)).integrate(steps, (UnivariateFunction) (func), lower, upper);
        return result;
    }


    UnivariateFunction func = null;
    UnivariateIntegrator integrator = null;

    private static final int steps = 1000000;

    class MyFunction implements UnivariateFunction {
        public double value(double x) {
            double y = Math.sqrt(1+Math.pow(3*a*x*x + 2*b*x + c ,2));        // 4.
            return y;
        }
        private static final double factor = 0.0001;
    }

}

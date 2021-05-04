package cn.wx;

import java.util.ArrayList;
import java.util.List;
import org.apache.commons.lang3.*;
import org.apache.commons.math3.*;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.ArrayRealVector;


public class interpolation {
    private static double[][] matrix;
    private static double[][] pola;

    public interpolation() throws Exception {
        //this.linear();
        //this.fit();
    }

    public static double[][] linear() throws Exception {
        /*遍历找缺失*/
        double[][] mat = readmat.getData();
        pola = new double[2][256];
        int index = 0;
        for (int i = 0; i < mat.length; i++) {
            if (mat[i][3] < 0.1) {
                int j = 0;
                while (mat[i + j][3] < 0.1) {
                    j++;
                }
                double b1 = (mat[i + j][2] - mat[i - 1][2]) / (mat[i + j][0] - mat[i - 1][0]);
                double a1 = (mat[i - 1][0] * mat[i + j][2] - mat[i - 1][2] * mat[i + j][0]) / (mat[i + j][0] - mat[i - 1][0]);

                double b2 = (mat[i + j][1] - mat[i - 1][1]) / (mat[i + j][0] - mat[i - 1][0]);
                double a2 = (mat[i - 1][0] * mat[i + j][1] - mat[i - 1][1] * mat[i + j][0]) / (mat[i + j][0] - mat[i - 1][0]);

                j = 0;
                while (mat[i + j][3] < 0.1) {
                    pola[0][index] = -a2 + b2 * mat[i + j][0];
                    pola[1][index] = -a1 + b1 * mat[i + j][0];

                    j++;
                    index++;
                }
            }
        }
        // System.out.print(len+1);
        return pola;
    }

    /*三次样条加弧长与时间的关系*/
    public static double[][] three() throws Exception {
        double[][] mat = readmat.getData();
        int len = 0;
        for (int i = 0; i < mat.length; i++) {
            if (Math.abs(mat[i][3] - 1) < 0.001)
                len++;
        }
        double[][] tran = new double[3][len];
        int ind = 0;
        for (int i = 0; i < mat.length; i++) {
            if (Math.abs(mat[i][3] - 1) < 0.001) {
                tran[0][ind] = mat[i][0];
                tran[1][ind] = mat[i][1];
                tran[2][ind] = mat[i][2];
                ind++;
            }
        }
        double t0, t1, t2, t3;
        double x0, x1, x2, x3;
        double y0, y1, y2, y3;
        pola = new double[2][256];
        int index = 0;
        int num0 = 0;
        for (int i = 0; i < mat.length; i++) {
            if (mat[i][3] < 0.1) {
                if (i >= 363) {
                    t0 = mat[352][0];
                    x0 = mat[352][1];
                    y0 = mat[352][2];
                    t1 = mat[354][0];
                    x1 = mat[354][1];
                    y1 = mat[354][2];
                    t2 = mat[362][0];
                    x2 = mat[362][1];
                    y2 = mat[362][2];
                    t3 = mat[365][0];
                    x3 = mat[365][1];
                    y3 = mat[365][2];
                } else if (i < 2) {
                    t0 = mat[0][0];
                    x0 = mat[0][1];
                    y0 = mat[0][2];
                    t1 = mat[2][0];
                    x1 = mat[2][1];
                    y1 = mat[2][2];
                    t2 = mat[3][0];
                    x2 = mat[3][1];
                    y2 = mat[3][2];
                    t3 = mat[4][0];
                    x3 = mat[4][1];
                    y3 = mat[4][2];
                } else {
                    int j = 2;
                    while (true) {
                        if (mat[i - j][3] > 0.9)
                            break;
                        j++;
                    }
                    t0 = mat[i - j][0];
                    x0 = mat[i - j][1];
                    y0 = mat[i - j][2];
                    t1 = mat[i - 1][0];
                    x1 = mat[i - 1][1];
                    y1 = mat[i - 1][2];
                    j = 1;
                    while (true) {
                        if (mat[i + j][3] > 0.9)
                            break;
                        j++;
                    }
                    t2 = mat[i + j][0];
                    x2 = mat[i + j][1];
                    y2 = mat[i + j][2];
                    while (true) {
                        if (mat[i + j + 1][3] > 0.9)
                            break;
                        j++;
                    }
                    t3 = mat[i + j + 1][0];
                    x3 = mat[i + j + 1][1];
                    y3 = mat[i + j + 1][2];
                }

              /*  RealMatrix coefficients =
                        new Array2DRowRealMatrix(new double[][] { { t0*t0*t0,t0*t0,t0,1,0,0,0,0,0,0,0,0 },
                                { t1*t1*t1,t1*t1,t1,1,0,0,0,0,0,0,0,0 }, { 3*t0,1,0,0,0,0,0,0,0,0,0,0 },
                                {3*t1*t1,2*t1,1,0,-3*t1*t1,-2*t1,-1,0,0,0,0,0},{3*t1,1,0,0,-3*t1,-1,0,0,0,0,0,0},
                                {0,0,0,0,t1*t1*t1,t1*t1,t1,1,0,0,0,0},{0,0,0,0,t2*t2*t2,t2*t2,t2,1,0,0,0,0},
                                {0,0,0,0,3*t2*t2,2*t2,1,0,-3*t2*t2,-2*t2,-1,0},{0,0,0,0,3*t2,1,0,0,-3*t2,-1,0,0},
                                {0,0,0,0,0,0,0,0,t2*t2*t2,t2*t2,t2,1},{0,0,0,0,0,0,0,0,t3*t3*t3,t3*t3,t3,1},
                                {0,0,0,0,0,0,0,0,3*t3,1,0,0}}, false);
                DecompositionSolver solver = new LUDecomposition(coefficients).getSolver();
                RealVector constants = new ArrayRealVector(new double[] { x0, x1, 0,0,0,x1,x2,0,0,x2,x3,0 }, false);
                RealVector solution = solver.solve(constants);*/


                RealMatrix coefficients2 =
                        new Array2DRowRealMatrix(new double[][]{{x0 * x0 * x0, x0 * x0, x0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
                                {x1 * x1 * x1, x1 * x1, x1, 1, 0, 0, 0, 0, 0, 0, 0, 0}, {3 * x0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                {3 * x1 * x1, 2 * x1, 1, 0, -3 * x1 * x1, -2 * x1, -1, 0, 0, 0, 0, 0}, {3 * x1, 1, 0, 0, -3 * x1, -1, 0, 0, 0, 0, 0, 0},
                                {0, 0, 0, 0, x1 * x1 * x1, x1 * x1, x1, 1, 0, 0, 0, 0}, {0, 0, 0, 0, x2 * x2 * x2, x2 * x2, x2, 1, 0, 0, 0, 0},
                                {0, 0, 0, 0, 3 * x2 * x2, 2 * x2, 1, 0, -3 * x2 * x2, -2 * x2, -1, 0}, {0, 0, 0, 0, 3 * x2, 1, 0, 0, -3 * x2, -1, 0, 0},
                                {0, 0, 0, 0, 0, 0, 0, 0, x2 * x2 * x2, x2 * x2, x2, 1}, {0, 0, 0, 0, 0, 0, 0, 0, x3 * x3 * x3, x3 * x3, x3, 1},
                                {0, 0, 0, 0, 0, 0, 0, 0, 3 * x3, 1, 0, 0}}, false);
                DecompositionSolver solver2 = new LUDecomposition(coefficients2).getSolver();
                RealVector constants2 = new ArrayRealVector(new double[]{y0, y1, 0, 0, 0, y1, y2, 0, 0, y2, y3, 0}, false);
                RealVector solution2 = solver2.solve(constants2);
                int j = 0;
                //System.out.println("========================================================");
                simpson sim = new simpson();

                double star = 0;
                while (mat[i + j][3] < 0.1) {
                    if (i >= 363) {
                        double dx;
                        double sumtime = t3 - t2;
                        double mytime = mat[i + j][0] - t2;

                        sim.setA(solution2.getEntry(8));
                        sim.setB(solution2.getEntry(9));
                        sim.setC(solution2.getEntry(10));
                        sim.setD(solution2.getEntry(11));
                        sim.setX0(x2);
                        sim.setX1(x3);
                        double sums = sim.running();
                        double ttemp = mytime * sums;
                        for (dx = x2 + 0.5; dx < x3; dx += 0.5) {
                            sim.setA(solution2.getEntry(8));
                            sim.setB(solution2.getEntry(9));
                            sim.setC(solution2.getEntry(10));
                            sim.setD(solution2.getEntry(11));
                            sim.setX0(x2);
                            sim.setX1(dx);
                            double mys = sim.running();
                            if (Math.abs(ttemp - mys * sumtime) < 1) {
                                pola[0][index] = dx;
                                break;
                            }
                        }
                        //System.out.println(pola[0][index]);
                        pola[1][index] = solution2.getEntry(8) * Math.pow(pola[0][index], 3)
                                + solution2.getEntry(9) * Math.pow(pola[0][index], 2) + solution2.getEntry(10) * pola[0][index] + solution2.getEntry(11);
                    } else if (i < 2) {
                        double dx;
                        double sumtime = t1 - t0;
                        double mytime = mat[i + j][0] - t0;

                        sim.setA(solution2.getEntry(0));
                        sim.setB(solution2.getEntry(1));
                        sim.setC(solution2.getEntry(2));
                        sim.setD(solution2.getEntry(3));
                        sim.setX0(x0);
                        sim.setX1(x1);
                        double sums = sim.running();
                        double ttemp = mytime * sums;
                        for (dx = x0 + 1; dx < x2; dx += 0.5) {
                            sim.setA(solution2.getEntry(0));
                            sim.setB(solution2.getEntry(1));
                            sim.setC(solution2.getEntry(2));
                            sim.setD(solution2.getEntry(3));
                            sim.setX0(x0);
                            sim.setX1(dx);
                            double mys = sim.running();
                            if (Math.abs(ttemp - mys * sumtime) < 0.06) {
                                pola[0][index] = dx;
                                break;
                            }
                        }
                        //System.out.println(pola[0][index]);
                        pola[1][index] = solution2.getEntry(0) * Math.pow(pola[0][index], 3)
                                + solution2.getEntry(1) * Math.pow(pola[0][index], 2) + solution2.getEntry(2) * pola[0][index] + solution2.getEntry(3);
                    } else {
                        double sumtime = t2 - t1;


                        sim.setA(solution2.getEntry(4));
                        sim.setB(solution2.getEntry(5));
                        sim.setC(solution2.getEntry(6));
                        sim.setD(solution2.getEntry(7));
                        sim.setX0(x1);
                        sim.setX1(x2);
                        if (x1 > x2) {
                            sim.setX0(x2);
                            sim.setX1(x1);
                        }
                        double sums = sim.running();
                        double mytime = mat[i + j][0] - t1;
                        double ttemp = sums / sumtime;
                        double dx;
                        double minabs = 100;
                        if (x1 < x2) {
                            if (star == 0) dx = x1 + 0.02;
                            else dx = star + 0.1;
                            for (; dx < x2; dx += 0.02) {
                                sim.setA(solution2.getEntry(4));
                                sim.setB(solution2.getEntry(5));
                                sim.setC(solution2.getEntry(6));
                                sim.setD(solution2.getEntry(7));
                                sim.setX0(x1);
                                sim.setX1(dx);
                                double mys = sim.running();
                                if (Math.abs(ttemp - mys / mytime) < minabs) {
                                    minabs = Math.abs(ttemp - mys / mytime);
                                    pola[0][index] = dx;
                                } else
                                    break;
                            }
                        } else {
                            if (star == 0) dx = x1 - 0.1;
                            else dx = star - 0.1;
                            for (; dx > x2; dx -= 0.1) {
                                sim.setA(solution2.getEntry(4));
                                sim.setB(solution2.getEntry(5));
                                sim.setC(solution2.getEntry(6));
                                sim.setD(solution2.getEntry(7));
                                sim.setX0(dx);
                                sim.setX1(x1);
                                double mys = sim.running();
                                if (Math.abs(ttemp - mys / mytime) < minabs) {
                                    minabs = Math.abs(ttemp - mys / mytime);
                                    pola[0][index] = dx;
                                } else
                                    break;
                            }
                        }
                        star = pola[0][index];
                        //System.out.println(pola[0][index]);
                        pola[1][index] = solution2.getEntry(4) * Math.pow(pola[0][index], 3)
                                + solution2.getEntry(5) * Math.pow(pola[0][index], 2) + solution2.getEntry(6) * pola[0][index] + solution2.getEntry(7);
                    }
                    index++;
                    j++;
                }
                i += j - 1;


            }
        }
        // System.out.print(len+1);
        return pola;
    }

    public static double[][] poly() throws Exception {
        double[][] mat = readmat.getData();
        //matrix =filter(4,matrix,0,false);
        double t0, t1, t2, t3, x0, x1, x2, x3, y0, y1, y2, y3;
        double a, b, c, d;
        pola = new double[2][256];
        int index = 0;
        for (int i = 0; i < mat.length; i++) {
            if (mat[i][3] < 0.1) {
                if (i >= 363) {
                    t0 = mat[352][0];
                    x0 = mat[352][1];
                    y0 = mat[352][2];
                    t1 = mat[354][0];
                    x1 = mat[354][1];
                    y1 = mat[354][2];
                    t2 = mat[362][0];
                    x2 = mat[362][1];
                    y2 = mat[362][2];
                    t3 = mat[365][0];
                    x3 = mat[365][1];
                    y3 = mat[365][2];
                } else if (i < 2) {
                    t0 = mat[0][0];
                    x0 = mat[0][1];
                    y0 = mat[0][2];
                    t1 = mat[2][0];
                    x1 = mat[2][1];
                    y1 = mat[2][2];
                    t2 = mat[3][0];
                    x2 = mat[3][1];
                    y2 = mat[3][2];
                    t3 = mat[4][0];
                    x3 = mat[4][1];
                    y3 = mat[4][2];
                } else {
                    int j = 2;
                    while (true) {
                        if (mat[i - j][3] > 0.9)
                            break;
                        j++;
                    }
                    t0 = mat[i - j][0];
                    x0 = mat[i - j][1];
                    y0 = mat[i - j][2];
                    t1 = mat[i - 1][0];
                    x1 = mat[i - 1][1];
                    y1 = mat[i - 1][2];
                    j = 1;
                    while (true) {
                        if (mat[i + j][3] > 0.9)
                            break;
                        j++;
                    }
                    t2 = mat[i + j][0];
                    x2 = mat[i + j][1];
                    y2 = mat[i + j][2];
                    while (true) {
                        if (mat[i + j + 1][3] > 0.9)
                            break;
                        j++;
                    }
                    t3 = mat[i + j + 1][0];
                    x3 = mat[i + j + 1][1];
                    y3 = mat[i + j + 1][2];
                }
                pola[0][index]=third(mat[i][0],t0,t1,t2,t3,x0,x1,x2,x3);
                pola[1][index]=third(mat[i][0],t0,t1,t2,t3,y0,y1,y2,y3);
                index++;
            }


        }
        return pola;
    }
    public static double third(double x,double x0,double x1,double x2,double x3,double y0,double y1,double y2,double y3) {
        double a,b,c,d;
        a = 0 - x1 * x2 * x3 * y0 / ((x0 - x1) * (x0 - x2) * (x0 - x3)) - x0 * x2 * x3 * y1 / ((x1 - x0) * (x1 - x2) * (x1 - x3)) - x0 * x1 * x3 * y2 / ((x2 - x0) * (x2 - x1) * (x2 - x3)) - x0 * x1 * x2 * y3 / ((x3 - x0) * (x3 - x1) * (x3 - x2));
        b = y0 * (x1 * x2 + x2 * x3 + x3 * x1) / ((x0 - x1) * (x0 - x2) * (x0 - x3)) + y1 * (x0 * x2 + x2 * x3 + x3 * x0) / ((x1 - x0) * (x1 - x2) * (x1 - x3)) + y2 * (x0 * x1 + x1 * x3 + x3 * x0) / ((x2 - x0) * (x2 - x1) * (x2 - x3)) + y3 * (x0 * x1 + x1 * x2 + x2 * x0) / ((x3 - x0) * (x3 - x1) * (x3 - x2));
        c = 0 - y0 * (x1 + x2 + x3) / ((x0 - x1) * (x0 - x2) * (x0 - x3)) - y1 * (x0 + x2 + x3) / ((x1 - x0) * (x1 - x2) * (x1 - x3)) - y2 * (x0 + x1 + x3) / ((x2 - x0) * (x2 - x1) * (x2 - x3)) - y3 * (x0 + x1 + x2) / ((x3 - x0) * (x3 - x1) * (x3 - x2));
        d = y0 / ((x0 - x1) * (x0 - x2) * (x0 - x3)) + y1 / ((x1 - x0) * (x1 - x2) * (x1 - x3)) + y2 / ((x2 - x0) * (x2 - x1) * (x2 - x3)) + y3 / ((x3 - x0) * (x3 - x1) * (x3 - x2));
        return a+b*x+c*x*x+d*x*x*x;
    }
}



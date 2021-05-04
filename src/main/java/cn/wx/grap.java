package cn.wx;


import java.awt.BasicStroke;
import java.awt.Color;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.DefaultXYDataset;
import java.awt.geom.Ellipse2D;

public class grap {
    private static double[][] matrix;
    private static double[][] filtered;
    private static double[][] method;

    public grap() throws Exception {
        this.showChart();
    }

    public void showChart() throws Exception {
        DefaultXYDataset xydataset = new DefaultXYDataset();
        matrix = readmat.filtering();

        method = interpolation.three();//linear,three,poly
        xydataset.addSeries("(x,y)", matrix);
        xydataset.addSeries("new(x,y)", method);
        //xydataset.addSeries("newnew(x,y)", fitted);
        JFreeChart chart = ChartFactory.createScatterPlot("", "x", "y",
                xydataset,
                PlotOrientation.VERTICAL,
                true,
                false,
                false);
        ChartFrame frame = new ChartFrame("", chart, true);
        chart.setBackgroundPaint(Color.white);


        XYPlot xyplot = (XYPlot) chart.getPlot();
        xyplot.setBackgroundPaint(new Color(255, 253, 246));
        xyplot.setDomainGridlinePaint(Color.gray); //网格线纵向颜色
        xyplot.setRangeGridlinePaint(Color.gray); //网格线横向颜色

        XYLineAndShapeRenderer xylineandshaperenderer = (XYLineAndShapeRenderer) xyplot.getRenderer();
        xylineandshaperenderer.setSeriesOutlinePaint(0, Color.WHITE);
        xylineandshaperenderer.setSeriesPaint(0, Color.blue);
        xylineandshaperenderer.setUseOutlinePaint(true);
        xylineandshaperenderer.setSeriesShape(0,  new Ellipse2D.Double(-1.5D, -1.5D, 3.0D, 3.0D));
        xylineandshaperenderer.setSeriesStroke(0,new BasicStroke(1F));


        xylineandshaperenderer.setSeriesOutlinePaint(1, Color.red);
        xylineandshaperenderer.setSeriesPaint(1, Color.red);
        xylineandshaperenderer.setSeriesShape(1,  new Ellipse2D.Double(-0.5D, -0.5D, 1.0D, 1.0D));
        xylineandshaperenderer.setSeriesStroke(1,new BasicStroke(1F));

        //BasicStroke realLine2 = new BasicStroke(1.0f);
        //xylineandshaperenderer.setSeriesStroke(1, realLine2);
        //xylineandshaperenderer.setSeriesPaint(1, Color.red);
        //xylineandshaperenderer.setSeriesLinesVisible(1,true);

        NumberAxis numberaxis = (NumberAxis) xyplot.getDomainAxis();
        numberaxis.setAutoRangeIncludesZero(false);
        numberaxis.setTickMarkInsideLength(2.0F);
        numberaxis.setTickMarkOutsideLength(0.0F);
        numberaxis.setAxisLineStroke(new BasicStroke(1.5f));

        frame.pack();
        frame.setVisible(true);

    }

    public static void main(String[] args) throws Exception {

        grap my = new grap();
    }
}




/*import com.jmatio.io.*;
import com.jmatio.types.MLDouble;

import java.io.IOException;
import java.util.ArrayList;

public class readmat {
    public static void main(String[] args) throws IOException {

        double[][] matTest=new double[][]{{1,2,3,4},{5,6,76,34}};//生成待存储的矩阵
        MLDouble mlDouble=new MLDouble("doubleArray",matTest);//doubleArray就是matlab中上述矩阵的标示符，load()之后，在matlab中使用doubleArray访问此矩阵
        ArrayList list=new ArrayList();//由于MatFileWriter()构造函数的参数为list类型，所以需要创建一个ArrayList
        list.add(mlDouble);
        new MatFileWriter("matTest.mat",list);//将矩阵写入到.mat文件中，文件名为matTest.mat
        System.out.println("mat writer done!");
    }
}*/
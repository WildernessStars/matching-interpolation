package cn.wx;

import com.jmatio.io.*;
import com.jmatio.types.MLDouble;
import java.util.ArrayList;
import com.jmatio.types.MLArray;
import com.jmatio.types.MLDouble;
import java.io.IOException;


public class readmat {
    public static double[][] getData() throws Exception {
        MatFileReader reader = new MatFileReader("data/homework.mat");
        MLArray mlArray = reader.getMLArray("coord");
        MLDouble d = (MLDouble) mlArray;
        int sizex = mlArray.getDimensions()[0];
        int sizey = mlArray.getDimensions()[1];
        //System.out.println(sizex+" "+sizey);
        double[][] matrix = (d.getArray());//只有jmatio v0.2版本中才有d.getArray方法

       for (int i = 0; i < sizex; i++) {
            for (int j = 0; j < sizey; j++) {
                System.out.print(matrix[i][j] + "\t");
            }

            System.out.println();
        }

        System.out.print(matrix.length + "\t");
        return matrix;
    }

 public static double[][] filtering() throws Exception {
     double[][] data = getData();

     int len = 0;
     for (int i = 0; i < data.length; i++) {
         if (Math.abs(data[i][3] - 1) < 0.001)
             len++;
     }
     double[][] matrix_tran = new double[2][len];
     int ind = 0;
     for (int i = 0; i < data.length; i++) {
         if (Math.abs(data[i][3] - 1) < 0.001) {
             matrix_tran[0][ind] = data[i][1];
             matrix_tran[1][ind] = data[i][2];
             ind++;
         }
     }


     //  System.out.print(matrix.length + "\t");
     return matrix_tran;
 }
}




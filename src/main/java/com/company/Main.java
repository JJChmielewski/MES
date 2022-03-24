package com.company;

import com.company.structs.*;

import java.io.*;
import java.util.Arrays;

public class Main {

    public static void main(String[] args) throws IOException {

        double H,B,k,alpha,tempOS,ro,cp,dt,temp0,tend;
        int nB,nH, numberOfPoints;

        H=0.1;
        B=0.1;
        nH=4;
        nB=4;
        k=25;
        numberOfPoints=3;
        alpha=300;
        tempOS=1200;
        ro=7800;
        cp=700;
        dt=50;
        temp0=100;
        tend=500;

        BufferedWriter bufferedWriter = new BufferedWriter(new OutputStreamWriter(new FileOutputStream("mes.csv"),"UTF-8"));


        Grid grid = new Grid(H,B,nH,nB);

        double[] T = new double[grid.nN];


        for(int i=0;i< grid.nN;i++){
            T[i] =temp0;
        }

        for(double i=dt; i<=tend;i+=dt){
            StringBuffer stringBuffer = new StringBuffer();

            System.out.println("Iteracja: "+i);

            grid = new Grid(H,B,nH,nB);

            grid.countLocalH(k,numberOfPoints);

            grid.countHbcAndP(alpha, tempOS, numberOfPoints);

            grid.countLocalC(ro, cp, numberOfPoints);

            grid.countGlobalHandCandP();

            grid.countGlobalHandPwithC( dt, T);


            //System.out.println(Arrays.toString(grid.elements));

            //grid.printGlobalHandC();

            T = grid.countT();

            System.out.println(Arrays.toString(T));

            //System.out.println(Arrays.toString(grid.elements));

            //System.out.println(Arrays.stream(T).max() + ", "+ Arrays.stream(T).min());

            stringBuffer.append(i);
            stringBuffer.append(";");
            stringBuffer.append(Arrays.stream(T).max().getAsDouble());
            stringBuffer.append(";");
            stringBuffer.append(Arrays.stream(T).min().getAsDouble());
            bufferedWriter.write(stringBuffer.toString());
            bufferedWriter.newLine();
        }
        bufferedWriter.flush();
        bufferedWriter.close();

    }
}

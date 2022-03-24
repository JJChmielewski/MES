package com.company.structs;

import java.util.Arrays;

public class Element {

    public int[] id;
    public int numberOfElement;
    public Node[] nodes;
    public double[][] Hbc;
    public double [] P;
    public double[][] H;
    public double[][] C;

    public Element() {
        this.id = new int[4];
        nodes = new Node[4];
        this.Hbc = new double[4][4];
        this.P = new double[4];
    }

    public Element(int[] id) {
        this.id = id;
    }

    public void countH(Element4_2D element4_2D, double k){

        this.H = new double[4][4];
        Gauss gauss = new Gauss(element4_2D.numberOfPoints);

        for(int pc=0;pc<element4_2D.numberOfPoints*element4_2D.numberOfPoints;pc++){
            Jakobian jakobian = new Jakobian(this, element4_2D,pc);
            element4_2D.countDerivativesXY(jakobian);
            for(int i=0;i<4;i++){
                for(int j=0;j<4;j++){
                    this.H[i][j] += jakobian.det * k * gauss.weights[pc/element4_2D.numberOfPoints] * gauss.weights[pc%element4_2D.numberOfPoints] * ( element4_2D.dNdx[pc][i]*element4_2D.dNdx[pc][j] + element4_2D.dNdy[pc][i]*element4_2D.dNdy[pc][j]);
                }
            }
        }

    }

    public void countC(double ro, double cp, int numberOfPoints){

        Gauss gauss = new Gauss(numberOfPoints);
        Element4_2D element4_2D = new Element4_2D(numberOfPoints);
        double[][] N = element4_2D.N;
        double[][] C = new double[4][4];

        for(int p=0;p<numberOfPoints*numberOfPoints;p++){
            Jakobian jakobian = new Jakobian(this,element4_2D,p);
            for(int x=0;x<4;x++){
                for(int y=0;y<4;y++){
                    C[x][y] += N[p][x] * N[p][y] *ro * cp * jakobian.det * gauss.weights[p%numberOfPoints] * gauss.weights[p/numberOfPoints];
                }
            }
        }

        this.C = C;

    }

    @Override
    public String toString() {
        String temp = "Element{" +
                "nodes=" + Arrays.toString(nodes);

        if(this.H != null){
            temp+="}, H={\n";
            for(int i=0;i<4;i++){
                temp += "\t"+ Arrays.toString(H[i]) +"\n";
            }
        }




        temp+="}, Hbc={\n";
        for(int i=0;i<4;i++){
            temp+=Arrays.toString(Hbc[i])+",\n";
        }


        temp+="P={\n";
        temp+=Arrays.toString(P)+",\n";

        temp+="}, numberOfElement=" + numberOfElement + '}';

        return temp;
    }
}

package com.company.structs;

import java.util.Arrays;

public class Element4_2D {

    public double[][] dNdKsi, dNdEta; //macierze funkcji kształtu lokalne
    public double[][] dNdx,dNdy; //macierze -||- globalne
    public double[][] N;

    public int numberOfPoints;

    public Element4_2D(int numberOfPoints) {

        Gauss gauss = new Gauss(numberOfPoints);

        this.numberOfPoints=numberOfPoints;
        this.dNdEta = new double[numberOfPoints*numberOfPoints][4];
        this.dNdKsi = new double[numberOfPoints*numberOfPoints][4];
        this.dNdx = null;
        this.dNdy = null;

        for(int i=0; i<numberOfPoints;i++){
            for(int j=0;j<numberOfPoints;j++){
                this.dNdKsi[i*(numberOfPoints)+j][0] = -1*(1-gauss.points[i])/4;
                this.dNdKsi[i*(numberOfPoints)+j][1] = (1-gauss.points[i])/4;
                this.dNdKsi[i*(numberOfPoints)+j][2] = -1*(1+gauss.points[i])/4; //zamiana indeksow 2 i 3 ze wzgledu na roznice w numeracji wezlow pomiedzy
                this.dNdKsi[i*(numberOfPoints)+j][3] = (1+gauss.points[i])/4;   //mna a prowadzacym
            }
        }

        for(int i=0; i<numberOfPoints;i++){
            for(int j=0;j<numberOfPoints;j++){
                this.dNdEta[i*(numberOfPoints)+j][0] = -1*(1-gauss.points[j])/4;
                this.dNdEta[i*(numberOfPoints)+j][1] = -1*(1+gauss.points[j])/4;
                this.dNdEta[i*(numberOfPoints)+j][2] = (1-gauss.points[j])/4; //sytuacja analogiczna jak powyzej
                this.dNdEta[i*(numberOfPoints)+j][3] = (1+gauss.points[j])/4;
            }
        }

        if(numberOfPoints ==2) {
            N = new double[4][4];

            N[0][0] = 0.25 * (1 - gauss.points[0]) * (1 - gauss.points[0]);
            N[0][1] = 0.25 * (1 + gauss.points[0]) * (1 - gauss.points[0]);
            N[0][2] = 0.25 * (1 - gauss.points[0]) * (1 + gauss.points[0]);
            N[0][3] = 0.25 * (1 + gauss.points[0]) * (1 + gauss.points[0]);

            N[1][0] = 0.25 * (1 - gauss.points[1]) * (1 - gauss.points[0]);
            N[1][1] = 0.25 * (1 + gauss.points[1]) * (1 - gauss.points[0]);
            N[1][2] = 0.25 * (1 - gauss.points[1]) * (1 + gauss.points[0]);
            N[1][3] = 0.25 * (1 + gauss.points[1]) * (1 + gauss.points[0]);

            N[2][0] = 0.25 * (1 - gauss.points[0]) * (1 - gauss.points[1]);
            N[2][1] = 0.25 * (1 + gauss.points[0]) * (1 - gauss.points[1]);
            N[2][2] = 0.25 * (1 - gauss.points[0]) * (1 + gauss.points[1]);
            N[2][3] = 0.25 * (1 + gauss.points[0]) * (1 + gauss.points[1]);

            N[3][0] = 0.25 * (1 - gauss.points[1]) * (1 - gauss.points[1]);
            N[3][1] = 0.25 * (1 + gauss.points[1]) * (1 - gauss.points[1]);
            N[3][2] = 0.25 * (1 - gauss.points[1]) * (1 + gauss.points[1]);
            N[3][3] = 0.25 * (1 + gauss.points[1]) * (1 + gauss.points[1]);
        }

        if(numberOfPoints ==3){
            N = new double[9][4];

            N[0][0] = 0.25 * (1 - gauss.points[0]) * (1 - gauss.points[0]);
            N[0][1] = 0.25 * (1 + gauss.points[0]) * (1 - gauss.points[0]);
            N[0][2] = 0.25 * (1 - gauss.points[0]) * (1 + gauss.points[0]);
            N[0][3] = 0.25 * (1 + gauss.points[0]) * (1 + gauss.points[0]);

            N[1][0] = 0.25 * (1 - gauss.points[1]) * (1 - gauss.points[0]);
            N[1][1] = 0.25 * (1 + gauss.points[1]) * (1 - gauss.points[0]);
            N[1][2] = 0.25 * (1 - gauss.points[1]) * (1 + gauss.points[0]);
            N[1][3] = 0.25 * (1 + gauss.points[1]) * (1 + gauss.points[0]);

            N[2][0] = 0.25 * (1 - gauss.points[2]) * (1 - gauss.points[0]);
            N[2][1] = 0.25 * (1 + gauss.points[2]) * (1 - gauss.points[0]);
            N[2][2] = 0.25 * (1 - gauss.points[2]) * (1 + gauss.points[0]);
            N[2][3] = 0.25 * (1 + gauss.points[2]) * (1 + gauss.points[0]);

            N[3][0] = 0.25 * (1 - gauss.points[0]) * (1 - gauss.points[1]);
            N[3][1] = 0.25 * (1 + gauss.points[0]) * (1 - gauss.points[1]);
            N[3][2] = 0.25 * (1 - gauss.points[0]) * (1 + gauss.points[1]);
            N[3][3] = 0.25 * (1 + gauss.points[0]) * (1 + gauss.points[1]);

            N[4][0] = 0.25 * (1 - gauss.points[1]) * (1 - gauss.points[1]);
            N[4][1] = 0.25 * (1 + gauss.points[1]) * (1 - gauss.points[1]);
            N[4][2] = 0.25 * (1 - gauss.points[1]) * (1 + gauss.points[1]);
            N[4][3] = 0.25 * (1 + gauss.points[1]) * (1 + gauss.points[1]);

            N[5][0] = 0.25 * (1 - gauss.points[2]) * (1 - gauss.points[1]);
            N[5][1] = 0.25 * (1 + gauss.points[2]) * (1 - gauss.points[1]);
            N[5][2] = 0.25 * (1 - gauss.points[2]) * (1 + gauss.points[1]);
            N[5][3] = 0.25 * (1 + gauss.points[2]) * (1 + gauss.points[1]);

            N[6][0] = 0.25 * (1 - gauss.points[0]) * (1 - gauss.points[2]);
            N[6][1] = 0.25 * (1 + gauss.points[0]) * (1 - gauss.points[2]);
            N[6][2] = 0.25 * (1 - gauss.points[0]) * (1 + gauss.points[2]);
            N[6][3] = 0.25 * (1 + gauss.points[0]) * (1 + gauss.points[2]);

            N[7][0] = 0.25 * (1 - gauss.points[1]) * (1 - gauss.points[2]);
            N[7][1] = 0.25 * (1 + gauss.points[1]) * (1 - gauss.points[2]);
            N[7][2] = 0.25 * (1 - gauss.points[1]) * (1 + gauss.points[2]);
            N[7][3] = 0.25 * (1 + gauss.points[1]) * (1 + gauss.points[2]);

            N[8][0] = 0.25 * (1 - gauss.points[2]) * (1 - gauss.points[2]);
            N[8][1] = 0.25 * (1 + gauss.points[2]) * (1 - gauss.points[2]);
            N[8][2] = 0.25 * (1 - gauss.points[2]) * (1 + gauss.points[2]);
            N[8][3] = 0.25 * (1 + gauss.points[2]) * (1 + gauss.points[2]);
        }
    }

    public void countDerivativesXY(Jakobian jakobian){

        this.dNdx = new double[this.numberOfPoints*this.numberOfPoints][4];
        this.dNdy = new double[this.numberOfPoints*this.numberOfPoints][4];

        jakobian.jakobian[0][0] *= (1/jakobian.det);
        jakobian.jakobian[0][1] *= (1/jakobian.det);
        jakobian.jakobian[1][0] *= (1/jakobian.det);
        jakobian.jakobian[1][1] *= (1/jakobian.det);

        //dla kazdego punktu calkowania obliczamy pochodna funkcji kształtu dla wspolrzednych globalnych za pomoca jakobianu
        for(int pc = 0;pc<numberOfPoints*numberOfPoints;pc++){
            for(int i=0;i<4;i++){
                this.dNdx[pc][i] = (jakobian.jakobian[0][0] * dNdKsi[pc][i]) + (jakobian.jakobian[0][1] * dNdEta[pc][i]);
                this.dNdy[pc][i] = (jakobian.jakobian[1][0] * dNdKsi[pc][i]) + (jakobian.jakobian[1][1] * dNdEta[pc][i]);
            }
        }
    }

    @Override
    public String toString() {

        String returned = "Element4_2D{" +
                "dNdKsi={\n";


        for(int i=0;i<this.numberOfPoints*numberOfPoints;i++){
            returned += "\t"+ Arrays.toString(dNdKsi[i]) +"\n";
        }
        returned+="}, dNdEta={\n";

        for(int i=0;i<this.numberOfPoints*numberOfPoints;i++){
            returned += "\t"+ Arrays.toString(dNdEta[i]) +"\n";
        }



        if(this.dNdx != null){
            returned+="}, dNdx={\n";
            for(int i=0;i<this.numberOfPoints*numberOfPoints;i++){
                returned += "\t"+ Arrays.toString(dNdx[i]) +"\n";
            }


        }

        if(this.dNdy != null){
            returned+="}, dNdy={\n";
            for(int i=0;i<this.numberOfPoints*numberOfPoints;i++){
                returned += "\t"+ Arrays.toString(dNdy[i]) +"\n";
            }


        }


        returned+="}}";

        return returned;
    }
}

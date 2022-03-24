package com.company.structs;


import java.util.Arrays;

import static java.lang.Math.sqrt;

public class Gauss {

    public double[] points;
    public double[] weights;

    public Gauss(int numberOfPoints) {

        this.points = new double[numberOfPoints];
        this.weights = new double[numberOfPoints];

        if(numberOfPoints == 2){
            this.points[0] = -1.0/sqrt(3);
            this.points[1] = 1.0/sqrt(3);
            this.weights[0] = 1.0;
            this.weights[1] = 1.0;
        }
        if(numberOfPoints == 3){
            this.points[0] = -sqrt(3.0/5.0);
            this.points[1] = 0;
            this.points[2] = sqrt(3.0/5.0);
            this.weights[0] = 5.0/9.0;
            this.weights[1] = 8.0/9.0;
            this.weights[2] = 5.0/9.0;
        }
    }


    public double function1D(double x){
        return 5*x*x+3*x+6;
    }

    public double function2D(double x, double y){
        return 5*x*x*y*y+3*x*y+6;
    }

    public double countIntegral1D(int numberOfPoints){

        double integral =0;

        for(int x=0;x<numberOfPoints;x++){
            integral += this.weights[x]*this.function1D(this.points[x]);
        }


        return integral;
    }

    public double countIntegral2D(int numberOfPoints){

        double integral =0;

        for(int x=0;x<numberOfPoints;x++){
            for(int y=0;y<numberOfPoints;y++){
                integral+=this.weights[x]*this.weights[y]*function2D(this.points[x],this.points[y]);
            }
        }
        return integral;
    }


    @Override
    public String toString() {
        return "Gauss{" +
                "points=" + Arrays.toString(points) +
                ", weights=" + Arrays.toString(weights) +
                '}';
    }
}

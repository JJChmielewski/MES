package com.company.structs;


import org.ejml.equation.Equation;
import org.ejml.simple.SimpleMatrix;

import java.util.Arrays;

public class Grid {

    public double H,B;

    public int nH,nB,nN,nE;

    public double dx,dy;

    public Node[] nodes;
    public Element[] elements;

    public double[][] globalH;
    public double[] globalP;
    public double[][] globalC;


    public Grid(double h, double b, int nH, int nB) {
        H = h;
        B = b;
        this.nH = nH;
        this.nB = nB;
        this.nN = nH*nB;
        this.nE = (nH-1)*(nB-1);
        this.dx = b/(nB-1);
        this.dy = h/(nH-1);
        this.globalH = new double[nN][nN];
        this.globalP = new double[nN];
        this.globalC = new double[nN][nN];

        this.nodes = new Node[nN];
        for(int i=0; i<nH;i++){
            for(int j=0;j<nB;j++){
                this.nodes[i*nB+j] = new Node();
                this.nodes[i*nB+j].numberOfNode=i*nB+j;
                this.nodes[i*nB+j].x=j*dx;
                this.nodes[i*nB+j].y=i*dy;
            }
        }

        this.elements = new Element[nE];
        for(int i=0; i<nH-1;i++){
            for(int j=0;j<nB-1;j++){
                this.elements[i*(nB-1)+j] = new Element();
                elements[i*(nB-1)+j].numberOfElement=i*(nB-1)+j;
                elements[i*(nB-1)+j].nodes[0] = this.nodes[i*nB+j];
                elements[i*(nB-1)+j].nodes[1] = this.nodes[i*nB+j+1];
                elements[i*(nB-1)+j].nodes[2] = this.nodes[(i+1)*nB+j];
                elements[i*(nB-1)+j].nodes[3] = this.nodes[(i+1)*nB+j+1];

            }
        }
    }

    public void countGlobalHandCandP(){

        globalH = new double[nN][nN];
        globalC = new double[nN][nN];
        globalP = new double[nN];

        for(int i=0;i<nE;i++){
            for(int x=0;x<4;x++){
                for(int y=0;y<4;y++){

                    globalH[elements[i].nodes[x].numberOfNode][elements[i].nodes[y].numberOfNode] += elements[i].H[x][y] +elements[i].Hbc[x][y];
                    globalC[elements[i].nodes[x].numberOfNode][elements[i].nodes[y].numberOfNode] += elements[i].C[x][y];
                }
                globalP[elements[i].nodes[x].numberOfNode] += elements[i].P[x];
            }
        }
    }

    public void countLocalH(double k, int numberOfPoints){

        Element4_2D element4_2D = new Element4_2D(numberOfPoints);

        for(int i=0;i< elements.length;i++){
            elements[i].countH(element4_2D,k);
        }
    }

    public void countLocalC(double ro, double cp, int numberOfPoints){
        for(int i=0;i< elements.length;i++){
            elements[i].countC(ro,cp,numberOfPoints);
        }
    }

    public void countHbcAndP(double alpha, double tempOfSurroundings, int numberOfPoints){
        Gauss gauss = new Gauss(numberOfPoints);
        double[][] N = new double[numberOfPoints][4];


        for(int i=0;i<nE;i++){
            elements[i].Hbc = new double[4][4];
            elements[i].P=new double[4];
        }

        //sciana 0, dolna
        for(int i=0; i<nB-1; i++){

            //funkcje kształtu liczone dla sciany
            for(int p=0;p<numberOfPoints;p++){
                N[p][0] = 0.25*(1-gauss.points[p])*(1-(-1));
                N[p][1] = 0.25*(1+gauss.points[p])*(1-(-1));
                N[p][2] = 0.25*(1-gauss.points[p])*(1+(-1));
                N[p][3] = 0.25*(1+gauss.points[p])*(1+(-1));
            }

            //wzory na prezentacji
            for(int x=0;x<numberOfPoints;x++){
                double det = dx/2; //x1 - x0
                for(int y=0;y<4;y++){
                    for(int z=0;z<4;z++){

                        elements[i].Hbc[y][z] += gauss.weights[x] * alpha *N[x][y] *N[x][z] * det;
                    }
                    elements[i].P[y] += alpha * gauss.weights[x] * N[x][y] * tempOfSurroundings*det;
                }

            }
        }

        //sciana 1, prawa
        for(int i=(nB-2); i<nE; i+=nB-1){

            for(int p=0;p<numberOfPoints;p++){
                N[p][0] = 0.25*(1-1)*(1-gauss.points[p]);
                N[p][1] = 0.25*(1+1)*(1-gauss.points[p]);
                N[p][2] = 0.25*(1-1)*(1+gauss.points[p]);
                N[p][3] = 0.25*(1+1)*(1+gauss.points[p]);
            }

            for(int x=0;x<numberOfPoints;x++){
                double det =dy/2;
                for(int y=0;y<4;y++){
                    for(int z=0;z<4;z++){

                        elements[i].Hbc[y][z] += gauss.weights[x] * alpha *N[x][y] *N[x][z] * det;
                    }
                    elements[i].P[y] += alpha * gauss.weights[x] * N[x][y] * tempOfSurroundings*det;
                }
            }
        }

        //sciana 2, górna
        for(int i=nE-1; i>nE-nB; i--){

            for(int p=0;p<numberOfPoints;p++){
                N[p][0] = 0.25*(1-gauss.points[p])*(1-(1));
                N[p][1] = 0.25*(1+gauss.points[p])*(1-(1));
                N[p][2] = 0.25*(1-gauss.points[p])*(1+(1));
                N[p][3] = 0.25*(1+gauss.points[p])*(1+(1));
            }

            for(int x=0;x<numberOfPoints;x++){
                double det = dx/2;
                for(int y=0;y<4;y++){
                    for(int z=0;z<4;z++){
                        elements[i].Hbc[y][z] += gauss.weights[x] * alpha *N[x][y] *N[x][z] * det;
                    }
                    elements[i].P[y] += alpha * gauss.weights[x] * N[x][y] * tempOfSurroundings*det;
                }
            }
        }

        //sciana 3, lewa
        for(int i=nE-(nB-1); i>=0; i-=nB-1){

            for(int p=0;p<numberOfPoints;p++){
                N[p][0] = 0.25*(1-(-1))*(1-gauss.points[p]);
                N[p][1] = 0.25*(1+(-1))*(1-gauss.points[p]);
                N[p][2] = 0.25*(1-(-1))*(1+gauss.points[p]);
                N[p][3] = 0.25*(1+(-1))*(1+gauss.points[p]);
            }


            for(int x=0;x<numberOfPoints;x++){
                double det = dy/2;

                for(int y=0;y<4;y++){
                    for(int z=0;z<4;z++){

                        elements[i].Hbc[y][z] += gauss.weights[x] * alpha *N[x][y] *N[x][z] * det;
                    }
                    elements[i].P[y] += alpha * gauss.weights[x] * N[x][y] * tempOfSurroundings*det;
                }
            }
        }
    }

    public void countGlobalHandPwithC(double dt, double[] T0){

        for(int x=0;x<nN;x++){
            for(int y=0;y<nN;y++){
                globalH[x][y] += (globalC[x][y] / dt);

                globalP[x] += (globalC[x][y] * T0[y]) /dt;
            }
        }
    }

    public void printGlobalHandC(){
        String temp=new String();

        for(int i=0;i<nN;i++){
            temp+=Arrays.toString(globalH[i])+"\n";
        }

        temp +="}, globalC={\n";

        for(int i=0;i<nN;i++){
            temp+=Arrays.toString(globalC[i])+"\n";
        }

        System.out.println(temp);
    }

    public void printGlobalP(){
        System.out.println(Arrays.toString(globalP));
    }

    public double[] countT(){

        double[] T = new double[nN];


        SimpleMatrix simpleMatrix1 = new SimpleMatrix(globalH);
        SimpleMatrix simpleMatrix2 = new SimpleMatrix(nN,1);

        for(int i=0;i<nN;i++){
            simpleMatrix2.setRow(i,0,globalP[i]);
        }


        Equation eq = new Equation();
        eq.alias(simpleMatrix1,"H",simpleMatrix2,"P");
        eq.process("t=P/H");
        //System.out.println((eq.lookupSimple("t")).get(0));

        try{
            SimpleMatrix solution = simpleMatrix1.solve(simpleMatrix2);
            for(int i=0;i<nN;i++){
                T[i] = solution.get(i,0);
            }
        }catch (Exception e){
            e.printStackTrace();
        }

        return T;


    }

    @Override
    public String toString() {

        System.out.println("Nodes:");
        for(Node temp:nodes){
            System.out.println("\t"+temp+"\n");
        }

        System.out.println("Elements:   ");
        for (Element temp: elements){
            System.out.println("\t"+temp+"\n");
        }

        String temp=new String();

        for(int i=0;i<nN;i++){
            temp+=Arrays.toString(globalH[i])+"\n";
        }

        return "Grid{" +
                "H=" + H +
                ", B=" + B +
                ", nH=" + nH +
                ", nB=" + nB +
                ", nN=" + nN +
                ", nE=" + nE +
                ", dx=" + dx +
                ", dy=" + dy +
                ", globalH={\n"+temp+"\n}"+
                '}';
    }
}

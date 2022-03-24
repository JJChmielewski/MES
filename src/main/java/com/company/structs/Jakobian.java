package com.company.structs;

import java.util.Arrays;

public class Jakobian {

    double[][] jakobian;
    double[][] jakobianInverse;

    double det;

    public Jakobian(Element element, Element4_2D element4_2D, int numberOfPoint){
        jakobian = new double[2][2];
        jakobianInverse = new double[2][2];

        Node[] nodes = element.nodes;
        //Node[] nodes = new Node[]{new Node(0.0,0.0), new Node(0.025,0.0), new Node(0.0, 0.025),new Node(0.025, 0.025)};

        jakobian[0][0] = (nodes[0].x * element4_2D.dNdKsi[numberOfPoint][0]) + (nodes[1].x * element4_2D.dNdKsi[numberOfPoint][1])+
                (nodes[2].x * element4_2D.dNdKsi[numberOfPoint][2]) + (nodes[3].x * element4_2D.dNdKsi[numberOfPoint][3]);

        jakobian[0][1] = (nodes[0].y * element4_2D.dNdKsi[numberOfPoint][0]) + (nodes[1].y * element4_2D.dNdKsi[numberOfPoint][1])+
                (nodes[2].y * element4_2D.dNdKsi[numberOfPoint][2]) + (nodes[3].y * element4_2D.dNdKsi[numberOfPoint][3]);

        jakobian[1][0] = (nodes[0].x * element4_2D.dNdEta[numberOfPoint][0]) + (nodes[1].x * element4_2D.dNdEta[numberOfPoint][1])+
                (nodes[2].x * element4_2D.dNdEta[numberOfPoint][2]) + (nodes[3].x * element4_2D.dNdEta[numberOfPoint][3]);

        jakobian[1][1] = (nodes[0].y * element4_2D.dNdEta[numberOfPoint][0]) + (nodes[1].y * element4_2D.dNdEta[numberOfPoint][1])+
                (nodes[2].y * element4_2D.dNdEta[numberOfPoint][2]) + (nodes[3].y * element4_2D.dNdEta[numberOfPoint][3]);

        jakobianInverse[0][0] = jakobian[1][1];
        jakobianInverse[0][1] = -1 * jakobian[0][1];
        jakobianInverse[1][0] = -1 * jakobian[1][0];
        jakobianInverse[1][1] = jakobian[0][0];

        det = (jakobian[0][0] * jakobian[1][1])-(jakobian[0][1] * jakobian[1][0]);
    }

    public Jakobian(int numberOfElement, int numberOfPoint, Element4_2D element4_2D, Grid grid){

        jakobian = new double[2][2];
        jakobianInverse = new double[2][2];

        Node[] nodes = grid.elements[numberOfElement].nodes;
        //Node[] nodes = new Node[]{new Node(0.0,0.0), new Node(0.025,0.0), new Node(0.0, 0.025),new Node(0.025, 0.025)};

        jakobian[0][0] = (nodes[0].x * element4_2D.dNdKsi[numberOfPoint][0]) + (nodes[1].x * element4_2D.dNdKsi[numberOfPoint][1])+
                (nodes[2].x * element4_2D.dNdKsi[numberOfPoint][2]) + (nodes[3].x * element4_2D.dNdKsi[numberOfPoint][3]);

        jakobian[0][1] = (nodes[0].y * element4_2D.dNdKsi[numberOfPoint][0]) + (nodes[1].y * element4_2D.dNdKsi[numberOfPoint][1])+
                (nodes[2].y * element4_2D.dNdKsi[numberOfPoint][2]) + (nodes[3].y * element4_2D.dNdKsi[numberOfPoint][3]);

        jakobian[1][0] = (nodes[0].x * element4_2D.dNdEta[numberOfPoint][0]) + (nodes[1].x * element4_2D.dNdEta[numberOfPoint][1])+
                (nodes[2].x * element4_2D.dNdEta[numberOfPoint][2]) + (nodes[3].x * element4_2D.dNdEta[numberOfPoint][3]);

        jakobian[1][1] = (nodes[0].y * element4_2D.dNdEta[numberOfPoint][0]) + (nodes[1].y * element4_2D.dNdEta[numberOfPoint][1])+
                (nodes[2].y * element4_2D.dNdEta[numberOfPoint][2]) + (nodes[3].y * element4_2D.dNdEta[numberOfPoint][3]);

        jakobianInverse[0][0] = jakobian[1][1];
        jakobianInverse[0][1] = -1 * jakobian[0][1];
        jakobianInverse[1][0] = -1 * jakobian[1][0];
        jakobianInverse[1][1] = jakobian[0][0];

        det = (jakobian[0][0] * jakobian[1][1])-(jakobian[0][1] * jakobian[1][0]);


    }

    @Override
    public String toString() {
        return "Jakobian{" +
                "jakobian=[\n" + Arrays.toString(jakobian[0]) +
                ",\n" + Arrays.toString(jakobian[1]) +
                "],\njakobianInverse=[\n" + Arrays.toString(jakobianInverse[0]) +
                ",\n"+Arrays.toString(jakobianInverse[1])+
                "],\ndet=" + det +
                '}';
    }
}

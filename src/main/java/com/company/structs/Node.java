package com.company.structs;

public class Node {

    public double x,y;
    public double temperature;
    public int numberOfNode;

    public Node() {
        this.x=0;
        this.y=0;
    }

    public Node(double x, double y) {
        this.x = x;
        this.y = y;
    }

    @Override
    public String toString() {
        return "Node{" +
                "x=" + x +
                ", y=" + y +
                ", numberOfNode=" + numberOfNode +
                '}';
    }
}

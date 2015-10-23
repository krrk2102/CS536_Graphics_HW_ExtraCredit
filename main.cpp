//
//  main.cpp
//  CG_hwEC
//
//  Created by Shangqi Wu on 14/11/19.
//  Copyright (c) 2014 Shangqi Wu. All rights reserved.
//

#include <iostream>
#include <stdlib.h>
#include <unistd.h>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <cmath>
#include <limits>

using namespace std;

// Claiming classes and structs will be used
class Point;
class Vertex;
class Line;
typedef struct Face_sturcture{
    int p1;
    int p2;
    int p3;
}Face;

// Claiming global variants will be used
int j=0, k=0, o=500, p=500;
double x=0, y=0, z=1, X=0, Y=0, Z=0, q=0, r=0, w=-1, Q=0, R=1, W=0, u=-0.7, v=-0.7, U=0.7, V=0.7, F=0.6, B=-0.6;
string f = "./bound-lo-sphere.smf", output = "out.xpm";
bool P = false; // If it is true, this program uses parallel projection, otherwise it uses perspective projection.
//bool b = false;
bool b = true;
//bool P = true; // For debugging use, testing parallel projection feature.
vector<vector<char> > pixel;

// Claiming functions will be used
string setheader();
string setend();
void help();
int rnd(double arg);
void outfile(string output);
void optana(int argc, char * const argv[]);
void readfile(string input, vector<Vertex> *vecver, vector<Face> *vecf);
vector<Vertex> projection(vector<Vertex> vecver);
int cstest(Line argl, int a, int b, int c, int d);
Line csclip(Line argl, int a, int b, int c, int d);
vector<double> cross_product(vector<double> input_vec1, vector<double> input_vec2);
double dot_product(vector<double> input_vec1, vector<double> input_vec2);

//--------------------------------------------------------------------------------------------------
class Point {
private:
    double px;
    double py;
    void trans(double m, double n) { // Tanslation
        double matrix[3][3] = {1, 0, m, 0, 1, n, 0, 0, 1};
        double pt[3] = {0, 0, 0};
        for (int i = 0; i < 3; i++) {
            pt[i] += matrix[i][0] * px;
            pt[i] += matrix[i][1] * py;
            pt[i] += matrix[i][2];
        }
        px=pt[0]; py=pt[1];
    }
    void scale(double sx,double sy) { // Scaling
        double matrix[3][3] = {sx, 0, 0, 0, sy, 0, 0, 0, 1};
        double pt[3] = {0, 0, 0};
        for (int i = 0; i < 3  ; i++) {
            pt[i] += matrix[i][0] * px;
            pt[i] += matrix[i][1] * py;
            pt[i] += matrix[i][2];
        }
        px=pt[0]; py=pt[1];
    }
public:
    void set(double argx, double argy) {
        px = argx;
        py = argy;
    }
    int getx() {
        return rnd(px);
    }
    int gety() {
        return rnd(py);
    }
    void toviewport(){
        if (P == false) {
            double zproj = abs(z/(B-z));
            double sx = (double)(o-j) / (double)(2*zproj);
            double sy = (double)(p-k) / (double)(2*zproj);
            trans(zproj, zproj);
            scale(sx, sy);
            trans((double)j, (double)k);
        }else{
            double sx = (double)(o-j) / (double)(2);
            double sy = (double)(p-k) / (double)(2);
            trans((double)1, (double)1);
            scale(sx, sy);
            trans((double)j, (double)k);
        }
    }
    double x(){
        return px;
    }
    double y(){
        return py;
    }
};

//-------------------------------------------------------------------------------------------------
class Vertex{
private:
    double verx;
    double very;
    double verz;
    int status; // Indicating positon of the vertex by Cohen-Sutherland code.
public:
    void set(double a, double b, double c){
        verx = a;
        very = b;
        verz = c;
    }
    int stats(){
        return status;
    }
    double getx() {
        return verx;
    }
    double gety() {
        return very;
    }
    double getz() {
        return verz;
    }
    void trivial_test(){
        status = 0b00000000;
        if (P == true) { // Trivial test for parallel projection.
            if (very > 1) { // Using 6-bit code for Cohen-Sutherland pretest
                status |= 0b000001;
            }else if (very < -1){
                status |= 0b000010;
            }
            if (verx > 1) {
                status |= 0b000100;
            }else if (verx < -1){
                status |= 0b001000;
            }
            if (verz < -1) {
                status |= 0b010000;
            }else if (verz >0){
                status |= 0b100000;
            }
        }else{
            double zmin = (z-F) / (B-z);
            if (very > -verz) {
                status |= 0b000001;
            }else if (very < verz){
                status |= 0b000010;
            }
            if (verx > -verz) {
                status |= 0b000100;
            }else if (verx < verz){
                status |= 0b001000;
            }
            if (verz < -1) {
                status |= 0b010000;
            }else if (verz > zmin){
                status |= 0b100000;
            }
        }
    }
    void matrix_multiply(double a[4][4]){
        double pt[4] = {0, 0, 0, 0};
        double temp = 0;
        for (int i = 0; i < 4; i++) {
            temp = temp + (a[i][0] * verx + a[i][1] * very);
            temp = temp + (a[i][2] * verz + a[i][3]);
            pt[i] = temp;
            temp = 0;
        }
        verx = pt[0]; very = pt[1]; verz = pt[2];
    }
    Point toProjection(){
        Point a;
        if (P == true) { // If implementing parallel projection, directly delete z coordinates.
            a.set(verx, very);
        }else{ // Implementing perspective projection to the desired plane.
            double zmin = (z-F) / (B-z);
            double zproj = z / (B-z);
            if (verz <= zmin) {
                double xtmp = verx * zproj / verz;
                double ytmp = very * zproj / verz;
                a.set(xtmp, ytmp);
            }else{ // Still generate a boundry point within the world window, to ensure continuity of the image.
                double xtmp = verx * zproj / zmin;
                double ytmp = very * zproj / zmin;
                a.set(xtmp, ytmp);
            }
        }
        return a;
    }
    vector<double> operator-(Vertex new_v) {
        vector<double> new_vec(3);
        new_vec[0] = verx - new_v.getx();
        new_vec[1] = very - new_v.gety();
        new_vec[2] = verz - new_v.getz();
        return new_vec;
    }
};

//--------------------------------------------------------------------------------------------------
class Line { // present the line by: y = slope*x + d
private:
    Point start;
    Point end;
    int xmax;
    int xmaxy;
    int xmin;
    int xminy;
    int ymax;
    int ymin;
    bool sd_exist;
    double slope;
    double d1;
    bool in;
    void cal() {
        if (start.getx() != end.getx()) {
            if (start.getx()>=end.getx()) {
                xmax = start.getx();
                xmaxy = start.gety();
                xmin = end.getx();
                xminy = end.gety();
            }else {
                xmax = end.getx();
                xmaxy = end.gety();
                xmin = start.getx();
                xminy = start.gety();
            }
            if (start.gety()>=end.gety()) { // set ymax and ymin
                ymax = start.gety();
                ymin = end.gety();
            }else {
                ymax = end.gety();
                ymin = start.gety();
            }
            slope = (start.y() - end.y()) / (start.x() - end.x());
            d1 = start.y() - (slope * start.x());
            sd_exist = true;
        }else {
            if (start.gety()>=end.gety()) {
                ymax = start.gety();
                ymin = end.gety();
                xmaxy = start.gety();
                xminy = end.gety();
            }else {
                ymax = end.gety();
                ymin = start.gety();
                xmaxy = end.gety();
                xminy = start.gety();
            }
            xmax = start.getx();
            xmin = end.getx();
            slope = numeric_limits<double>::max();
            d1 = numeric_limits<double>::max();
            sd_exist = false;
        }
    }
public:
    void setall(double argx1, double argy1, double argx2, double argy2) {
        start.set(argx1, argy1);
        end.set(argx2, argy2);
        cal();
    }
    void setp(Point a, Point b) {
        start = a;
        end = b;
        cal();
    }
    bool sd() {
        return sd_exist;
    }
    Point gets() {
        return start;
    }
    Point gete() {
        return end;
    }
    Point cal_inter_wlr(float x) { // return the intersaction point of x=x
        Point tmp;
        if (sd_exist) {
            double y = slope*x +d1;
            tmp.set(x, y);
        }else {
            tmp.set(numeric_limits<double>::max(), numeric_limits<double>::max());
        }
        return tmp;
    }
    Point cal_inter_wtb(float y) { // return the intersaction point of y=y
        Point tmp;
        if (sd_exist && slope!=0) {
            double x = (y - d1) / slope;
            tmp.set(x, y);
        }else if (slope == 0) {
            tmp.set(numeric_limits<int>::max(), y);
        }else {
            tmp.set(start.x(), y);
        }
        return tmp;
    }
    void setin(bool ifin) {
        in = ifin;
    }
    bool getin() {
        return in;
    }
    int getxmin() {
        return xmin;
    }
    int getxmax() {
        return xmax;
    }
    int getymin() {
        return ymin;
    }
    int getymax() {
        return ymax;
    }
    int getxmaxy() {
        return xmaxy;
    }
    int getxminy() {
        return xminy;
    }
    void showall() { // Calculate all the points of the line.
        if (sd_exist == true) {
            // Fllowing codes are of Bresenham Algorithm, presented in L-02_Lines.pdf. Codes are modifiied for this cpp file.
            int dx, dy, D, x, y;
            dx = xmax - xmin;
            dy = ymax - ymin;
            if (0<slope && slope<1) {
                D = 2*dy - dx;
                y = ymin;
                for (x = xmin; x <= xmax; x++) {
                    pixel[y][x] = '+';
                    if (D <= 0) {
                        D += 2*dy;
                    }else {
                        D += 2*(dy - dx);
                        y++;
                    }
                }
            }else if (slope > 1) {
                D = 2*dx - dy;
                x = xmin;
                for (y = ymin; y <= ymax; y++) {
                    pixel[y][x] = '+';
                    if (D <= 0) {
                        D += 2*dx;
                    }else {
                        D += 2*(dx - dy);
                        x++;
                    }
                }
            }else if (-1<slope && slope<0) {
                D = 2*dy - dx;
                y = ymax;
                for (x = xmin; x <= xmax; x++) {
                    pixel[y][x] = '+';
                    if (D <= 0) {
                        D += 2*dy;
                    }else {
                        D += 2*(dy - dx);
                        y--;
                    }
                }
            }else if (slope == 1) {
                y = ymin;
                for (x = xmin; x <= xmax; x++) {
                    pixel[y][x]='+';
                    y++;
                }
            }else if (slope == -1) {
                y = ymax;
                for (x = xmin; x <= xmax; x++) {
                    pixel[y][x] = '+';
                    y--;
                }
            }else if (slope == 0) {
                y = ymin;
                for (x = xmin; x <= xmax; x++) {
                    pixel[y][x] = '+';
                }
            }
            else { // i.e., slope<-1
                D = 2*dx - abs(dy);
                x = xmin;
                for (y = ymax; y >= ymin; y--) {
                    pixel[y][x] = '+';
                    if (D <= 0) {
                        D += 2*dx;
                    }else {
                        D += 2*(dx - dy);
                        x++;
                    }
                }
            }
        }else if (sd_exist == false) { // for vertical lines
            int x = xmin;
            for (int y = ymin; y <= ymax; y++) {
                pixel[y][x] = '+';
            }
        }
    }
};

// Main Function------------------------------------------------------------------------------------
int main(int argc, char * argv[]) {
    // Prepare pixels and input instructions.
    int height = 501;
    int width = 501;
    pixel.resize(height);
    for (int i = 0; i < height; i++) {
        pixel[i].resize(width);
        for (int j = 0; j < width; j++) {
            pixel[i][j] = '-';
        }
    }
    // Analyzing input instructions.
    optana(argc, argv);
    Vertex viewpoint;
    viewpoint.set(X, Y, Z);
    f = "/Users/wushangqi/cghw4/bound-lo-sphere.smf"; // dir only for debug*****************************************
    vector<Vertex> vecver;
    vector<Face> vecf;
    // Read the file.
    readfile(f, &vecver, &vecf);
    // Processing 3D Projection.
    vecver = projection(vecver);
    vector<Point> vecp;
    int psize = (int)vecver.size();
    int fsize = (int)vecf.size();
    // Translate #d Vertex into Points in the view plane.
    for (int i = 0; i < psize; i++) {
        vecver[i].trivial_test();
        vecp.push_back(vecver[i].toProjection());
        vecp[i].toviewport();
    }
    // Prepare to draw edges of each face.
    vector<Line> vecl;
    for (int i = 0; i < fsize; i++) { // Analysing trivial test result first.
        int s1 = vecver[vecf[i].p1].stats() & vecver[vecf[i].p2].stats();
        int s2 = vecver[vecf[i].p2].stats() & vecver[vecf[i].p3].stats();
        int s3 = vecver[vecf[i].p3].stats() & vecver[vecf[i].p1].stats();
        if (s1==0 || s2==0 || s3==0) { // If the face is completely outside the view volume, then its edges will not be drawed.
            if (b == true) { // If implement backface culliing, calculate its normalizing vector.
                vector<double> v1 = vecver[vecf[i].p2] - vecver[vecf[i].p1];
                vector<double> v2 = vecver[vecf[i].p3] - vecver[vecf[i].p1];
                // Since smf vertex are in order of counter-clockwise, calculated normal vector pointing towards screen.
                vector<double> norm = cross_product(v1, v2);
                double dir;
                if (P == false) {
                    vector<double> vp = viewpoint - vecver[vecf[i].p1];
                    // If dot peoduct of normal vector and vector pointing from face to view point is positive, this face is front face adn should not be clipped.
                    dir = dot_product(norm, vp);
                } else {
                    vector<double> vp(3);
                    vp[0]=0; vp[1]=0; vp[2] = Z - vecver[vecf[i].p1].getz();
                    dir = dot_product(norm, vp);
                }
                if (dir > 0) {
                    Line buff;
                    buff.setp(vecp[vecf[i].p1], vecp[vecf[i].p2]);
                    vecl.push_back(buff);
                    buff.setp(vecp[vecf[i].p1], vecp[vecf[i].p3]);
                    vecl.push_back(buff);
                    buff.setp(vecp[vecf[i].p2], vecp[vecf[i].p3]);
                    vecl.push_back(buff);
                }
            }else { // If do not implement backface culling, draw every edge of each face.
                Line buff;
                buff.setp(vecp[vecf[i].p1], vecp[vecf[i].p2]);
                vecl.push_back(buff);
                buff.setp(vecp[vecf[i].p1], vecp[vecf[i].p3]);
                vecl.push_back(buff);
                buff.setp(vecp[vecf[i].p2], vecp[vecf[i].p3]);
                vecl.push_back(buff);
            }
        }
    }
    // Clip and draw lines in the viewport.
    int lsize = (int)vecl.size();
    for (int i = 0; i < lsize; i++) {
        vecl[i] = csclip(vecl[i], j, k, o, p);
        if (vecl[i].getin() == true) {
            vecl[i].showall();
        }
    }
    // Write output xpm file.
    output = "/Users/wushangqi/out.xpm"; // dir only for debug******************************************************
    outfile(output);
    // Display output xpm image automatically.
    string shell = "display " + output;
    system(shell.c_str());
    return 0;
}

//---------------------------------------------------------------------------------------------------
double dot_product(vector<double> input_vec1, vector<double> input_vec2) {
    double dot = input_vec1[0]*input_vec2[0] + input_vec1[1]*input_vec2[1] + input_vec1[2]*input_vec2[2];
    return dot;
}

//---------------------------------------------------------------------------------------------------
vector<double> cross_product(vector<double> input_vec1, vector<double> input_vec2) {
    vector<double> new_vec(3);
    new_vec[0] = input_vec1[1]*input_vec2[2] - input_vec1[2]*input_vec2[1];
    new_vec[1] = input_vec1[2]*input_vec2[0] - input_vec1[0]*input_vec2[2];
    new_vec[2] = input_vec1[0]*input_vec2[1] - input_vec1[1]*input_vec2[0];
    return new_vec;
}

//---------------------------------------------------------------------------------------------------
vector<Vertex> projection(vector<Vertex> vecver){
    // Prepare to process vertex.
    int versize = (int)vecver.size();
    vector<double> rz(3); // VPN / |VPN|
    double rzsqrt = sqrt(q*q + r*r + w*w);
    rz[0] = q / rzsqrt;
    rz[1] = r / rzsqrt;
    rz[2] = w / rzsqrt;
    vector<double> rx(3); // VUP x Rz first
    rx[0] = R*rz[2] - W*rz[1];
    rx[1] = W*rz[0] - Q*rz[2];
    rx[2] = Q*rz[1] - R*rz[0];
    double rxsqrt = sqrt(rx[0]*rx[0] + rx[1]*rx[1] + rx[2]*rx[2]); // then normalizing it
    rx[0] = rx[0] / rxsqrt;
    rx[1] = rx[1] / rxsqrt;
    rx[2] = rx[2] / rxsqrt;
    vector<double> ry(3); // Implementing Rz x Rx
    ry[0] = rz[1]*rx[2] - rz[2]*rx[1];
    ry[1] = rz[2]*rx[0] - rz[0]*rx[2];
    ry[2] = rz[0]*rx[1] - rz[1]*rx[0];
    double shx = (0.5*(U+u)-x) / z;
    double shy = (0.5*(V+v)-y) / z;
    // Starting to process 3D projection.
    if (P == true) { // Implementing Parallel Projection
        double tpar1 = -(U+u) / 2.0;
        double tpar2 = -(V+v) / 2.0;
        double spar1 = 2.0 / (U-u);
        double spar2 = 2.0 / (V-v);
        double spar3 = 1.0 / (F-B);
        double R[4][4] = {{rx[0], rx[1], rx[2], 0}, {ry[0], ry[1], ry[2], 0}, {rz[0], rz[1], rz[2], 0}, {0, 0, 0, 1}};
        double Tvrp[4][4] = {{1, 0, 0, -X}, {0, 1, 0, -Y}, {0, 0, 1, -Z}, {0, 0, 0, 1}};
        double SHpar[4][4] = {{1, 0, shx, 0}, {0, 1, shy, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
        double Spar[4][4] = {{spar1, 0, 0, 0}, {0, spar2, 0, 0}, {0, 0, spar3, 0}, {0, 0, 0, 1}};
        double Tpar[4][4] = {{1, 0, 0, tpar1}, {0, 1, 0, tpar2}, {0, 0, 1, -F}, {0, 0, 0, 1}};
        double tmp1[4][4] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        double tmp2[4][4] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        double tmp3[4][4] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        double Npar[4][4] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        double temp = 0;
        for (int i = 0; i < 4; i++) { // Pipeline Step 1: Translate VRP to the origin.
            for (int j = 0; j < 4; j++) { // Pipeline Step 2: Rotate so VPN becomes z, VUP becomes y, and u becomes x.
                temp = temp + (R[i][0] * Tvrp[0][j] + R[i][1] * Tvrp[1][j]);
                temp = temp + (R[i][2] * Tvrp[2][j] + R[i][3] * Tvrp[3][j]);
                tmp1[i][j] = temp;
                temp = 0;
            }
        }
        for (int i = 0; i < 4; i++) { // Pipeline Step 3: Shear to make direction of the projection become parallel to z.
            for (int j = 0; j < 4; j++) {
                temp = temp + (SHpar[i][0] * tmp1[0][j] + SHpar[i][1] * tmp1[1][j]);
                temp = temp + (SHpar[i][2] * tmp1[2][j] + SHpar[i][3] * tmp1[3][j]);
                tmp2[i][j] = temp;
                temp = 0;
            }
        }
        for (int i = 0; i < 4; i++) { // Pipeline Step 4: Translate and scale into a canonical view volume.
            for (int j = 0; j < 4; j++) {
                temp = temp + (Tpar[i][0] * tmp2[0][j] + Tpar[i][1] * tmp2[1][j]);
                temp = temp + (Tpar[i][2] * tmp2[2][j] + Tpar[i][3] * tmp2[3][j]);
                tmp3[i][j] = temp;
                temp = 0;
            }
        }
        for (int i = 0; i < 4; i++) {  // Implementing last scaling step to generate normalizing matrix.
            for (int j = 0; j < 4; j++) {
                temp = temp + (Spar[i][0] * tmp3[0][j] + Spar[i][1] * tmp3[1][j]);
                temp = temp + (Spar[i][2] * tmp3[2][j] + Spar[i][3] * tmp3[3][j]);
                Npar[i][j] = temp;
                temp = 0;
            }
        }
        for (int i = 0; i < versize; i++) { // Using this matrix to process coordinates of all vertex.
            vecver[i].matrix_multiply(Npar);
            vecver[i].trivial_test(); // Implemeting simple view volume rejection test for each vertex.
        }
    }else{ // Implementing Perspective Projection
        double sperx = 2.0*(-z)/((U-u)*(-z+B));
        double spery = 2.0*(-z)/((V-v)*(-z+B));
        double sperz = -1.0/(-z+B);
        double R[4][4] = {{rx[0], rx[1], rx[2], 0}, {ry[0], ry[1], ry[2], 0}, {rz[0], rz[1], rz[2], 0}, {0, 0, 0, 1}};
        double Tvrp[4][4] = {{1, 0, 0, -X}, {0, 1, 0, -Y}, {0, 0, 1, -Z}, {0, 0, 0, 1}};
        double SHper[4][4] = {{1, 0, shx, 0}, {0, 1, shy, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
        double Sper[4][4] = {{sperx, 0, 0, 0}, {0, spery, 0, 0}, {0, 0, sperz, 0}, {0, 0, 0, 1}};
        double Tprp[4][4] = {{1, 0, 0, -x}, {0, 1, 0, -y}, {0, 0, 1, -z}, {0, 0, 0, 1}};
        double tmp1[4][4] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        double tmp2[4][4] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        double tmp3[4][4] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        double Nper[4][4] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        double temp = 0;
        for (int i = 0; i < 4; i++) { // Pipeline Step 1: Translate VRP to the origin.
            for (int j = 0; j < 4; j++) {
                temp = temp + (R[i][0] * Tvrp[0][j] + R[i][1] * Tvrp[1][j]);
                temp = temp + (R[i][2] * Tvrp[2][j] + R[i][3] * Tvrp[3][j]);
                tmp1[i][j] = temp;
                temp = 0;
            }
        }
        for (int i = 0; i < 4; i++) { // Pipeline Step 2: Rotate so VPN becomes z, VUP becomes y, and u becomes x.
            for (int j = 0; j < 4; j++) {
                temp = temp + (Tprp[i][0] * tmp1[0][j] + Tprp[i][1] * tmp1[1][j]);
                temp = temp + (Tprp[i][2] * tmp1[2][j] + Tprp[i][3] * tmp1[3][j]);
                tmp2[i][j] = temp;
                temp = 0;
            }
        }
        for (int i = 0; i < 4; i++) { // Pipeline Step 3: Translate COP to origin.
            for (int j = 0; j < 4; j++) {
                temp = temp + (SHper[i][0] * tmp2[0][j] + SHper[i][1] * tmp2[1][j]);
                temp = temp + (SHper[i][2] * tmp2[2][j] + SHper[i][3] * tmp2[3][j]);
                tmp3[i][j] = temp;
                temp = 0;
            }
        }
        for (int i = 0; i < 4; i++) { // Pipeline Step 5: Scale into a canonical view volume for clippiing.
            for (int j = 0; j < 4; j++) {
                temp = temp + (Sper[i][0] * tmp3[0][j] + Sper[i][1] * tmp3[1][j]);
                temp = temp + (Sper[i][2] * tmp3[2][j] + Sper[i][3] * tmp3[3][j]);
                Nper[i][j] = temp;
                temp = 0;
            }
        } // Generated the normalizing matrix for all coordinates calculation of every vertex.
        for (int i = 0; i < versize; i++) { // Using this matrix to process coordinates of all vertex.
            vecver[i].matrix_multiply(Nper);
            vecver[i].trivial_test(); // Implementing simple view volume rejection test.
        }
    }
    return vecver;
}

//---------------------------------------------------------------------------------------------------
int cstest(Line argl, int a, int b, int c, int d){ // Return 1 or 2 or 3 if the line is completely visible. The simple test of C-S algorithm
    int xs = argl.getxmin();
    int ys = argl.getxminy();
    int xe = argl.getxmax();
    int ye = argl.getxmaxy();
    if (((b<=ys&&ys<=d) && (b<=ye&&ye<=d)) && ((xs<a) && (xe>c))){
        return 1; //Lines go from WL to WR.
    }else if ((a<=xs&&xs<=c) && (a<=xe&&xe<=c) && (ys<b && ye>d)){
        return 2; //Lines go from WB to WT.
    }else if ((a<=xs&&xs<=c) && (a<=xe&&xe<=c) && (b<=ys&&ys<=d) && (b<=ye&&ye<=d)){
        return 3; //Lines begin and end within the world window.
    }else return 0; //Lines cannot pass the simple exam.
}

//---------------------------------------------------------------------------------------------------
Line csclip(Line argl, int a, int b, int c, int d){ //Futher test and clip the line by using Cohen-Sutherland algorithm.
    int xleft = argl.getxmin();
    int xlefty = argl.getxminy();
    int xright = argl.getxmax();
    int xrighty = argl.getxmaxy();
    int ya = argl.cal_inter_wlr((float)a).gety(); // intersaction with WL(x=a)
    int yc = argl.cal_inter_wlr((float)c).gety(); // intersaction with WR(x=c)
    int xb = argl.cal_inter_wtb((float)b).getx(); // intersaction with WB(y=b)
    int xd = argl.cal_inter_wtb((float)d).getx(); // intersaction with WT(y=d)
    int flag = cstest(argl, a, b, c, d);
    if (flag == 0){ //If the line has not passed the simple test, take it to do the following complex test.
        if (xleft < a && xright > c) {
            if ((b<=ya&&ya<=d) && (b<=yc&&yc<=d)){
                argl.setall(a, ya, c, yc);
                argl.setin(true);
            }else if ((b<=ya&&ya<=d) && yc<b){
                argl.setall(a, ya, xb, b);
                argl.setin(true);
            }else if ((b<=ya&&ya<=d) && yc>d){
                argl.setall(a, ya, xd, d);
                argl.setin(true);
            }else if (ya<b && (b<=yc&&yc<=d)){
                argl.setall(xb, b, c, yc);
                argl.setin(true);
            }else if (ya>d && (b<=yc&yc<=d)){
                argl.setall(xd, d, c, yc);
                argl.setin(true);
            }else if (ya<b && yc>d){
                argl.setall(xb, b, xd, d);
                argl.setin(true);
            }else if (ya>d && yc<b){
                argl.setall(xd, d, xb, b);
                argl.setin(true);
            }else argl.setin(false);
        }else if ((a<=xleft && xleft<=c) && xright>c){
            if ((b<=xlefty&&xlefty<=d) && (b<=yc&&yc<=d)){
                argl.setall(xleft, xlefty, c, yc);
                argl.setin(true);
            }else if ((b<=xlefty&xlefty<=d) && yc<b){
                argl.setall(xleft, xlefty, xb, b);
                argl.setin(true);
            }else if ((b<=xlefty&&xlefty<=d) && yc>d){
                argl.setall(xleft, xlefty, xd, d);
                argl.setin(true);
            }else if (xlefty<b && (b<=yc&&yc<=d)){
                argl.setall(xb, b, c, yc);
                argl.setin(true);
            }else if (xlefty<b && yc>d){
                argl.setall(xb, b, xd, d);
                argl.setin(true);
            }else if (xlefty>d && (b<=yc&&yc<=d)){
                argl.setall(xd, d, c, yc);
                argl.setin(true);
            }else if (xlefty>d && yc<b){
                argl.setall(xd, d, xb, b);
                argl.setin(true);
            }else argl.setin(false);
        }else if (xleft<a && (a<=xright && xright<=c)){
            if ((b<=ya&&ya<=d) && (b<=xrighty&&xrighty<=d)){
                argl.setall(a, ya, xright, xrighty);
                argl.setin(true);
            }else if (ya<b && (b<=xrighty&&xrighty<=d)){
                argl.setall(xb, b, xright, xrighty);
                argl.setin(true);
            }else if (ya>d && (b<=xrighty&&xrighty<=d)){
                argl.setall(xd, d, xright, xrighty);
                argl.setin(true);
            }else if ((b<=ya&&ya<=d) && xrighty<b){
                argl.setall(a, ya, xb, b);
                argl.setin(true);
            }else if ((b<=ya&&ya<=d) && xrighty>d){
                argl.setall(a, ya, xd, d);
                argl.setin(true);
            }else if (ya<b && xrighty>d){
                argl.setall(xb, b, xd, d);
                argl.setin(true);
            }else if (ya>d &&xrighty<b){
                argl.setall(xd, d, xb, b);
                argl.setin(true);
            }else argl.setin(false);
        }else if ((a<=xleft&&xleft<=c) && (a<=xright&&xright<=c)){
            if (xlefty<b && (b<=xrighty&&xrighty<=d)){
                argl.setall(xb, b, xright, xrighty);
                argl.setin(true);
            }else if (xlefty<b && xrighty>d){
                argl.setall(xb, b, xd, d);
                argl.setin(true);
            }else if (xlefty>d && xrighty<b){
                argl.setall(xd, d, xb, b);
                argl.setin(true);
            }else if (xlefty>d && (b<=xrighty&&xrighty<=d)){
                argl.setall(xd, d, xright, xrighty);
                argl.setin(true);
            }else if ((b<=xlefty&&xlefty<=d) && xrighty<b){
                argl.setall(xleft, xlefty, xb, b);
                argl.setin(true);
            }else if ((b<=xlefty&&xlefty<=d) && xrighty>d){
                argl.setall(xleft, xlefty, xd, d);
                argl.setin(true);
            }else if ((b<=xlefty&&xlefty<=d) && (b<=xrighty&&xrighty<=d)){
                argl.setin(true);
            }else argl.setin(false);
        }else argl.setin(false);
    }else if (flag == 1){
        argl.setall(a, ya, c, yc);
        argl.setin(true);
    }else if (flag == 2){
        argl.setall(xb, b, xd, d);
        argl.setin(true);
    }else if (flag == 3){ // Lines totally within the window do not need to be clipped.
        argl.setin(true);
    }else argl.setin(false); // Mark and ignore all lines out the world window.
    return argl;
}

//---------------------------------------------------------------------------------------------------
void readfile(string input, vector<Vertex> *vecver, vector<Face> *vecf){
    // Start to read file.
    ifstream infile(f.c_str());
    if (!infile) {
        cout<<"Cannot open your file, please check your input path."<<endl;
        abort();
    }
    bool storev = false;
    bool storef = false;
    string str;
    vector<double> buff;
    while (infile) {
        infile>>str;
        if (str.compare("v") == 0) {
            storev = true;
        }else if (str.compare("f") == 0){
            storef = true;
        }else if (storev == true){
            buff.push_back(atof(str.c_str()));
            if (buff.size() == 3) {
                Vertex ver;
                ver.set(buff[0], buff[1], buff[2]);
                vecver->push_back(ver);
                storev = false;
                buff.clear();
            }
        }else if (storef == true){
            buff.push_back(atoi(str.c_str()));
            if (buff.size() == 3) {
                Face fac;
                fac.p1 = (int)buff[0] - 1;
                fac.p2 = (int)buff[1] - 1;
                fac.p3 = (int)buff[2] - 1;
                vecf->push_back(fac);
                storef = false;
                buff.clear();
            }
        }
    }
    infile.close();
}

//---------------------------------------------------------------------------------------------------
void outfile(string output){
    ofstream out(output.c_str());
    string line = "";
    if (!out) {
        cout<<"Cannot write an output file, please check your output path."<<endl;
    }
    out<<setheader()<<endl;
    int height = 501;
    int width = 501;
    //cout<<setheader()<<endl;
    for (int i = height - 1; i >= 0; i--) {
        for (int j = 0; j< width; j++) {
            line += pixel[i][j];
        }
        line = "\"" + line + "\"";
        if (i != 0) {
            line = line + ",";
        }
        out<<line<<endl;
        //cout<<line<<endl;
        line.clear();
    }
    out<<setend()<<endl;
    //cout<<setend()<<endl;
    out.close();
    
}

//---------------------------------------------------------------------------------------------------
string setheader() {
    stringstream tmp;
    int w = 501;
    int h = 501;
    tmp<<w;
    string intw;
    tmp>>intw;
    tmp.clear();
    tmp<<h;
    string inth;
    tmp>>inth;
    string str = "/* XPM */\nstatic char *CG_hw4[] = {\n/* width height num_colors chars_per_pixel */\n\"" + intw + " " + inth +" 2 1\",\n/* colors */\n\"- c #ffffff\",\n\"+ c #000000\",\n/* pixels */";
    return str;
}

//---------------------------------------------------------------------------------------------------
string setend() {
    string str = "};";
    return str;
}

//---------------------------------------------------------------------------------------------------
int rnd(double arg){
    if (arg >= 0) {
        return (int)(arg + 0.5);
    }else{
        return (int)(arg - 0.5);
    }
}

//---------------------------------------------------------------------------------------------------
void optana(int argc, char * const argv[]){
    // analyze input option and set default options
    int opt;
    while ((opt = getopt(argc, argv, "f:j:k:o:p:x:y:z:X:Y:Z:q:r:w:Q:R:W:u:v:U:V:F:B:Pbh"))!= -1) {
        switch (opt) {
            case 'f':
                f = optarg;
                break;
            case 'j':
                j = atoi(optarg);
                break;
            case 'k':
                k = atoi(optarg);
                break;
            case 'o':
                o = atoi(optarg);
                break;
            case 'p':
                p = atoi(optarg);
                break;
            case 'x':
                x = atof(optarg);
                break;
            case 'y':
                y = atof(optarg);
                break;
            case 'z':
                z = atof(optarg);
                break;
            case 'X':
                X  = atof(optarg);
                break;
            case 'Y':
                Y = atof(optarg);
                break;
            case 'Z':
                Z = atof(optarg);
                break;
            case 'q':
                q = atof(optarg);
                break;
            case 'r':
                r = atof(optarg);
                break;
            case 'w':
                w = atof(optarg);
                break;
            case 'Q':
                Q = atof(optarg);
                break;
            case 'R':
                R = atof(optarg);
                break;
            case 'W':
                W = atof(optarg);
                break;
            case 'u':
                u = atof(optarg);
                break;
            case 'v':
                v = atof(optarg);
                break;
            case 'U':
                U = atof(optarg);
                break;
            case 'V':
                V = atof(optarg);
                break;
            case 'F':
                F = atof(optarg);
                break;
            case 'B':
                B = atof(optarg);
                break;
            case 'P':
                P = true;
                break;
            case 'b':
                b = true;
                break;
            case 'h':
                help();
                abort();
                break;
            default:cout<<"Your input commands may not be correct. Please enter -h for help."<<endl;
                abort();
                break;
        }
    }
    if (j > o) {
        int tmp;
        tmp = j;
        j = o;
        o = tmp;
    }
    if (k > p) {
        int tmp;
        tmp = k;
        k = p;
        p = tmp;
    }
    if ((o>500) || (p>500)) {
        cout<<"The viewport cannot be larger than xpm image."<<endl;
        abort();
    }
}

//-------------------------------------------------------------------
void help(){
    cout<<"This program implements z-buffering with near and far plane clipping. The output XPM image has up to 60 colors, which are determined by distance from screen. Deeper the vertex inside the screen, darker the pixel displays."<<endl;
    cout<<"Output XPM image will be saved in the same working path with name of out.xpm"<<endl;
    cout<<"Following are the features and input options available of this program."<<endl;
    cout<<"[-f] Followed by the name of the first SMF model. Its base color is \"ff0000\" (red)."<<endl;
    cout<<"[-g] Followed by the name of the second SMF model. Its base color is \"00ff00\" (green)."<<endl;
    cout<<"[-i] Followed by the name of the third SMF model. Its base color is \"0000ff\" (blue)."<<endl;
    cout<<"[-j] The next argument is an integer lower bound in the x dimension of the view port in screen coordinates"<<endl;
    cout<<"[-k] The next argument is an integer in lower bound in the y dimension of the view port in screen coordinates"<<endl;
    cout<<"The next argument is an integer in upper bound in the x dimension of the view port in screen coordinates"<<endl;
    cout<<"[-p] The next argument is an integer in upper bound in the y dimension of the view port in screen coordinates"<<endl;
    cout<<"[-x] floating point x of Projection Reference Point (PRP) in VRC coordinates"<<endl;
    cout<<"[-y] floating point y of Projection Reference Point (PRP) in VRC coordinates"<<endl;
    cout<<"[-z] floating point z of Projection Reference Point (PRP) in VRC coordinates"<<endl;
    cout<<"[-X] floating point x of View Reference Point (VRP) in world coordinates"<<endl;
    cout<<"[-Y] floating point y of View Reference Point (VRP) in world coordinates"<<endl;
    cout<<"[-Z] floating point z of View Reference Point (VRP) in world coordinates"<<endl;
    cout<<"[-q] floating point x of View Plane Normal vector (VPN) in world coordinates"<<endl;
    cout<<"[-r] floating point y of View Plane Normal vector (VPN) in world coordinates"<<endl;
    cout<<"[-w] floating point z of View Plane Normal vector (VPN) in world coordinates"<<endl;
    cout<<"[-Q] floating point x of View Up Vector (VUP) in world coordinates"<<endl;
    cout<<"[-R] floating point y of View Up Vector (VUP) in world coordinates"<<endl;
    cout<<"[-W] floating point z of View UP Vector (VUP) in world coordinates"<<endl;
    cout<<"[-u] floating point u min of the VRC window in VRC coordinates"<<endl;
    cout<<"[-v] floating point v min of the VRC window in VRC coordinates"<<endl;
    cout<<"[-U] floating point u max of the VRC window in VRC coordinates"<<endl;
    cout<<"[-V] floating point v max of the VRC window in VRC coordinates"<<endl;
    cout<<"[-P] Use parallel projection. If this flag is not present, use perspective projection."<<endl;
    cout<<"[-F] Followed by the floating point coordinate of the Front (Near) plane in VRC coordinates."<<endl;
    cout<<"[-B] Followed by the floating point coordinate of the Back (Far) plane in VRC coordinates."<<endl;
    cout<<"[-b] Usw backface culling, If this flag is not present, do not culling any faces inside view volume."<<endl;
    cout<<"The default values are: "<<endl;
    cout<<"./CG_hw5 -f bound-sprellpsd.smf -j 0 -k 0 -o 500 -p 500 -x 0.0 -y 0.0 -z 1.0 -X 0.0 -Y 0.0 -Z 0.0 -q 0.0 -r 0.0 -w -1.0 -Q 0.0 -R 1.0 -W 0.0 -u -0.7 -v -0.7 -U 0.7 -V 0.7 -F 0.6 -B -0.6 > out.xpm"<<endl;
}
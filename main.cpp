#include <iostream>
#include <Eigen/Dense>
#include "src/domain/QuadMesh.h"
#include "src/domain/FeSpace.h"
#include "src/problems/LaplaceProblem.h"
#include <math.h>

using namespace std;

double analyticalSolution (double x, double y){
    //return 0.;
    //return (y*(y-1.)*x*(x-1.));
    return (x-x*x)*exp(3*x)*sin(2*M_PI*y);
    //return sin(2*M_PI*x)*sin(2*M_PI*y);
}

double f (double x, double y){
    //return x*x*y*y;
    //return -2*(y*(y-1.)+x*(x-1.));
    //return x*x*y*y;
    return exp(3*x)*sin(2*M_PI*y)*(4*M_PI*M_PI*(x-x*x)+9*x*x+3*x-4);
    //return 8*M_PI*M_PI*sin(2*M_PI*x)*sin(2*M_PI*y);
    //return 1.;
}

int main() {
    LaplaceProblem problem("parameters", f);
    problem.solve();
    problem.printOnFile("test");
    problem.compareWithExactSolution(analyticalSolution);
    problem.computeEnergyError(analyticalSolution);

    return 0;
}
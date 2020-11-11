

#include "AdvectionDiffusionProblem.h"
#include "../operators/InteriorPenalityOperator.h"

#define EPS 1e-6

template <typename T>
std::vector<T> linspace(T a, T b, size_t N) {
    T h = (b - a) / static_cast<T>(N-1);
    std::vector<T> xs(N);
    typename std::vector<T>::iterator x;
    T val;
    for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
        *x = val;
    return xs;
}

AdvectionDiffusionProblem::AdvectionDiffusionProblem(FeSpace fespace,
                               double gamma,
                               double tau,
                               double diffCoeff,
                               double advCoeffX,
                               double advCoeffY,
                               std::function<double (double, double)> f): fespace(fespace),
                                                                          gamma(gamma),
                                                                          tau(tau),
                                                                          diffCoeff(diffCoeff),
                                                                          advCoeffX(advCoeffX),
                                                                          advCoeffY(advCoeffY),
                                                                          f(f){
    solved = false;
}

void AdvectionDiffusionProblem::solve() {
    SystemMatrix systemMatrix(fespace);

    DiffusionOperator diffusionOperator(fespace, 0.);
    diffusionOperator.updateSystemMatrix(systemMatrix);

    AdvectionOperator advectionOperator(fespace, 1., 1.);
    advectionOperator.updateSystemMatrix(systemMatrix);

    InteriorPenalityOperator interiorPenality(fespace);
    interiorPenality.updateSystemMatrix(systemMatrix);

    InteriorPenalitySymmetric interiorPenalitySymmetric(fespace, tau);
    interiorPenalitySymmetric.updateSystemMatrix(systemMatrix);

    StabilizerOperator stabilizerOperator(fespace, gamma);
    stabilizerOperator.updateSystemMatrix(systemMatrix);

    TransportEdgeOperator transportEdgeOperator(fespace, 1., 1.);
    transportEdgeOperator.updateSystemMatrix(systemMatrix);

    systemMatrix.buildSparseMatrix();
    //std::cout<<matrix<<std::endl;
    RightHandSideFunction fProj (fespace, f, 10);
    auto fCoeff =fProj.getSystemVector();
    auto matrix = systemMatrix.getMatrix();
    //std::cout<<fCoeff<<std::endl;
    std::cout<< "Non zeros "<< matrix.nonZeros()<<std::endl;
    std::cout<< "Size "<< static_cast<double>(matrix.cols())*matrix.rows()<<std::endl;
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double,Eigen::RowMajor> > solver;
    solver.compute(matrix);
    solutionCoeff = solver.solve(fCoeff);
    solved = true;
}

void AdvectionDiffusionProblem::printOnScreen() {
    quadMesh mesh = fespace.getMesh();
    std::pair<double, double> xDomain = mesh.getDomainX();
    std::pair<double, double> yDomain = mesh.getDomainY();

    std::vector<double> xCoord = linspace(xDomain.first+ EPS, xDomain.second - EPS, 20);
    std::vector<double> yCoord = linspace(yDomain.first+ EPS, yDomain.second - EPS, 20);

    SolutionEvaluation eval(fespace);
    Eigen::MatrixXd evaluation = eval.evaluate(solutionCoeff, xCoord, yCoord);
    std::cout << evaluation << std::endl;

}

double AdvectionDiffusionProblem::compareWithExactSolution(std::function<double(double, double)> exactSolution) {
    Error error (exactSolution, solutionCoeff);
    double l2error = error.computeL2Error(fespace);
    std::cout<< "L2 error: " << l2error << std::endl;
    return l2error;
}

void AdvectionDiffusionProblem::printOnFile(std::string fileName) {
    quadMesh mesh = fespace.getMesh();
    std::pair<double, double> xDomain = mesh.getDomainX();
    std::pair<double, double> yDomain = mesh.getDomainY();

    unsigned int numPoints = 100;

    std::vector<double> xCoord = linspace(xDomain.first + EPS, xDomain.second - EPS, numPoints);
    std::vector<double> yCoord = linspace(yDomain.first + EPS, yDomain.second - EPS, numPoints);

    SolutionEvaluation eval(fespace);
    Eigen::MatrixXd evaluation = eval.evaluate(solutionCoeff, xCoord, yCoord);

    unsigned int numCols = numPoints;

    /*
    std::ofstream file("../output/" + fileName + ".txt");
    file << numCols << " ";
    for (int i = 0; i < numPoints; ++i) {
        file << xCoord[i] << " ";
    }
    file << std::endl;

    for (int i = 0; i < numPoints; ++i) {
        for (int j = 0; j < numPoints; ++j) {
            if (j == 0) {
                file << yCoord[i] << " ";
            }
            file << evaluation(i, j) << " ";
        }
        file << std::endl;
    }
    file << std::endl;

    file.close();
    */

    std::ofstream file2("../output/" + fileName + ".dat");
    //file << numCols << " ";
    for (int k = 0; k < numPoints; ++k) {
    for (int i = 0; i < numPoints; ++i) {
        file2 << xCoord[k] << " ";
        file2 << yCoord[i] << " ";

            file2 << evaluation(k, i) << " ";

        file2 << std::endl;
    }
    file2 << std::endl;
    }
    file2 << std::endl;

    file2.close();


}

AdvectionDiffusionProblem::AdvectionDiffusionProblem(std::string fileName,std::function<double (double, double)> f):
    f(f)
    {
    std::vector<double> parameters;
    std::string param_name;
    double param_value;
    std::ifstream file("../input/" + fileName + ".txt");
    while ( file >> param_name >> param_value )
    {
        parameters.push_back(param_value);
    }
    std::pair<double,double> domainX = std::make_pair(parameters[0],parameters[1]);
    std::pair<double,double> domainY = std::make_pair(parameters[2],parameters[3]);
    int nX= static_cast<int>(parameters[4]);
    int nY=static_cast<int>(parameters[5]);
    int maxDegree = static_cast<int>(parameters[6]);
    int localDegree =static_cast<int>(parameters[7]);
    diffCoeff = parameters[8];
    advCoeffX = parameters[11];
    advCoeffY = parameters[12];
    gamma = parameters[9];
    tau = parameters[10];

    quadMesh mesh(domainX, domainY, nX, nY);
    Eigen::MatrixXi elemDegrees = localDegree * Eigen::MatrixXi::Ones(2, mesh.getNumberElements());
    fespace=FeSpace(maxDegree, mesh, elemDegrees);
    solved =false;
}

double AdvectionDiffusionProblem::computeEnergyError(std::function<double(double, double)> exactSolution) {
    Error error (exactSolution, solutionCoeff);
    double energyError = error.computeEnergyError(fespace, gamma);
    std::cout<< "Energy error: " << energyError << std::endl;
    return energyError;
}

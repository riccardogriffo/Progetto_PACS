//
// Created by Alberto Chiappa on 26/10/17.
//

#include "GaussLobattoQuadrature.h"

std::pair<Eigen::VectorXd, Eigen::VectorXd> GaussLobattoQuadrature::computeNodesAndWeights(unsigned int numNodes,
                                                                                           std::pair<double, double> domain) {
/* *****************************************************************************
!    LOB_R8 generates a Gauss-Lobatto quadrature rule.
!   The modified Golub-Welsch alghoritm is used. See Golub 73.
!   Discussion:
!   Given numNodes (total number of quadrature nodes, internal plus endpoints),
!   this routine generates the numNodes-point Gauss-Lobatto-Legendre quadrature formula.
!
!      Integral over(-1,1) of f(X) dX
!      = w(0) f(-1) + sum from k=1 to k=n of w(k) f(x(k))
!      + w(n+1) f(1) + R(n;f).
!
!   The nodes are returned as csi_GLL(k) = x(k), the weights as wgt_GLL(k)
!   = w(k), k=0,...,numNodes-1.
!   The first numNodes-1 recursion coefficients a(k), b(k), k = 0,...,numNodes-2, for the measure
!   dlambda(x) = dx are are computed atrcording to the Golub-Welsch alghoritm by
!   calling coeflege( numNodes-1, a(1:numNodes-1), b(1:numNodes-1) ).
!   Then the modified numNodes-1-th recursion coefficients a(numNodes-1) and b(numNodes-1) are computed
!   acccording to the alghoritm descibed in Golub 1973.
!   Finally, the modified tridiagonal Jacoby matrix is computed and
!   the numNodes nodes and numNodes weights are computed in terms of the
!   eigenvalues and first component of the normalized eigenvectors of
!   the slightly modified Jacobi matrix of order  numNodes.
!
!   Author:
!   Modified by G. Tumolo 29 June 2011 on the basis of an existing code developed
!   on 10 August 2007 by Walter Gautschi
!
!   References:
!
!   Gene H. Gulob,
!   Some modified matrix eigenvalue problems
!   SIAM Review, Vol.15, No.2, April 1973 pp 318-334 and in particular pp 333-334.
!
!   Walter Gautschi,
!   Algorithm 726:
!   ORTHPOL - A Package of Routines for Generating Orthogonal
!   Polynomials and Gauss-Type Quadrature Rules,
!   ACM Transactions on Mathematical Software,
!   Volume 20, Number 1, March 1994, pages 21-62.
!
!   Parameters:
!
!   Int numNodes, the total ( numNodes-2 interior +  2 endoints) number of
!   quadrature nodes in the Gauss-Lobatto formula.
!
!   double csi_GLL(numNodes), the nodes in increasing order.
!   csi_GLL(k) = x(k), k=0,...,numNodes-1
!
!   double wgt_GLL(numNodes), the GLL weights
!   wgt_GLL(k) = w(k), k=0,...,numNodes-1
!
!
!   double left = -1.0, right = 1.0, the prescribed left and right
!   endpoints x(0) and x(numNodes-1) of the Gauss-Lobatto formula
!
     */

    Eigen::VectorXd csi_GLL(numNodes);
    Eigen::VectorXd wgt_GLL(numNodes);

    double left = -1., right = 1.;	//For Gauss-Legendre-Lobatto

    Eigen::VectorXd a(numNodes);
    Eigen::VectorXd b(numNodes), temp(numNodes);
    Eigen :: MatrixXd JacM=Eigen::MatrixXd::Zero(numNodes,numNodes);     //USARE MATRICE SPARSA!!!      //initialization to zero
    Eigen::VectorXi index(numNodes);
    double det, p0l, p0r, p1l, p1r, pm1l, pm1r;

    LegendreCoefficients(numNodes-1, a, b);

    p0l = 0.;
    p0r = 0.;
    p1l = 1.;
    p1r = 1.;

    for(int k = 0; k < numNodes-1; k++){
        pm1l = p0l;
        p0l = p1l;
        pm1r = p0r;
        p0r = p1r;
        p1l = (left - a(k)) * p0l - b(k) * pm1l;
        p1r = (right - a(k)) * p0r - b(k) * pm1r;
    }
    det = p1l * p0r - p1r * p0l;
    a(numNodes-1) = (left * p1l * p0r - right * p1r * p0l) / det;
    b(numNodes-1) = (right - left) * p1l * p1r / det;


    JacM(0,0) = a(0);
    for(int i = 1; i < numNodes; i++){
        JacM(i, i) = a(i);
        JacM(i-1, i) = sqrt(b(i));
        JacM(i, i-1) = JacM(i-1, i);
    }

    Eigen::ComplexEigenSolver<Eigen::MatrixXd> ces(JacM);    //calcolo autovettori e autovalori
    csi_GLL = ces.eigenvalues().real();
    wgt_GLL = ces.eigenvectors().row(0).real();
    wgt_GLL.normalized();    // normalized

    for(int i = 0; i < numNodes; i++){
        wgt_GLL(i) = 2. * pow(wgt_GLL(i), 2); // 2 is the interval size
    }

    unsigned int n{0};

    std::generate(index.data(), index.data()+numNodes, [&]{ return n++; });      //inizializzo il vettore degli indici

    std::sort( index.data(), index.data()+numNodes, [&](int i1, int i2) -> bool { return csi_GLL(i1) < csi_GLL(i2); } );    //trova gli indici del vettore ordinato

    std::sort(csi_GLL.data(), csi_GLL.data()+numNodes);         //sort csi_GLL

    temp = wgt_GLL;
    for(int i = 0; i < numNodes; i++){
        wgt_GLL(i) = temp(index(i));  // metto i pesi in corrispondenza dei nodi
    };

    csi_GLL = (domain.second-domain.first)/2*csi_GLL + (domain.second+domain.first)/2*Eigen::VectorXd::Ones(numNodes);
    wgt_GLL = (domain.second-domain.first)/2*wgt_GLL;


    return std::make_pair(csi_GLL, wgt_GLL);
}

void GaussLobattoQuadrature::LegendreCoefficients(unsigned const int M, Eigen::VectorXd &a, Eigen::VectorXd &b) {
//   Computes coefficients a and b of the three terms recurrence formula for
//   Legendre polynomials and stores them in a, b vectors of double of size M.
//   Remark: b(M-1) is computed but NOT used, since the Jacobi matrix, associated
//          to Legendre polynomials, is defined as the MxM tridiagonal matrix
//           with a(0:M-1) as the main diagonal and with b(0:M-2) as the first
//           subdiagonal and supradiagonal.

    b(0) = 2.;
    a.setZero();
    for(int k = 1; k < M; k++){
        b(k) = 1. / (4.- 1. / (k * k));         //scrivere elevamento al quadrato
    }
}

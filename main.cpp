#include <iostream>
#include <cstdlib>
#include "SimplexAbstract.hpp"
#include <math.h>
#include <list>
#include <time.h>
#include "SimplexAlpha.hpp"
#include "Covering_CORE.hpp"
//#include "OptimizedCovering.hpp"
//#include "Lex.hpp"
#include <mpi.h>
#include <stdio.h>
#include <cstdio>
#include "LocalSearch.hpp"

using namespace std;

Complex komplex;
SimplicialMap map0, map1;
ComplexProduct KxK;
SubComplexJ L;

Covering covering;

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

int main(int argc, char **argv) {

        srand(time(NULL));

        int M = 1000;
        double r = 0.1;
        int Ns = 1;

        int maxInt = 2;
        int numOfMaxSimplex = 3;
        komplex.initComplex(numOfMaxSimplex, maxInt + 1);

        komplex.K.A[0][0].initSimplex(2, 0, 1);
        komplex.K.A[0][1].initSimplex(2, 1, 2);
        komplex.K.A[0][2].initSimplex(2, 0, 2);

        KxK.initComplexProduct(komplex);
        KxK.escMaximalSimplices();
        map1.projection1(KxK);
        map0.projection2(KxK);

        L.initSubComplexJ(18);
        int counting = 0;

        for (int i = 0; i < KxK.listOfFacets.m; i++) {
                for (int j = 0; j < KxK.listOfFacets.rowLength.getA(0, i); j++) {
                        L.initA(counting, KxK.listOfFacets.A[i][j]);
                        counting += 1;
                }
        }
	

	 L.initZero_Skeleton();

	 L.escSubComplexJ();
	
         komplex.initAdjMat();
         komplex.graph.addWeight(0, 1, 1);
         komplex.graph.addWeight(0, 2, 1);
         komplex.graph.addWeight(1, 2, 1);

	
	covering.initCovering(L, komplex, map1, map0, M, r, argc, argv);	
	covering.runCovering(L, komplex, map1, map0, M, r, argc, argv);
	covering.endCovering();

        cout << "\n\n\n";

    return 0;
}

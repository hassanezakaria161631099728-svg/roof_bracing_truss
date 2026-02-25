#include "expressions.h"
#include "fem.h"
#include "io.h"
#include <vector>
#include <iostream>
int main()
{
double h = 3;
int n = 3;
double L = 4;
double E = 2.1e11;
double Aver = 0.02;
double Iver = 8e-5;
double Ahor = 0.015;
double Ihor = 6e-5;
double Adiag = 0.01;
double Idiag = 1e-5;
double q = 0;
// 1. generate model
std::vector<Node> nodes = generateNodes(h, n, L);
std::vector<Element> elements = generateElements(n, E, Aver, Iver,Ahor, Ihor, Adiag, Idiag, q);
std::vector<DOF> dofs = generateDOFs(nodes);
std::vector<int> bc = extractBC(dofs);
int totalDOF = dofs.size();
//2. compute geometry
for (auto& e: elements)
    computeGeometry(e, nodes);
//3. initialize global system
std::vector<std::vector<double>>K(totalDOF,std::vector<double>(totalDOF,0.0));
std::vector<double>F(totalDOF,0.0);

//4.assemble global matrices
for (const auto& e: elements)
{
    double kg[6][6];
    double fe[6];
    globalStiffness(e, kg);
    equivalentLoad(e, fe);
    assembleGlobal(K, kg, e.n1, e.n2);
    assembleLoad(F, fe, e.n1, e.n2);

}
//  nodal loads
F[3]+= -1000;
F[6]+= -1000;
F[9]+= -1000;

//5. partition DOFs
std::vector<int> freeDOF,fixedDOF;
    partitionDOF(totalDOF,bc,freeDOF,fixedDOF);
//6. build reduced system
    std::vector<std::vector<double>> Kff;
    std::vector<double> Ff;

    buildReducedSystem(K,F,freeDOF,Kff,Ff);
//7. solve
    std::vector<double> Uf(freeDOF.size());
    solveSystem(Kff,Ff,Uf);
//8. expand full displacement vector
    std::vector<double> U;
    expandDisplacements(totalDOF,freeDOF,Uf,U);
//9. compute reactions
    std::vector<double> R;
    computeReactions(K,F,U,fixedDOF,R);
//10. outputs
plotStructure(nodes, elements, 800, 600);
    std::string filename = "results.txt";
writeresults(filename, nodes, elements, dofs, U, R, fixedDOF);
    return 0;
}


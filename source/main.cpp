#include "expressions.h"
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
std::vector<Node> nodes = generateNodes(h, n, L);
std::vector<Element> elements = generateElements(n, E, Aver, Iver,Ahor, Ihor, Adiag, Idiag, q);
std::vector<DOF> dofs = generateDOFs(nodes);
//plotStructure(nodes, elements, 800, 600);
writeNodesToFile("results.txt",nodes,elements,dofs);
    return 0;
}

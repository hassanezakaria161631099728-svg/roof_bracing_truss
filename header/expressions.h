#pragma once
#include <vector>
#include <string>
struct Node
{int id; double x; double y;};
std::vector<Node> generateNodes(double h, int n, double L);
struct Element {int id; int n1; int n2;
double E;
double A;
double I;
double q;
std::string load_type; // "local" or "globalS"
double L,c,s;
};
std::vector<Element> generateElements(int n, double E, double Aver, double Iver,
                                      double Ahor, double Ihor, double Adiag,double Idiag, double q);
struct DOF
{
    int id;         // global DOF number
    int nodeId;     // node it belongs to
    std::string type;  // "ux", "uy", "rz"
    bool constrained;
};
std::vector<DOF> generateDOFs(const std::vector<Node>& nodes);

void plotStructure(const std::vector<Node>& nodes,
                   const std::vector<Element>& elements,
                   int windowWidth,
                   int windowHeight);

std::vector<int> extractBC(const std::vector<DOF>& dofs);

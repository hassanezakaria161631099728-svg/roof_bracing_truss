#ifndef FEM_H
#define FEM_H
#include "expressions.h"
#include <vector>
#include <string>




// Geometry
void computeGeometry(Element &e, const std::vector<Node> &nodes);

// Element routines
void localStiffness(const Element &e, double k[6][6]);
void transformationMatrix(const Element &e, double T[6][6]);
void globalStiffness(const Element &e, double kg[6][6]);
void equivalentLoad(const Element &e, double fe[6]);

// Assembly
void assembleGlobal(std::vector<std::vector<double>>& K,double kg[6][6],int n1,int n2);
void assembleLoad(std::vector<double>& F,double fe[6],int n1,int n2);

// DOF handling
void partitionDOF(int totalDOF,const std::vector<int>& constrained,
                  std::vector<int>& freeDOF,std::vector<int>& fixedDOF);

void buildReducedSystem(const std::vector<std::vector<double>>& K,
                        const std::vector<double>& F,
                        const std::vector<int>& freeDOF,
                        std::vector<std::vector<double>>& Kff,
                        std::vector<double>& Ff);

void expandDisplacements(int totalDOF,
                         const std::vector<int>& freeDOF,
                         const std::vector<double>& Uf,
                         std::vector<double>& U);

void computeReactions(const std::vector<std::vector<double>>& K,
                      const std::vector<double>& F,
                      const std::vector<double>& U,
                      const std::vector<int>& fixedDOF,
                      std::vector<double>& R);

// Linear solver
void solveSystem(const std::vector<std::vector<double>>& K,
                 const std::vector<double>& F,
                 std::vector<double>& U);

void elementInternalForces(
    const Element& e,
    const std::vector<double>& U,
    double fl[6]);

#endif


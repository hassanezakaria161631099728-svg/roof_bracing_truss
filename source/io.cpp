#include "io.h"
#include <fstream>
#include <vector>
#include <string>
#include "expressions.h"
#include <iomanip>
#include "fem.h"
void writeresults(const std::string& filename,
                      const std::vector<Node>& nodes,
                      const std::vector<Element>& elements,
                      const std::vector<DOF>& dofs,
                      const std::vector<double>& U,
                      const std::vector<double>& R,
                    const std::vector<int>& fixedDOF
                      )

{
    std::ofstream file(filename);

    if(!file)
        return;

    file << "=========================\n";
    file << "        NODE LIST        \n";
    file << "=========================\n\n";

    file << std::setw(10) << "ID"
         << std::setw(15) << "X"
         << std::setw(15) << "Y"
         << "\n";

    file << "--------------------------------------\n";

    for(const auto& node : nodes)
    {
        file << std::setw(10) << node.id
             << std::setw(15) << std::fixed << std::setprecision(2) << node.x
             << std::setw(15) << std::fixed << std::setprecision(2) << node.y
             << "\n";
    }
    file << "=========================\n";
    file << "        ELEMENT CONNECTIVITY        \n";
    file << "=========================\n\n";

    file << std::setw(10) << "ID"
         << std::setw(15) << "Node 1"
         << std::setw(15) << "Node 2"
         << std::setw(15) << "E"
         << std::setw(15) << "A"
         << std::setw(15) << "I"
         << std::setw(15) << "q"
         << std::setw(15) << "load type"
         << "\n";

    file << "--------------------------------------\n";

    for(const auto& elem : elements)
    {
        file << std::setw(10) << elem.id
             << std::setw(15) << std::fixed << std::setprecision(2) << elem.n1
             << std::setw(15) << std::fixed << std::setprecision(2) << elem.n2
             << std::setw(15) << std::scientific << std::setprecision(2) << elem.E
             << std::setw(15) << std::fixed << std::setprecision(2) << elem.A
             << std::setw(15) << std::scientific << std::setprecision(2) << elem.I
             << std::setw(15) << std::fixed << std::setprecision(2) << elem.q
             << std::setw(15) << elem.load_type

             << "\n";
    }

    file << "=========================\n";
    file << "        DOF LIST        \n";
    file << "=========================\n\n";

    file << std::setw(10) << "ID"
         << std::setw(15) << "node ID"
         << std::setw(15) << "freedom type"
         << std::setw(15) << "freedom situation"
         << "\n";

    file << "--------------------------------------\n";

    for(const auto& dof : dofs)
    {
        file << std::setw(10) << dof.id
             << std::setw(15) << dof.nodeId
             << std::setw(15) << dof.type
             << std::setw(15) << (dof.constrained? "Constrained" : "Free")
             << "\n";
    }

    file << "=========================================\n";
    file << "        NODAL DISPLACEMENTS [cm]\n";
    file << "=========================================\n\n";

    file << std::setw(10) << "Node"
        << std::setw(18) << "ux"
        << std::setw(18) << "uy"
        << std::setw(18) << "rotz\n";

    file << "--------------------------------------------------------------\n";
    int nNodes= nodes.size();

    for(int i=0;i<nNodes;i++)
    {
        file << std::setw(10) << i+1
            << std::setw(18) << U[3*i]
            << std::setw(18) << U[3*i+1]
            << std::setw(18) << U[3*i+2]
            << "\n";
    }

    file << "\n";

    file << "=========================================\n";
    file << "        SUPPORT REACTIONS [kg & m]\n";
    file << "=========================================\n\n";

    file << std::setw(10) << "Node"
        << std::setw(18) << "Rx"
        << std::setw(18) << "Ry"
        << std::setw(18) << "Mz\n";

    file << "--------------------------------------------------------------\n";

    std::vector<bool> isFixed(3*nodes.size(),false);
    for(int d : fixedDOF)
        isFixed[d] = true;

    for(int i=0;i<nNodes;i++)
    {
        if(isFixed[3*i] || isFixed[3*i+1] || isFixed[3*i+2])
        {
            file << std::setw(10) << i+1
                << std::setw(18) << R[3*i]
                << std::setw(18) << R[3*i+1]
                << std::setw(18) << R[3*i+2]
                << "\n";
        }
    }

    file << "\n";

 file << "ELEMENT FORCES (LOCAL AXIS) [kg & m]\n";
    file << "--------------------------------------------------------------\n";
    file << std::setw(10) << "Element"
        << std::setw(10) << "Node"
        << std::setw(15) << "N"
        << std::setw(15) << "V"
        << std::setw(15) << "M\n";
    file << "--------------------------------------------------------------\n";

    for(size_t i=0;i<elements.size();i++)
    {
        double fl[6];
        elementInternalForces(elements[i],U,fl);

        // Node 1 end
        file << std::setw(10) << i+1
            << std::setw(10) << elements[i].n1
            << std::setw(15) << fl[0]
            << std::setw(15) << fl[1]
            << std::setw(15) << fl[2]
            << "\n";

        // Node 2 end
        file << std::setw(10) << i+1
            << std::setw(10) << elements[i].n2
            << std::setw(15) << fl[3]
            << std::setw(15) << fl[4]
            << std::setw(15) << fl[5]
            << "\n";
    }

    file.close();
}



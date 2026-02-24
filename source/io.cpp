#include "io.h"
#include <fstream>
#include <iomanip>

void writeNodesToFile(const std::string& filename,
                      const std::vector<Node>& nodes,
                      const std::vector<Element>& elements,
                      const std::vector<DOF> dofs)

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
         << "\n";

    file << "--------------------------------------\n";

    for(const auto& elem : elements)
    {
        file << std::setw(10) << elem.id
             << std::setw(15) << std::fixed << std::setprecision(2) << elem.node1
             << std::setw(15) << std::fixed << std::setprecision(2) << elem.node2
             << std::setw(15) << std::scientific << std::setprecision(2) << elem.E
             << std::setw(15) << std::fixed << std::setprecision(2) << elem.A
             << std::setw(15) << std::scientific << std::setprecision(2) << elem.I
             << std::setw(15) << std::fixed << std::setprecision(2) << elem.q

             << "\n";
    }

    file << "=========================\n";
    file << "        DOF LIST        \n";
    file << "=========================\n\n";

    file << std::setw(10) << "ID"
         << std::setw(15) << "node ID"
         << std::setw(15) << "Y"
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

    file.close();
}

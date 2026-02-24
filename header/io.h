#pragma once
#include <vector>
#include <string>
#include "expressions.h"

void writeNodesToFile(const std::string& filename,
                      const std::vector<Node>& nodes,
                      const std::vector<Element>& elements,
                      const std::vector<DOF> dofs);

#include "expressions.h"
#include <vector>
#include <string>
#include <cmath>
#include <graphics.h>

std::vector<Node> generateNodes(double h, int n, double L)
{
    std::vector<Node> nodes;

    for(int col = 0; col < 2; ++col)
    {
        double xCoord = (col == 0) ? 0.0 : L;

        for(int i = 0; i <= n; ++i)
        {
            Node node;
            node.id = col*(n+1) + i + 1;
            node.x  = xCoord;
            node.y  = i * h;

            nodes.push_back(node);
        }
    }

    return nodes;
}


std::vector<Element> generateElements(int n, double E, double Aver, double Iver,
                                      double Ahor, double Ihor, double Adiag, double Idiag, double q)
{
    std::vector<Element> elements;
    int id = 1;

    for(int i = 1; i <= n; ++i)
    {
        int leftBottom  = i;
        int leftTop     = i + 1;
        int rightBottom = n + 1 + i;
        int rightTop    = n + 1 + i + 1;

        // ---- Left column ----
        elements.push_back({
            id++, leftBottom, leftTop,
            E, Aver, Iver,
            q, "local"
        });

        // ---- Right column ----
        elements.push_back({
            id++, rightBottom, rightTop,
            E, Aver, Iver,
            q, "local"
        });

        // ---- Beam ----
        elements.push_back({
            id++, leftTop, rightTop,
            E, Ahor, Ihor,
            q, "local"
        });

        // ---- Diagonal 1 ----
        elements.push_back({
            id++, leftBottom, rightTop,
            E, Adiag, Idiag,
            q, "local"
        });

        // ---- Diagonal 2 ----
        elements.push_back({
            id++, leftTop, rightBottom,
            E, Adiag, Idiag,
            q, "local"
        });
    }

    return elements;
}

std::vector<DOF> generateDOFs(const std::vector<Node>& nodes)
{
    std::vector<DOF> dofs;
    int id = 1;

    for(const auto& node : nodes)
    {
        bool isGround = (node.y == 0.0);

        // UX
        dofs.push_back({
            id++,
            node.id,
            "ux",
            isGround   // constrained if ground
        });

        // UY
        dofs.push_back({
            id++,
            node.id,
            "uy",
            isGround   // constrained if ground
        });

        // RZ
        dofs.push_back({
            id++,
            node.id,
            "rz",
            false      // rotation always free
        });
    }

    return dofs;
}

void plotStructure(const std::vector<Node>& nodes,
                   const std::vector<Element>& elements,
                   int windowWidth,
                   int windowHeight)
{
    initwindow(windowWidth, windowHeight, "2D Frame");

    int scale = 60;
    int margin = 100;

    // -------------------------
    // DRAW ELEMENTS
    // -------------------------
    for(const auto& elem : elements)
    {
        const Node& n1 = nodes[elem.node1 - 1];
        const Node& n2 = nodes[elem.node2 - 1];

        int x1 = margin + n1.x * scale;
        int y1 = windowHeight - (margin + n1.y * scale);

        int x2 = margin + n2.x * scale;
        int y2 = windowHeight - (margin + n2.y * scale);

        line(x1, y1, x2, y2);

        // Element ID at midpoint
        int midX = (x1 + x2) / 2;
        int midY = (y1 + y2) / 2;

        std::string label = std::to_string(elem.id);
        outtextxy(midX, midY, const_cast<char*>(label.c_str()));
    }

    // -------------------------
    // DRAW NODES
    // -------------------------
    for(const auto& node : nodes)
    {
        int x = margin + node.x * scale;
        int y = windowHeight - (margin + node.y * scale);

        circle(x, y, 5);   // draw node

        std::string label;

        if(std::abs(node.y) < 1e-9)
            label = "p" + std::to_string(node.id);
        else
            label = "n" + std::to_string(node.id);

        outtextxy(x + 8, y - 8, const_cast<char*>(label.c_str()));
    }

    getch();
    closegraph();
}

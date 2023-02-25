#ifndef _IO_H_
#define _IO_H_

#include "solver.h"

class IO
{
private:
    Solver *solp;

public:
    IO(Solver &solver) : solp(&solver) {}
    void Print(NodeE &node, std::ostream &out = std::cout);
    void Print(EdgeE &edge, std::ostream &out = std::cout);
    void Print(ElemE &elem, std::ostream &out = std::cout);
    void Print(NodeU &node, std::ostream &out = std::cout);
    void Print(EdgeU &edge, std::ostream &out = std::cout);
    void Print(ElemU &elem, std::ostream &out = std::cout);
    void Print(SubElem &sub, std::ostream &out = std::cout);
    void Print(SubEdge &sub, std::ostream &out = std::cout);
    void PrintNodeE(const char *file);
    void PrintEdgeE(const char *file);
    void PrintElemE(const char *file);
    void PrintNodeU(const char *file);
    void PrintEdgeU(const char *file);
    void PrintElemU(const char *file);
    void PrintRefBuf(const char *file);
    void PrintQuadRule(const char *file);
};

#endif
#include "io.h"

void IO::Print(NodeE &node, std::ostream &out)
{
    out << ' ' << std::setw(12) << node.pos.x
        << ' ' << std::setw(12) << node.pos.y;
    out << '\n';
}

void IO::Print(NodeU &node, std::ostream &out)
{
    out << ' ' << std::setw(6) << node.par
        << ' ' << std::setw(12) << node.pos.x
        << ' ' << std::setw(12) << node.pos.y;
    out << '\n';
}

void IO::Print(EdgeE &edge, std::ostream &out)
{
    out << ' ' << std::setw(6) << edge.beg
        << ' ' << std::setw(6) << edge.end;
    out << '\n';
}

void IO::Print(EdgeU &edge, std::ostream &out)
{
    out << ' ' << std::setw(6) << edge.beg
        << ' ' << std::setw(6) << edge.end
        << '\n';
}

void IO::Print(ElemE &elem, std::ostream &out)
{
    int i;

    out << "vertices:\n";
    for (i = 0; i < 4; i++)
    {
        out << ' ' << std::setw(6) << elem.vert_buf[i];
    }
    out << '\n';
    out << "edges:\n";
    for (i = 0; i < 4; i++)
    {
        out << ' ' << std::setw(6) << elem.edge_buf[i];
    }
    out << '\n';
    out << "DoFs:\n";
    for (i = 0; i < elem.ndof; i++)
    {
        out << std::setw(12) << elem.udof[i] << '\n';
    }
    out << '\n';
}

void IO::Print(ElemU &elem, std::ostream &out)
{
    int i;
    out << "vertices:\n";
    for (i = 0; i < 4; i++)
    {
        out << ' ' << std::setw(6) << elem.vert_buf[i];
    }
    out << '\n';
    out << "edges:\n";
    for (i = 0; i < 4; i++)
    {
        out << ' ' << std::setw(6) << elem.edge_buf[i];
    }
    out << '\n';
}

void IO::Print(SubElem &sub, std::ostream &out)
{
    out << ' ' << std::setw(6) << sub.par
        << ' ' << std::setw(6) << sub.nseg
        << ' ' << std::setw(6) << sub.seg[0]
        << ' ' << std::setw(6) << sub.seg[1]
        << ' ' << std::setw(6) << sub.seg[2]
        << ' ' << std::setw(6) << sub.seg[3]
        << ' ' << std::setw(8) << sub.x0
        << ' ' << std::setw(8) << sub.x1
        << ' ' << std::setw(8) << sub.y0
        << ' ' << std::setw(8) << sub.y1
        << '\n';
}

void IO::Print(POI &poi, std::ostream &out)
{
    out << ' ' << std::setw(6) << poi.par
        << ' ' << std::setw(6) << poi.type
        << ' ' << std::setw(8) << std::setprecision(8) << poi.pos.x
        << ' ' << std::setw(8) << std::setprecision(8) << poi.pos.y
        << '\n';
}

void IO::Print(Segment &seg, std::ostream &out)
{
    out << ' ' << std::setw(6) << seg.beg
        << ' ' << std::setw(6) << seg.end
        << '\n';
}

void IO::PrintNodeE(const char *file)
{
    std::ofstream out;
    out.open(file);
    out << std::fixed << std::setprecision(8);
    for (int i = 0; i < solp->N_node; i++)
    {
        out << ' ' << std::setw(6) << i;
        Print(solp->node_e[i], out);
    }
    out.close();
}

void IO::PrintEdgeE(const char *file)
{
    std::ofstream out;
    out.open(file);
    for (int i = 0; i < solp->N_edge; i++)
    {
        out << ' ' << std::setw(6) << i;
        Print(solp->edge_e[i], out);
    }
    out.close();
}

void IO::PrintElemE(const char *file)
{
    std::ofstream out;
    out.open(file);
    for (int i = 0; i < solp->N_elem; i++)
    {
        out << "[ElemE " << i << "]\n";
        Print(solp->elem_e[i], out);
    }
    out.close();
}

void IO::PrintNodeU(const char *file)
{
    std::ofstream out;
    out.open(file);
    out << std::fixed << std::setprecision(8);
    for (int i = 0; i < solp->N_node; i++)
    {
        out << ' ' << std::setw(6) << i;
        Print(solp->node_u[i], out);
    }
    out.close();
}

void IO::PrintEdgeU(const char *file)
{
    std::ofstream out;
    out.open(file);
    for (int i = 0; i < solp->N_edge; i++)
    {
        out << ' ' << std::setw(6) << i;
        Print(solp->edge_u[i], out);
    }
    out.close();
}

void IO::PrintElemU(const char *file)
{
    std::ofstream out;
    out.open(file);
    for (int i = 0; i < solp->N_elem; i++)
    {
        out << "[ElemU " << i << "]\n";
        Print(solp->elem_u[i], out);
    }
    out.close();
}

void IO::PrintPOIBufU(const char *file)
{
    std::ofstream out;
    out.open(file);
    for (int i = 0; i < solp->N_elem; i++)
    {
        out << "[ElemU " << i << "]\n";
        out << "POI_BUF: " << '\n';
        for (int j = 0; j < solp->elem_u[i].npoi; j++)
        {
            out << ' ' << std::setw(6) << j;
            Print(solp->elem_u[i].poi_buf[j], out);
        }
    }
    out.close();
}

void IO::PrintSegmentBufU(const char *file)
{
    std::ofstream out;
    out.open(file);
    for (int i = 0; i < solp->N_elem; i++)
    {
        out << "[ElemU " << i << "]\n";
        out << "SEG_BUF: " << '\n';
        for (int j = 0; j < solp->elem_u[i].nseg; j++)
        {
            out << ' ' << std::setw(6) << j;
            Print(solp->elem_u[i].seg_buf[j], out);
        }
    }
    out.close();
}

void IO::PrintSubElemBufU(const char *file)
{
    std::ofstream out;
    out.open(file);
    for (int i = 0; i < solp->N_elem; i++)
    {
        out << "[ElemU " << i << "]\n";
        out << "SUB_BUF: " << '\n';
        for (int j = 0; j < solp->elem_u[i].nsub; j++)
        {
            out << ' ' << std::setw(6) << j;
            Print(solp->elem_u[i].sub_buf[j], out);
        }
    }
    out.close();
}

void IO::PrintPOIBufE(const char *file)
{
    std::ofstream out;
    out.open(file);
    for (int i = 0; i < solp->N_elem; i++)
    {
        out << "[ElemE " << i << "]\n";
        out << "POI_BUF: " << '\n';
        for (int j = 0; j < solp->elem_e[i].npoi; j++)
        {
            out << ' ' << std::setw(6) << j;
            Print(solp->elem_e[i].poi_buf[j], out);
        }
    }
    out.close();
}

void IO::PrintSegmentBufE(const char *file)
{
    std::ofstream out;
    out.open(file);
    for (int i = 0; i < solp->N_elem; i++)
    {
        out << "[ElemE " << i << "]\n";
        out << "SEG_BUF: " << '\n';
        for (int j = 0; j < solp->elem_e[i].nseg; j++)
        {
            out << ' ' << std::setw(6) << j;
            Print(solp->elem_e[i].seg_buf[j], out);
        }
    }
    out.close();
}

void IO::PrintSubElemBufE(const char *file)
{
    std::ofstream out;
    out.open(file);
    for (int i = 0; i < solp->N_elem; i++)
    {
        out << "[ElemE " << i << "]\n";
        out << "SUB_BUF: " << '\n';
        for (int j = 0; j < solp->elem_e[i].nsub; j++)
        {
            out << ' ' << std::setw(6) << j;
            Print(solp->elem_e[i].sub_buf[j], out);
        }
    }
    out.close();
}

void IO::PrintQuadRule(const char *file)
{
    std::ofstream out;
    out.open(file);
    out << std::fixed << std::setprecision(8);
    for (int i = 0; i < solp->N_elem; i++)
    {
        out << "[ElemE " << i << "]\n";
        out << "BBOX: " << '\n';
        out << ' ' << std::setw(12) << solp->elem_e[i].x0
            << ' ' << std::setw(12) << solp->elem_e[i].x1
            << ' ' << std::setw(12) << solp->elem_e[i].y0
            << ' ' << std::setw(12) << solp->elem_e[i].y1
            << '\n';
        out << "QUAD_E: " << '\n';
        for (int j = 0; j < solp->elem_e[i].nsub; j++)
        {
            for (int k = 0; k < solp->N_GL_2D; k++)
            {
                out << ' ' << std::setw(12) << solp->elem_e[i].sub_buf[j].qp[k].x
                    << ' ' << std::setw(12) << solp->elem_e[i].sub_buf[j].qp[k].y
                    << ' ' << std::setw(12) << solp->elem_e[i].sub_buf[j].qp[k].w
                    << '\n';
            }
        }
        out << "QUAD_U: " << '\n';
        for (int j = 0; j < solp->elem_u[i].nsub; j++)
        {
            for (int k = 0; k < solp->N_GL_2D; k++)
            {
                out << ' ' << std::setw(12) << solp->elem_u[i].sub_buf[j].qp[k].x
                    << ' ' << std::setw(12) << solp->elem_u[i].sub_buf[j].qp[k].y
                    << ' ' << std::setw(12) << solp->elem_u[i].sub_buf[j].qp[k].w
                    << '\n';
            }
        }
    }
    out.close();
}
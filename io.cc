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
    out << ' ' << std::setw(8) << sub.p0.x
        << ' ' << std::setw(8) << sub.p0.y
        << ' ' << std::setw(8) << sub.p1.x
        << ' ' << std::setw(8) << sub.p1.y
        << '\n';
}

void IO::Print(SubEdge &sub, std::ostream &out)
{
    out << ' ' << std::setw(8) << sub.p0.x
        << ' ' << std::setw(8) << sub.p0.y
        << ' ' << std::setw(8) << sub.p1.x
        << ' ' << std::setw(8) << sub.p1.y
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
        Print(solp->node_e_buf[i], out);
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
        Print(solp->edge_e_buf[i], out);
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
        Print(solp->elem_e_buf[i], out);
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
        Print(solp->node_u_buf[i], out);
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
        Print(solp->edge_u_buf[i], out);
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
        Print(solp->elem_u_buf[i], out);
    }
    out.close();
}

void IO::PrintRefBuf(const char *file)
{
    std::ofstream out;
    out.open(file);
    out << "SUB_ELEM: " << '\n';
    for (int i = 0; i < solp->ref_elem_buf.nsub; i++)
    {
        out << ' ' << std::setw(6) << i;
        Print(solp->ref_elem_buf.sub_buf[i], out);
    }
    out << "SUB_EDGE_B: " << '\n';
    for (int i = 0; i < solp->ref_edge_b_buf.nsub; i++)
    {
        out << ' ' << std::setw(6) << i;
        Print(solp->ref_edge_b_buf.sub_buf[i], out);
    }
    out << "SUB_EDGE_R: " << '\n';
    for (int i = 0; i < solp->ref_edge_r_buf.nsub; i++)
    {
        out << ' ' << std::setw(6) << i;
        Print(solp->ref_edge_r_buf.sub_buf[i], out);
    }
    out << "SUB_EDGE_T: " << '\n';
    for (int i = 0; i < solp->ref_edge_t_buf.nsub; i++)
    {
        out << ' ' << std::setw(6) << i;
        Print(solp->ref_edge_t_buf.sub_buf[i], out);
    }
    out << "SUB_EDGE_L: " << '\n';
    for (int i = 0; i < solp->ref_edge_l_buf.nsub; i++)
    {
        out << ' ' << std::setw(6) << i;
        Print(solp->ref_edge_l_buf.sub_buf[i], out);
    }

#ifdef FAST_SLDG
    out << std::fixed << std::setprecision(8);
    out << "QP_ELEM_E_BUF: " << '\n';
    for (int i = 0; i < solp->ref_elem_buf.nsub; i++)
    {
        for (int k = 0; k < solp->N_GL_2D; k++)
        {
            out << ' ' << std::setw(13) << solp->ref_elem_buf.qp_e_buf[i][k].x
                << ' ' << std::setw(13) << solp->ref_elem_buf.qp_e_buf[i][k].y
                << ' ' << std::setw(13) << solp->ref_elem_buf.qp_e_buf[i][k].w
                << '\n';
        }
    }

    out << "QP_ELEM_U_BUF: " << '\n';
    for (int i = 0; i < solp->ref_elem_buf.nsub; i++)
    {
        for (int k = 0; k < solp->N_GL_2D; k++)
        {
            out << ' ' << std::setw(13) << solp->ref_elem_buf.qp_u_buf[i][k].x
                << ' ' << std::setw(13) << solp->ref_elem_buf.qp_u_buf[i][k].y
                << ' ' << std::setw(13) << solp->ref_elem_buf.qp_u_buf[i][k].w
                << '\n';
        }
    }

    out << "V_QP_ELEM_E_BUF: " << '\n';
    for (int i = 0; i < solp->ref_elem_buf.nsub; i++)
    {
        for (int j = 0; j < solp->N_dof_loc; j++)
        {
            for (int k = 0; k < solp->N_GL_2D; k++)
            {
                out << ' ' << std::setw(6) << i
                    << ' ' << std::setw(6) << j
                    << ' ' << std::setw(6) << k
                    << ' ' << std::setw(13) << solp->ref_elem_buf.v_qp_e_buf[i][j][k] << '\n';
            }
        }
    }

    out << "V_QP_ELEM_U_BUF: " << '\n';
    for (int i = 0; i < solp->ref_elem_buf.nsub; i++)
    {
        for (int j = 0; j < solp->N_dof_loc; j++)
        {
            for (int k = 0; k < solp->N_GL_2D; k++)
            {
                out << ' ' << std::setw(6) << i
                    << ' ' << std::setw(6) << j
                    << ' ' << std::setw(6) << k
                    << ' ' << std::setw(13) << solp->ref_elem_buf.v_qp_u_buf[i][j][k] << '\n';
            }
        }
    }

    out << "QP_EDGE_B_E_BUF: " << '\n';
    for (int i = 0; i < solp->ref_edge_b_buf.nsub; i++)
    {
        for (int k = 0; k < solp->N_GL_1D; k++)
        {
            out << ' ' << std::setw(13) << solp->ref_edge_b_buf.qp_e_buf[i][k].x
                << ' ' << std::setw(13) << solp->ref_edge_b_buf.qp_e_buf[i][k].y
                << ' ' << std::setw(13) << solp->ref_edge_b_buf.qp_e_buf[i][k].w
                << '\n';
        }
    }

    out << "QP_EDGE_B_U_BUF: " << '\n';
    for (int i = 0; i < solp->ref_edge_b_buf.nsub; i++)
    {
        for (int k = 0; k < solp->N_GL_1D; k++)
        {
            out << ' ' << std::setw(13) << solp->ref_edge_b_buf.qp_u_buf[i][k].x
                << ' ' << std::setw(13) << solp->ref_edge_b_buf.qp_u_buf[i][k].y
                << ' ' << std::setw(13) << solp->ref_edge_b_buf.qp_u_buf[i][k].w
                << '\n';
        }
    }

    out << "V_QP_EDGE_B_E_BUF: " << '\n';
    for (int i = 0; i < solp->ref_edge_b_buf.nsub; i++)
    {
        for (int j = 0; j < solp->N_dof_loc; j++)
        {
            for (int k = 0; k < solp->N_GL_1D; k++)
            {
                out << ' ' << std::setw(6) << i
                    << ' ' << std::setw(6) << j
                    << ' ' << std::setw(6) << k
                    << ' ' << std::setw(13) << solp->ref_edge_b_buf.v_qp_e_buf[i][j][k][0]
                    << ' ' << std::setw(13) << solp->ref_edge_b_buf.v_qp_e_buf[i][j][k][1] << '\n';
            }
        }
    }

    out << "V_QP_EDGE_B_U_BUF: " << '\n';
    for (int i = 0; i < solp->ref_edge_b_buf.nsub; i++)
    {
        for (int j = 0; j < solp->N_dof_loc; j++)
        {
            for (int k = 0; k < solp->N_GL_1D; k++)
            {
                out << ' ' << std::setw(6) << i
                    << ' ' << std::setw(6) << j
                    << ' ' << std::setw(6) << k
                    << ' ' << std::setw(13) << solp->ref_edge_b_buf.v_qp_u_buf[i][j][k][0]
                    << ' ' << std::setw(13) << solp->ref_edge_b_buf.v_qp_u_buf[i][j][k][1]
                    << '\n';
            }
        }
    }

    out << "QP_EDGE_R_E_BUF: " << '\n';
    for (int i = 0; i < solp->ref_edge_r_buf.nsub; i++)
    {
        for (int k = 0; k < solp->N_GL_1D; k++)
        {
            out << ' ' << std::setw(13) << solp->ref_edge_r_buf.qp_e_buf[i][k].x
                << ' ' << std::setw(13) << solp->ref_edge_r_buf.qp_e_buf[i][k].y
                << ' ' << std::setw(13) << solp->ref_edge_r_buf.qp_e_buf[i][k].w
                << '\n';
        }
    }

    out << "QP_EDGE_R_U_BUF: " << '\n';
    for (int i = 0; i < solp->ref_edge_r_buf.nsub; i++)
    {
        for (int k = 0; k < solp->N_GL_1D; k++)
        {
            out << ' ' << std::setw(13) << solp->ref_edge_r_buf.qp_u_buf[i][k].x
                << ' ' << std::setw(13) << solp->ref_edge_r_buf.qp_u_buf[i][k].y
                << ' ' << std::setw(13) << solp->ref_edge_r_buf.qp_u_buf[i][k].w
                << '\n';
        }
    }

    out << "V_QP_EDGE_R_E_BUF: " << '\n';
    for (int i = 0; i < solp->ref_edge_r_buf.nsub; i++)
    {
        for (int j = 0; j < solp->N_dof_loc; j++)
        {
            for (int k = 0; k < solp->N_GL_1D; k++)
            {
                out << ' ' << std::setw(6) << i
                    << ' ' << std::setw(6) << j
                    << ' ' << std::setw(6) << k
                    << ' ' << std::setw(13) << solp->ref_edge_r_buf.v_qp_e_buf[i][j][k][0]
                    << ' ' << std::setw(13) << solp->ref_edge_r_buf.v_qp_e_buf[i][j][k][1]
                    << '\n';
            }
        }
    }

    out << "V_QP_EDGE_R_U_BUF: " << '\n';
    for (int i = 0; i < solp->ref_edge_r_buf.nsub; i++)
    {
        for (int j = 0; j < solp->N_dof_loc; j++)
        {
            for (int k = 0; k < solp->N_GL_1D; k++)
            {
                out << ' ' << std::setw(6) << i
                    << ' ' << std::setw(6) << j
                    << ' ' << std::setw(6) << k
                    << ' ' << std::setw(13) << solp->ref_edge_r_buf.v_qp_u_buf[i][j][k][0]
                    << ' ' << std::setw(13) << solp->ref_edge_r_buf.v_qp_u_buf[i][j][k][1]
                    << '\n';
            }
        }
    }

    out << "QP_EDGE_T_E_BUF: " << '\n';
    for (int i = 0; i < solp->ref_edge_t_buf.nsub; i++)
    {
        for (int k = 0; k < solp->N_GL_1D; k++)
        {
            out << ' ' << std::setw(13) << solp->ref_edge_t_buf.qp_e_buf[i][k].x
                << ' ' << std::setw(13) << solp->ref_edge_t_buf.qp_e_buf[i][k].y
                << ' ' << std::setw(13) << solp->ref_edge_t_buf.qp_e_buf[i][k].w
                << '\n';
        }
    }

    out << "QP_EDGE_T_U_BUF: " << '\n';
    for (int i = 0; i < solp->ref_edge_t_buf.nsub; i++)
    {
        for (int k = 0; k < solp->N_GL_1D; k++)
        {
            out << ' ' << std::setw(13) << solp->ref_edge_t_buf.qp_u_buf[i][k].x
                << ' ' << std::setw(13) << solp->ref_edge_t_buf.qp_u_buf[i][k].y
                << ' ' << std::setw(13) << solp->ref_edge_t_buf.qp_u_buf[i][k].w
                << '\n';
        }
    }

    out << "V_QP_EDGE_T_E_BUF: " << '\n';
    for (int i = 0; i < solp->ref_edge_t_buf.nsub; i++)
    {
        for (int j = 0; j < solp->N_dof_loc; j++)
        {
            for (int k = 0; k < solp->N_GL_1D; k++)
            {
                out << ' ' << std::setw(6) << i
                    << ' ' << std::setw(6) << j
                    << ' ' << std::setw(6) << k
                    << ' ' << std::setw(13) << solp->ref_edge_t_buf.v_qp_e_buf[i][j][k][0]
                    << ' ' << std::setw(13) << solp->ref_edge_t_buf.v_qp_e_buf[i][j][k][1]
                    << '\n';
            }
        }
    }

    out << "V_QP_EDGE_T_U_BUF: " << '\n';
    for (int i = 0; i < solp->ref_edge_t_buf.nsub; i++)
    {
        for (int j = 0; j < solp->N_dof_loc; j++)
        {
            for (int k = 0; k < solp->N_GL_1D; k++)
            {
                out << ' ' << std::setw(6) << i
                    << ' ' << std::setw(6) << j
                    << ' ' << std::setw(6) << k
                    << ' ' << std::setw(13) << solp->ref_edge_t_buf.v_qp_u_buf[i][j][k][0]
                    << ' ' << std::setw(13) << solp->ref_edge_t_buf.v_qp_u_buf[i][j][k][1]
                    << '\n';
            }
        }
    }

    out << "QP_EDGE_L_E_BUF: " << '\n';
    for (int i = 0; i < solp->ref_edge_l_buf.nsub; i++)
    {
        for (int k = 0; k < solp->N_GL_1D; k++)
        {
            out << ' ' << std::setw(13) << solp->ref_edge_l_buf.qp_e_buf[i][k].x
                << ' ' << std::setw(13) << solp->ref_edge_l_buf.qp_e_buf[i][k].y
                << ' ' << std::setw(13) << solp->ref_edge_l_buf.qp_e_buf[i][k].w
                << '\n';
        }
    }

    out << "QP_EDGE_L_U_BUF: " << '\n';
    for (int i = 0; i < solp->ref_edge_l_buf.nsub; i++)
    {
        for (int k = 0; k < solp->N_GL_1D; k++)
        {
            out << ' ' << std::setw(13) << solp->ref_edge_l_buf.qp_u_buf[i][k].x
                << ' ' << std::setw(13) << solp->ref_edge_l_buf.qp_u_buf[i][k].y
                << ' ' << std::setw(13) << solp->ref_edge_l_buf.qp_u_buf[i][k].w
                << '\n';
        }
    }

    out << "V_QP_EDGE_L_E_BUF: " << '\n';
    for (int i = 0; i < solp->ref_edge_l_buf.nsub; i++)
    {
        for (int j = 0; j < solp->N_dof_loc; j++)
        {
            for (int k = 0; k < solp->N_GL_1D; k++)
            {
                out << ' ' << std::setw(6) << i
                    << ' ' << std::setw(6) << j
                    << ' ' << std::setw(6) << k
                    << ' ' << std::setw(13) << solp->ref_edge_l_buf.v_qp_e_buf[i][j][k][0]
                    << ' ' << std::setw(13) << solp->ref_edge_l_buf.v_qp_e_buf[i][j][k][1]
                    << '\n';
            }
        }
    }

    out << "V_QP_EDGE_L_U_BUF: " << '\n';
    for (int i = 0; i < solp->ref_edge_l_buf.nsub; i++)
    {
        for (int j = 0; j < solp->N_dof_loc; j++)
        {
            for (int k = 0; k < solp->N_GL_1D; k++)
            {
                out << ' ' << std::setw(6) << i
                    << ' ' << std::setw(6) << j
                    << ' ' << std::setw(6) << k
                    << ' ' << std::setw(13) << solp->ref_edge_l_buf.v_qp_u_buf[i][j][k][0]
                    << ' ' << std::setw(13) << solp->ref_edge_l_buf.v_qp_u_buf[i][j][k][1]
                    << '\n';
            }
        }
    }
#endif
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
        out << ' ' << std::setw(12) << solp->elem_e_buf[i].p0.x
            << ' ' << std::setw(12) << solp->elem_e_buf[i].p1.x
            << ' ' << std::setw(12) << solp->elem_e_buf[i].p0.y
            << ' ' << std::setw(12) << solp->elem_e_buf[i].p1.y
            << '\n';
        out << "QUAD_E: " << '\n';
        for (int j = 0; j < solp->ref_elem_buf.nsub; j++)
        {
            for (int k = 0; k < solp->N_GL_2D; k++)
            {
                out << ' ' << std::setw(6) << i
                    << ' ' << std::setw(12) << solp->elem_e_buf[i].qr.qp_buf[j][k].x
                    << ' ' << std::setw(12) << solp->elem_e_buf[i].qr.qp_buf[j][k].y
                    << ' ' << std::setw(12) << solp->elem_e_buf[i].qr.qp_buf[j][k].w
                    << '\n';
            }
        }
        out << "QUAD_U: " << '\n';
        for (int j = 0; j < solp->ref_elem_buf.nsub; j++)
        {
            for (int k = 0; k < solp->N_GL_2D; k++)
            {
                out << ' ' << std::setw(6) << solp->elem_u_buf[i].qr.par_buf[j]
                    << ' ' << std::setw(12) << solp->elem_u_buf[i].qr.qp_buf[j][k].x
                    << ' ' << std::setw(12) << solp->elem_u_buf[i].qr.qp_buf[j][k].y
                    << ' ' << std::setw(12) << solp->elem_u_buf[i].qr.qp_buf[j][k].w
                    << '\n';
            }
        }
    }
    out.close();
}
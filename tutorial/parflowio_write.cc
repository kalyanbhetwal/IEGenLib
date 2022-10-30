#include "iegenlib.h"
#include <utility>
#include <fstream>
#include <iostream>

using iegenlib::Computation;
using namespace std;

int main(int argc, char **argv){


    vector < pair<string, string> > dataReads;
    vector < pair<string, string> > dataWrites;
    Computation parflowio_w;

    parflowio_w.addDataSpace("nsg", "int");
    parflowio_w.addDataSpace("x", "int");
    parflowio_w.addDataSpace("y", "int");
    parflowio_w.addDataSpace("z", "int");
    parflowio_w.addDataSpace("nx", "int");
    parflowio_w.addDataSpace("ny", "int");
    parflowio_w.addDataSpace("nz", "int");
    parflowio_w.addDataSpace("rx", "int");
    parflowio_w.addDataSpace("ry", "int");
    parflowio_w.addDataSpace("rz", "int");
    parflowio_w.addDataSpace("errcheck", "int");
    parflowio_w.addDataSpace("x_overlap", "int");
    parflowio_w.addDataSpace("clip_x", "int");
    parflowio_w.addDataSpace("extent_x", "int");
    parflowio_w.addDataSpace("qq", "int");
    parflowio_w.addDataSpace("tmp", "uint64_t");
    parflowio_w.addDataSpace("buf", "uint64_t");

    parflowio_w.addDataSpace("m_p", "int");
    parflowio_w.addDataSpace("m_q", "int");
    parflowio_w.addDataSpace("m_r", "int");

    parflowio_w.addDataSpace("fp", "file*");
    parflowio_w.addDataSpace("m_X", "double");
    parflowio_w.addDataSpace("m_Y", "double");
    parflowio_w.addDataSpace("m_Z", "double");
    parflowio_w.addDataSpace("m_nx", "int");
    parflowio_w.addDataSpace("m_ny", "int");
    parflowio_w.addDataSpace("m_nz", "int");
    parflowio_w.addDataSpace("m_dX", "double");
    parflowio_w.addDataSpace("m_dY", "double");
    parflowio_w.addDataSpace("m_dZ", "double");
    parflowio_w.addDataSpace("m_numSubgrids", "double");
    parflowio_w.addDataSpace("max_x_extent", "int");

    parflowio_w.addDataSpace("byte_offsets", "long*");
    parflowio_w.addDataSpace("sg_count", "long long");
    parflowio_w.addDataSpace("x_extent", "int");

    /*
    parflowio_w.addDataSpace("nsg_x", "int");
    parflowio_w.addDataSpace("nsg_y", "int");
    parflowio_w.addDataSpace("nsg_z", "int");
    */



    Stmt s0("std::FILE* fp = std::fopen(filename.c_str(), 'wb');",
            "{[0]}",
            "{[0]->[0]}",
            {{"filename", "{[0]->[0]}"}},
            {{"fp", "{[0]->[0]}"}});

    parflowio_w.addStmt(&s0);


    Stmt s1("m_numSubgrids = m_p * m_q * m_r;    WRITEDOUBLE(m_X,fp);    WRITEDOUBLE(m_Y,fp);    WRITEDOUBLE(m_Z,fp);    WRITEINT(m_nx,fp);    WRITEINT(m_ny,fp);    WRITEINT(m_nz,fp);    WRITEDOUBLE(m_dX,fp);    WRITEDOUBLE(m_dY,fp);    WRITEDOUBLE(m_dZ,fp);    WRITEINT(m_numSubgrids,fp);max_x_extent =calcExtent(m_nx,m_p,0);",
            "{[0]}",
            "{[0]->[1]}",
            {{"m_p", "{[0]->[1]}"},
             {"m_q", "{[0]->[1]}"},
             {"m_r", "{[0]->[1]}"},
             {"m_X", "{[0]->[1]}"},
             {"m_Y", "{[0]->[1]}"},
             {"m_Z", "{[0]->[1]}"},
             {"m_nx", "{[0]->[1]}"},
             {"m_ny", "{[0]->[1]}"},
             {"m_nz", "{[0]->[1]}"},
             {"m_dX", "{[0]->[1]}"},
             {"m_dY", "{[0]->[1]}"},
             {"m_dZ", "{[0]->[1]}"},
             {"m_numSubgrids", "{[0]->[1]}"}},
            {
                {"m_numSubgrids", "{[0]->[1]}"},
                {"max_x_extent", "{[0]->[1]}"},
                {"fp", "{[0]->[1]}"},
             });

    parflowio_w.addStmt(&s1);


    /// how to add this statement
    ///     std::vector<double> writeBuf(max_x_extent);

    /// also these three statements
/*
    int nsg=0;
    byte_offsets[0] = 0;
    long long sg_count = 1;
*/


    Stmt s2("nsg=0; byte_offsets[0]=0; sg_count=1;",
            "{[0]}",
            "{[0]->[2]}",
            {{"filename", "{[0]->[2]}"}},
            {{"fp", "{[0]->[2]}"}});

    parflowio_w.addStmt(&s2);

    /*
     *     for(int nsg_z=0;nsg_z<m_r;nsg_z++){
        for(int nsg_y=0;nsg_y<m_q;nsg_y++) {
            for (int nsg_x = 0;nsg_x<m_p;nsg_x++)
     */



    Stmt s3("x = m_X + calcOffset(m_nx,m_p,nsg_x);y = m_Y + calcOffset(m_ny,m_q,nsg_y);z = m_Z + calcOffset(m_nz,m_r,nsg_z);// x,y,z of lower lefthand cornerWRITEINT(x, fp);WRITEINT(y, fp);WRITEINT(z, fp);x_extent =calcExtent(m_nx,m_p,nsg_x);WRITEINT(x_extent, fp);WRITEINT(calcExtent(m_ny,m_q,nsg_y), fp);WRITEINT(calcExtent(m_nz,m_r,nsg_z), fp);WRITEINT(1, fp);WRITEINT(1, fp);WRITEINT(1, fp);",
            "{[nsg_z, nsg_y, nsg_x]: 0<= nsg_z< m_r &&  0<= nsg_y< m_q &&  0<= nsg_x< m_p  }",
            "{[nsg_z, nsg_y, nsg_x]->[0,nsg_z,0,nsg_y,0,nsg_x,0]}",
            {
                        {"m_X", "{[nsg_z, nsg_y, nsg_x]->[0]}"},
                        {"m_Y", "{[nsg_z, nsg_y, nsg_x]->[0]}"},
                        {"m_Z", "{[nsg_z, nsg_y, nsg_x]->[0]}"},
                        {"m_nx", "{[nsg_z, nsg_y, nsg_x]->[0]}"},
                        {"m_p", "{[nsg_z, nsg_y, nsg_x]->[0]}"},
                        {"nsg_x", "{[nsg_z, nsg_y, nsg_x]->[0]}"},
                        {"m_ny", "{[nsg_z, nsg_y, nsg_x]->[0]}"},
                        {"m_q", "{[nsg_z, nsg_y, nsg_x]->[0]}"},
                        {"nsg_y", "{[nsg_z, nsg_y, nsg_x]->[0]}"},
                        {"m_nz", "{[nsg_z, nsg_y, nsg_x]->[0]}"},
                        {"m_r", "{[nsg_z, nsg_y, nsg_x]->[0]}"},
                        {"nsg_z", "{[nsg_z, nsg_y, nsg_x]->[0]}"},
                        },
            {{"fp", "{[nsg_z, nsg_y, nsg_x]->[0]}"}});

    parflowio_w.addStmt(&s3);


    /*
     *   for(iz=calcOffset(m_nz,m_r,nsg_z); iz < calcOffset(m_nz,m_r,nsg_z+1);iz++){
                    for(iy=calcOffset(m_ny,m_q,nsg_y); iy < calcOffset(m_ny,m_q,nsg_y+1);iy++){
                     uint64_t* buf = (uint64_t*)&(m_data[iz*m_nx*m_ny+iy*m_nx+calcOffset(m_nx,m_p,nsg_x)]);
     */


    Stmt s4("buf = (uint64_t*)&(m_data[iz*m_nx*m_ny+iy*m_nx+calcOffset(m_nx,m_p,nsg_x)]);",
            "{[nsg_z, nsg_y, nsg_x, iz, iy]: 0< nsg_z< m_r &&  0< nsg_y< m_q &&  0< nsg_x< m_p  && calcOffset(m_nz,m_r,nsg_z)  <=iz< calcOffset(m_nz,m_r,nsg_z+1) && calcOffset(m_ny,m_q,nsg_y) <=iz< calcOffset(m_ny,m_q,nsg_y+1) }",
            "{[nsg_z, nsg_y, nsg_x, iz, iy]->[0,nsg_z,0,nsg_y,0,nsg_x,0,iz,0,iy,0]}",
            {
                {"m_data", "{[nsg_z, nsg_y, nsg_x, iz, iy]->[0]}"},
             },
            {{"buf", "{[nsg_z, nsg_y, nsg_x, iz, iy]->[0]}"}});

    parflowio_w.addStmt(&s4);


    /*
     *  for(j=0;j<x_extent;j++){
            uint64_t tmp = buf[j];
            tmp = bswap64(tmp);
            writeBuf[j] = *(double*)(&tmp);
        }

     */

    Stmt s5("  tmp = buf[j]; tmp = bswap64(tmp); writeBuf[j] = *(double*)(&tmp);",
            "{[nsg_z, nsg_y, nsg_x, iz, iy,0,j,0]: 0< nsg_z< m_r &&  0< nsg_y< m_q &&  0< nsg_x< m_p   && calcOffset(m_nz,m_r,nsg_z)  <=iz< calcOffset(m_nz,m_r,nsg_z+1) && calcOffset(m_ny,m_q,nsg_y) <=iz< calcOffset(m_ny,m_q,nsg_y+1) && 0<=j<x_extent }",
            "{[nsg_z, nsg_y, nsg_x, iz, iy,j]->[0,nsg_z,0,nsg_y,0,nsg_x,0,iz,0,iy,1,j,0]}",
            {
                    {"buf", "{[nsg_z, nsg_y, nsg_x, iz, iy,j]->[j]}"},
                    {"tmp", "{[nsg_z, nsg_y, nsg_x, iz, iy,j]->[0]}"},
            },
            {{"tmp", "{[nsg_z, nsg_y, nsg_x, iz, iy,j]->[0]}"},
             {"writeBuf", "{[nsg_z, nsg_y, nsg_x, iz, iy,j]->[j]}"}
            });

    parflowio_w.addStmt(&s5);



    Stmt s6("written = fwrite(writeBuf.data(),sizeof(double),x_extent,fp);",
            "{[nsg_z, nsg_y, nsg_x, iz, iy]: 0< nsg_z< m_r &&  0< nsg_y< m_q &&  0< nsg_x< m_p  && calcOffset(m_nz,m_r,nsg_z) <=iz< calcOffset(m_nz,m_r,nsg_z+1) && calcOffset(m_ny,m_q,nsg_y) <=iz< calcOffset(m_ny,m_q,nsg_y+1) }",
            "{[nsg_z, nsg_y, nsg_x, iz, iy]->[0,nsg_z,0,nsg_y,0,nsg_x,0,iz,0,iy,2]}",
            {
                    {"writeBuf", "{[nsg_z, nsg_y, nsg_x, iz, iy]->[0]}"},
                    {"x_extent", "{[nsg_z, nsg_y, nsg_x, iz, iy]->[0]}"},
            },
            {
                {"written", "{[nsg_z, nsg_y, nsg_x, iz, iy]->[0]}"},
                {"fp", "{[nsg_z, nsg_y, nsg_x, iz, iy]->[0]}"},
             });

    parflowio_w.addStmt(&s6);



    Stmt s7("byte_offsets[sg_count] = ftell(fp); sg_count++;",
            "{[nsg_z, nsg_y, nsg_x]: 0<= nsg_z< m_r &&  0<= nsg_y< m_q &&  0<= nsg_x< m_p  }",
            "{[nsg_z, nsg_y, nsg_x]->[0,nsg_z,0,nsg_y,0,nsg_x,1]}",
            {
                    {"fp", "{[nsg_z, nsg_y, nsg_x]->[0]}"},
                    {"sg_count", "{[nsg_z, nsg_y, nsg_x]->[0]}"}
            },
            {
                     {"byte_offsets", "{[nsg_z, nsg_y, nsg_x]->[sg_count]}"},
                     {"sg_count", "{[nsg_z, nsg_y, nsg_x]->[sg_count]}"}
            });

    parflowio_w.addStmt(&s7);



    Stmt s8("nsg++;",
            "{[nsg_z]: 0<= nsg_z< m_r}",
            "{[nsg_z]->[0,nsg_z,1]}",
            {
                    {"nsg", "{[nsg_z]->[0]}"}
            },
            {
                    {"nsg", "{[nsg_z]->[0]}"}
            });

    parflowio_w.addStmt(&s8);

    //Calling
    parflowio_w.finalize();
    std:: cout << parflowio_w.toDotString() << '\n';

}

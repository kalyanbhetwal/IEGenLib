#include "iegenlib.h"
#include <utility>
#include <fstream>
#include <iostream>

using iegenlib::Computation;
using namespace std;

int main(int argc, char **argv){


    vector < pair<string, string> > dataReads;
    vector < pair<string, string> > dataWrites;
    Computation parflowio;

    parflowio.addDataSpace("nsg","int");
    parflowio.addDataSpace("x","int");
    parflowio.addDataSpace("y","int");
    parflowio.addDataSpace("z","int");
    parflowio.addDataSpace("nx","int");
    parflowio.addDataSpace("ny","int");
    parflowio.addDataSpace("nz","int");
    parflowio.addDataSpace("rx","int");
    parflowio.addDataSpace("ry","int");
    parflowio.addDataSpace("rz","int");
    parflowio.addDataSpace("errcheck","int");
    parflowio.addDataSpace("x_overlap","int");
    parflowio.addDataSpace("clip_x","int");
    parflowio.addDataSpace("extent_x","int");
    parflowio.addDataSpace("qq","int");
    parflowio.addDataSpace("tmp","uint64_t");
    parflowio.addDataSpace("buf","uint64_t");
    parflowio.addDataSpace("m_nx","uint64_t");
    parflowio.addDataSpace("m_ny","uint64_t");


    Stmt s0("  READINT(x,m_fp,errcheck);READINT(y,m_fp,errcheck);READINT(z,m_fp,errcheck);READINT(nx,m_fp,errcheck);READINT(ny,m_fp,errcheck);READINT(nz,m_fp,errcheck);READINT(rx,m_fp,errcheck);READINT(ry,m_fp,errcheck);READINT(rz,m_fp,errcheck);",
            "{[nsg] : 0 <= nsg < m_numSubgrids}",
            "{[nsg]->[0, nsg, 0]}",
            {
                    {"m_fp", "{[nsg] -> [0]}"},
            },
            {
                    {"m_fp", "{[nsg] -> [0]}"},
                    {"x", "{[nsg] -> [0]}"},
                    {"y" ,"{[nsg]->[0]}"},
                    {"nx" ,"{[nsg]->[0]}"},
                    {"ny" ,"{[nsg]->[0]}"},
                    {"nz" ,"{[nsg]->[0]}"},
                    {"rz" ,"{[nsg]->[0]}"},
                    {"ry" ,"{[nsg]->[0]}"},
                    {"rz" ,"{[nsg]->[0]}"}
            });

    parflowio.addStmt(&s0);

// Statement 1
// long long qq = z*m_nx*m_ny + y*m_nx + x;
// long long k,i,j;

    Stmt s1("qq = z*m_nx*m_ny + y*m_nx + x",
            "{[nsg] : 0 <= nsg < m_numSubgrids}",
            "{[nsg]->[0, nsg, 1]}",
            {
                    {"m_ny", "{[nsg] -> [0]}"},
                    {"m_nx", "{[nsg] -> [0]}"},
                    {"x", "{[nsg] -> [0]}"},
                    {"y", "{[nsg] -> [0]}"},
                    {"z", "{[nsg] -> [0]}"}
            },
            {
                    {"qq" ,"{[nsg]->[0]}"}
            });

    parflowio.addStmt(&s1);



    // statement 2
    /*
      for (k=0; k<nz; k++){
        for(i=0;i<ny;i++){
          // read full "pencil"
          long long index = qq+k*m_nx*m_ny+i*m_nx;
          uint64_t* buf = (uint64_t*)&(m_data[index]);
          int read_count = fread(buf,8,nx,m_fp);
          if(read_count != nx){
              perror("Error Reading Data, File Ended Unexpectedly");
              return 1;
          }
    */

    Stmt s2("index = qq+k*m_nx*m_ny+i*m_nx; buf = (uint64_t*)&(m_data[index]);read_count = fread(buf,8,nx,m_fp);",
            "{[nsg,k,i] : 0 <= k < nz && 0<=i<ny && 0 <= nsg < m_numSubgrids}",
            "{[nsg,k,i]->[0, nsg,2 , k, 0,i,0 ]}",
            {
                    {"m_fp", "{[nsg] -> [0]}"},
                    {"m_data", "{[nsg,k,i] -> [0]}"},
                    {"qq", "{[nsg,k,i] -> [0]}"}
            },
            {   {"m_fp", "{[nsg] -> [0]}"},
                {"buf", "{[nsg,k,i] -> [0]}"}
            });

    parflowio.addStmt(&s2);

// statement
/*
    // handle byte order
    // uint64_t* buf = (uint64_t*)&(m_data[index]);
    for(j=0;j<nx;j++){
        uint64_t tmp = buf[j];
        tmp = bswap64(tmp);
        m_data[index+j] = *(double*)(&tmp);
*/

    Stmt s12(" tmp = buf[j];  tmp = bswap64(tmp);  m_data[index+j] = *(double*)(&tmp);",
             "{[nsg,k,i,j] :0 <= k < nz && 0<=i<ny && 0 <= nsg < m_numSubgrids && 0<=j<nx }",
             "{[nsg,k,i,j]->[0, nsg, 2, k, 0,i,1 ,j,0]}",
             {
                     {"buf","{[nsg,k,i,j]->[j]}" }

             },
             {
                     {"tmp", "{[nsg,k,i,j] -> [0]}"},
                     {"m_data", "{[nsg,k,i,j] -> [0]}"}
             });

    parflowio.addStmt(&s12);


    // Stmt sy("sum+=test(x,y,z)",
    //   "[z,y,x]: 0<=z<NZ && 0<=y<=NY && 0<=x<NX",
    //   "{[z,y,z]->[1, z,0,y,0,x,0]}",
    //   {
    //     {"m_fp", "{[0] -> [0]}"},
    //     {"nx", "{[0] -> [0]}"},
    //     {"ny", "{[0] -> [0]}"},
    //     {"nz", "{[0] -> [0]}"},
    //     {"x_overlap","{[0] -> [0]}" },
    //     {"y_overlap","{[0] -> [0]}" },

    //   },
    //   { });


    // parflowio.addStmt(&sy);



    //Calling
    parflowio.finalize();
    std:: cout << parflowio.toDotString()<<'\n';

    //cout << parflowio.codeGen();
    return 0;
}

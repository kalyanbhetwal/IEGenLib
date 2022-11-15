/*!
 * \file SSA_test.cc
 *
 * \brief Tests for the dominanceTree class.
 *
 * \date Started: 08/02/22
 *
 * \authors Kalyan Bhetwal
 *
 * Copyright (c) 2020, Boise State University <br>
 * All rights reserved. <br>
 * See ../../COPYING for details. <br>
 */
#include "SSA.h"
#include "Computation.h"
#include <gtest/gtest.h>
#include "set_relation/set_relation.h"
#include <utility>
#include <vector>
#include <map>
#include "code_gen/parser/parser.h"
#include "omega/Relation.h"
using namespace SSA;
using namespace std;


TEST(SSATest123, DominanceTreeTEST111){


    Computation * comp = new Computation();
    comp->addDataSpace("x", "int");

    //s1
    //    iegenlib::Set* s1 = new iegenlib::Set("{[0]}");
    comp->addStmt(new Stmt (
            "x = 1;",
            "{[0]}",
            "{[0]->[0]}",
            {},
            {{"x", "{[0]->[0]}"}}
    ));

    //s2
    // iegenlib::Set* s2 = new iegenlib::Set("{[1]}");

    comp->addStmt(new Stmt (
            "x = 2;",
            "{[0]}",
            "{[0]->[1]}",
            {{"x", "{[0]->[0]}"}},
            {{"x", "{[0]->[0]}"}}
    ));

//    //s3
//    //iegenlib::Set* s3 = new iegenlib::Set("{[2,t,0]: 0<=t<M}");
//
//    comp->addStmt(new Stmt (
//            "x=2;",
//            "{[t]:0<=t<M }",
//            "{[t]->[2,t,0]}",
//            {},
//            {{"x", "{[t]->[t]}"}}
//    ));

    //s4
    //iegenlib::Set* s4 = new iegenlib::Set("{[2,t,1,p,0]: 0<=t<M && p> 10 }");

    comp->addStmt(new Stmt (
            "x = 2;",
            "{[t,p]:0<=t<M && p> 10}",
            "{[t,p]->[2,t,0,p,0]}",
            {{"x", "{[t,p]->[t,p]}"}},
            {{"x", "{[t,p]->[t,p]}"}}
    ));

    //s5
    //iegenlib::Set* s5 = new iegenlib::Set("{[2,t,1,p,1,q,0]:0<=t<M && p>10 && q>10}");

    comp->addStmt(new Stmt (
            "x = 2;",
            "{[t,p,q]:0<=t<M && p>10 && q>10}",
            "{[t,p,q]->[2,t,0,p,1,q,0]}",
            {{"x", "{[t,p,q]->[t,p,q]}"}},
            {{"x", "{[t,p,q]->[t,p,q]}"}}
    ));

    //s6
    //iegenlib::Set* s6 = new iegenlib::Set("{[2,t,1,p,1,n,0]:0<=t<M && p>10 && n<=10 }");

    comp->addStmt(new Stmt (
            "x = 2;",
            "{[t,p,v]:0<=t<M && p>10 && v<=10}",
            "{[t,p,v]->[2,t,0,p,1,v,0]}",
            {{"x", "{[t,p,v]->[t,p,v]}"}},
            {{"x", "{[t,p,v]->[t,p,v]}"}}
    ));

    //s7
    //iegenlib::Set* s7 = new iegenlib::Set("{[2,t,1,p,2]: 0<=t<M && p>10}");

    comp->addStmt(new Stmt (
            "x = 2;",
            "{[t,p]:0<=t<M && p>10}",
            "{[t,p]->[2,t,0,p,2]}",
            {{"x", "{[t,p]->[t,p]}"}},
            {{"x", "{[t,p]->[t,p]}"}}
    ));

    //s8
    //iegenlib::Set* s8 = new iegenlib::Set("{[2,t,2,m,0]: 0<=t<M && m<=10}");


    comp->addStmt(new Stmt (
            "x = 2;",
            "{[t,u]:0<=t<M && u<=10}",
            "{[t,u]->[2,t,0,u,0]}",
            {{"x", "{[t,u]->[t,u]}"}},
            {{"x", "{[t,u]->[t,u]}"}}
    ));

    //s9
    //iegenlib::Set* s9 = new iegenlib::Set("{[2,t,3]: 0<=t<M}");

    comp->addStmt(new Stmt (
            "x = 2;",
            "{[t]:0<=t<M}",
            "{[t]->[2,t,1]}",
            {{"x", "{[t]->[t]}"}},
            {{"x", "{[t]->[t]}"}}
    ));

    //s10
    //iegenlib::Set* s10 = new iegenlib::Set("{[2,t,4,s,0,r,0]: 0<=t<M && 0<=s<S && r>10}");

    comp->addStmt(new Stmt (
            "x = 2;",
            "{[t,s,r]:0<=t<M && s>10 && r>20}",
            "{[t,s,r]->[2,t,3,s,0,r,0]}",
            {{"x", "{[t,s,r]->[t,s,r]}"}},
            {{"x", "{[t,s,r]->[t,s,r]}"}}
    ));
    //s11
    //iegenlib::Set* s11 = new iegenlib::Set("{[2,t,5]: 0<=t<M}");

    comp->addStmt(new Stmt (
            "x = 2;",
            "{[t]:0<=t<M}",
            "{[t]->[2,t,4]}",
            {{"x", "{[t]->[t]}"}},
            {{"x", "{[t]->[t]}"}}
    ));
    //s12
    //iegenlib::Set* s12 = new iegenlib::Set("{[3]}");

    comp->addStmt(new Stmt (
            "x = 2;",
            "{[0]}",
            "{[0]->[3]}",
            {{"x", "{[0]->[0]}"}},
            {{"x", "{[0]->[0]}"}}
    ));
  //SSA::generateSSA(comp);
    comp->finalize();
    std:: cout << comp->toDotString();
    EXPECT_EQ(1,1);

}

TEST(SSATest,AllPossiblePathsTest){
    Computation * comp = new Computation();
    comp->addDataSpace("x", "int");
    //s1
    //    iegenlib::Set* s1 = new iegenlib::Set("{[0]}");
    comp->addStmt(new Stmt (
            "x = 1;",
            "{[0]}",
            "{[0]->[0]}",
            {},
            {{"x", "{[0]->[0]}"}}
    ));
    //s2
    // iegenlib::Set* s2 = new iegenlib::Set("{[1]}");
    comp->addStmt(new Stmt (
            "x = 2;",
            "{[0]}",
            "{[0]->[1]}",
            {{"x", "{[0]->[0]}"}},
            {{"x", "{[0]->[0]}"}}
    ));

    //s4
    //iegenlib::Set* s4 = new iegenlib::Set("{[2,t,1,p,0]: 0<=t<M && p> 10 }");
    comp->addStmt(new Stmt (
            "x = 2;",
            "{[t,p]:0<=t<M && p> 10}",
            "{[t,p]->[2,t,0,p,0]}",
            {{"x", "{[t,p]->[t,p]}"}},
            {{"x", "{[t,p]->[t,p]}"}}
    ));
    //s5
    //iegenlib::Set* s5 = new iegenlib::Set("{[2,t,1,p,1,q,0]:0<=t<M && p>10 && q>10}");
    comp->addStmt(new Stmt (
            "x = 2;",
            "{[t,p,q]:0<=t<M && p>10 && q>10}",
            "{[t,p,q]->[2,t,0,p,1,q,0]}",
            {{"x", "{[t,p,q]->[t,p,q]}"}},
            {{"x", "{[t,p,q]->[t,p,q]}"}}
    ));

    //s6
    //iegenlib::Set* s6 = new iegenlib::Set("{[2,t,1,p,1,n,0]:0<=t<M && p>10 && n<=10 }");
    comp->addStmt(new Stmt (
            "x = 2;",
            "{[t,p,v]:0<=t<M && p>10 && v<=10}",
            "{[t,p,v]->[2,t,0,p,1,v,0]}",
            {{"x", "{[t,p,v]->[t,p,v]}"}},
            {{"x", "{[t,p,v]->[t,p,v]}"}}
    ));

    //s7
    //iegenlib::Set* s7 = new iegenlib::Set("{[2,t,1,p,2]: 0<=t<M && p>10}");
    comp->addStmt(new Stmt (
            "x = 2;",
            "{[t,p]:0<=t<M && p>10}",
            "{[t,p]->[2,t,0,p,2]}",
            {{"x", "{[t,p]->[t,p]}"}},
            {{"x", "{[t,p]->[t,p]}"}}
    ));

    //s8
    //iegenlib::Set* s8 = new iegenlib::Set("{[2,t,2,m,0]: 0<=t<M && m<=10}");
    comp->addStmt(new Stmt (
            "x = 2;",
            "{[t,u]:0<=t<M && u<=10}",
            "{[t,u]->[2,t,0,u,0]}",
            {{"x", "{[t,u]->[t,u]}"}},
            {{"x", "{[t,u]->[t,u]}"}}
    ));

    //s9
    //iegenlib::Set* s9 = new iegenlib::Set("{[2,t,3]: 0<=t<M}");
    comp->addStmt(new Stmt (
            "x = 2;",
            "{[t]:0<=t<M}",
            "{[t]->[2,t,1]}",
            {{"x", "{[t]->[t]}"}},
            {{"x", "{[t]->[t]}"}}
    ));

    //s10
    //iegenlib::Set* s10 = new iegenlib::Set("{[2,t,4,s,0,r,0]: 0<=t<M && 0<=s<S && r>10}");
    comp->addStmt(new Stmt (
            "x = 2;",
            "{[t,s,r]:0<=t<M && s>10 && r>20}",
            "{[t,s,r]->[2,t,3,s,0,r,0]}",
            {{"x", "{[t,s,r]->[t,s,r]}"}},
            {{"x", "{[t,s,r]->[t,s,r]}"}}
    ));
    //s11
    //iegenlib::Set* s11 = new iegenlib::Set("{[2,t,5]: 0<=t<M}");
    comp->addStmt(new Stmt (
            "x = 2;",
            "{[t]:0<=t<M}",
            "{[t]->[2,t,4]}",
            {{"x", "{[t]->[t]}"}},
            {{"x", "{[t]->[t]}"}}
    ));
    //s12
    //iegenlib::Set* s12 = new iegenlib::Set("{[3]}");
    comp->addStmt(new Stmt (
            "x = 2;",
            "{[0]}",
            "{[0]->[3]}",
            {{"x", "{[0]->[0]}"}},
            {{"x", "{[0]->[0]}"}}
    ));
    Node * node = createScheduleTree(comp);
    node->calc_all_backward_paths();
    ASSERT_EQ(SSA::Member::possiblePaths.size(),11);

    ASSERT_EQ(SSA::Member::possiblePaths[0].size(),11);
}


TEST(SSATest,AllPossiblePathsBug){
    Computation * comp = new Computation();
    comp->addDataSpace("x", "int");
    //s1
    //    iegenlib::Set* s1 = new iegenlib::Set("{[0]}");
    comp->addStmt(new Stmt (
            "x = 1;",
            "{[0]}",
            "{[0]->[0]}",
            {},
            {{"x", "{[0]->[0]}"}}
    ));
    //s2
    // iegenlib::Set* s2 = new iegenlib::Set("{[1]}");
    comp->addStmt(new Stmt (
            "x = 2;",
            "{[0]}",
            "{[0]->[1]}",
            {{"x", "{[0]->[0]}"}},
            {{"x", "{[0]->[0]}"}}
    ));

    //s4
    //iegenlib::Set* s4 = new iegenlib::Set("{[2,t,1,p,0]: 0<=t<M && p> 10 }");
    comp->addStmt(new Stmt (
            "x = 2;",
            "{[t,p]:0<=t<M && p> 10}",
            "{[t,p]->[2,t,0,p,0]}",
            {{"x", "{[t,p]->[t,p]}"}},
            {{"x", "{[t,p]->[t,p]}"}}
    ));
    //s5
    //iegenlib::Set* s5 = new iegenlib::Set("{[2,t,1,p,1,q,0]:0<=t<M && p>10 && q>10}");
    comp->addStmt(new Stmt (
            "x = 2;",
            "{[t,p,q]:0<=t<M && p>10 && q>10}",
            "{[t,p,q]->[2,t,0,p,1,q,0]}",
            {{"x", "{[t,p,q]->[t,p,q]}"}},
            {{"x", "{[t,p,q]->[t,p,q]}"}}
    ));

    //s6
    //iegenlib::Set* s6 = new iegenlib::Set("{[2,t,1,p,1,n,0]:0<=t<M && p>10 && n<=10 }");
    comp->addStmt(new Stmt (
            "x = 2;",
            "{[t,p,v]:0<=t<M && p>10 && v<=10}",
            "{[t,p,v]->[2,t,0,p,1,v,0]}",
            {{"x", "{[t,p,v]->[t,p,v]}"}},
            {{"x", "{[t,p,v]->[t,p,v]}"}}
    ));

    //s7
    //iegenlib::Set* s7 = new iegenlib::Set("{[2,t,1,p,2]: 0<=t<M && p>10}");
    comp->addStmt(new Stmt (
            "x = 2;",
            "{[t,p]:0<=t<M && p>10}",
            "{[t,p]->[2,t,0,p,2]}",
            {{"x", "{[t,p]->[t,p]}"}},
            {{"x", "{[t,p]->[t,p]}"}}
    ));

    //s8
    //iegenlib::Set* s8 = new iegenlib::Set("{[2,t,2,m,0]: 0<=t<M && m<=10}");
    comp->addStmt(new Stmt (
            "x = 2;",
            "{[t,u]:0<=t<M && u<=10}",
            "{[t,u]->[2,t,0,u,0]}",
            {{"x", "{[t,u]->[t,u]}"}},
            {{"x", "{[t,u]->[t,u]}"}}
    ));

    //s9
    //iegenlib::Set* s9 = new iegenlib::Set("{[2,t,3]: 0<=t<M}");
    comp->addStmt(new Stmt (
            "x = 2;",
            "{[t]:0<=t<M}",
            "{[t]->[2,t,1]}",
            {{"x", "{[t]->[t]}"}},
            {{"x", "{[t]->[t]}"}}
    ));

    //s10
    //iegenlib::Set* s10 = new iegenlib::Set("{[2,t,4,s,0,r,0]: 0<=t<M && 0<=s<S && r>10}");
    comp->addStmt(new Stmt (
            "x = 2;",
            "{[t,s,r]:0<=t<M && s>10 && r>20}",
            "{[t,s,r]->[2,t,3,s,0,r,0]}",
            {{"x", "{[t,s,r]->[t,s,r]}"}},
            {{"x", "{[t,s,r]->[t,s,r]}"}}
    ));
    //s11
    //iegenlib::Set* s11 = new iegenlib::Set("{[2,t,5]: 0<=t<M}");
    comp->addStmt(new Stmt (
            "x = 2;",
            "{[t]:0<=t<M}",
            "{[t]->[2,t,4]}",
            {{"x", "{[t]->[t]}"}},
            {{"x", "{[t]->[t]}"}}
    ));
    //s12
    //iegenlib::Set* s12 = new iegenlib::Set("{[3]}");
    comp->addStmt(new Stmt (
            "x = 2;",
            "{[0]}",
            "{[0]->[3]}",
            {{"x", "{[0]->[0]}"}},
            {{"x", "{[0]->[0]}"}}
    ));
    Node * node = createScheduleTree(comp);
    node->calc_all_backward_paths();

    for(std::map<Stmt*, std::vector<Stmt*>>::const_iterator it = SSA::Member::possiblePaths.begin(); it != SSA::Member::possiblePaths.end(); ++it)
    {
       std::cout <<"Stmt " << it->first->getExecutionSchedule()->prettyPrintString()<<std::endl;
        for(auto v: it->second){
            std::cout<<"Back paths " << v->getExecutionSchedule()->prettyPrintString()<<std::endl;
        }
       std::cout<<"----------------------------------------"<<std::endl;
    }
    ASSERT_EQ(SSA::Member::possiblePaths.size(),11);

    ASSERT_EQ(SSA::Member::possiblePaths[0].size(),11);
}

TEST(SSATest, SSARenaming1) {

    Computation* computation = new Computation();

    computation->addParameter("foo", "int");
    computation->addDataSpace("bar", "int");

    Stmt* s1 = new Stmt("bar = foo;",
                        "{[0]}",
                        "{[0]->[0]}",
                        {{"foo", "{[0]->[0]}"}},
                        {{"bar", "{[0]->[0]}"}}
    );
    computation->addStmt(s1);

    Stmt* s2 =new Stmt("foo = bar + 1",
                       "{[0]}", "{[0]->[1]}",
                       {{"bar", "{[0]->[0]}"}},
                       {{"foo", "{[0]->[0]}"}}
    );
    computation->addStmt(s2);

    computation->finalize();
    //std:: cout << codeGen;
    EXPECT_EQ("a","a");
    delete computation;
}
TEST(SSATest, MTTKRP){

    // this one is dense
    //void mttkrp(int I,int J, int K, int R,double *X,
    //               double *A, double *B, double *C) {
    // for (i = 0; i < I; i++)
    //   for (j = 0; j < J; j++)
    //     for (k = 0; k < K; k++)
    //       for (r = 0; r < R; r++)
    //         A[i,r] += X[i,j,k]*B[j,r]*C[k,r];
    ///
    vector < pair<string, string> > dataReads;
    vector < pair<string, string> > dataWrites;
    Computation mttkrp;
    mttkrp.addDataSpace("X","double*");
    mttkrp.addDataSpace("A","double*");
    mttkrp.addDataSpace("B","double*");
    mttkrp.addDataSpace("C","double*");
    Stmt *s0 = new Stmt("x=1", "{[i,j,k,r] : 0 <= i < I and 0<=j<J and 0<=k<K and 0<=r<R}",
                       "{[i,j,k,r]->[0,i,0,j,0,k,0,r,0]}",
                       {},
                       { {"A", "{[i,k,l,j]->[i,j]}"},
                         {"B", "{[i,k,l,j]->[i,k,l]}"},
                         {"D", "{[i,k,l,j]->[l,j]}"},
                         {"C", "{[i,k,l,j]->[k,j]}"}});
    mttkrp.addStmt(s0);
    Stmt *s1 = new Stmt("A(i,r) += X(i,j,k)*B(j,r)*C(k,r)",
                        "{[i,j,k,r] : 0 <= i < I and 0<=j<J and 0<=k<K and 0<=r<R}",
                        "{[i,j,k,r]->[1,i,0,j,0,k,0,r,0]}",
                        {
                                // data reads
                                {"A", "{[i,k,l,j]->[i,j]}"},
                                {"B", "{[i,k,l,j]->[i,k,l]}"},
                                {"D", "{[i,k,l,j]->[l,j]}"},
                                {"C", "{[i,k,l,j]->[k,j]}"},
                        },
                        {
                                // data writes
                                {"A", "{[i,k,l,j]->[i,j]}"},
                        });

    mttkrp.addStmt(s1);
//    Stmt *s2 = new Stmt("x=2", "{[0]}", "{[0]->[2]}",{},{});
//    mttkrp.addStmt(s2);

    mttkrp.finalize();
    std:: cout << mttkrp.toDotString();

    EXPECT_EQ("1","1");

}

TEST(SSATest, DISABLED_Parflowio){


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
                    {"rz" ,"{[nsg]->[0]}"},
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
                    {"qq", "{[nsg,k,i] -> [0]}"},
                    {"m_fp", "{[nsg] -> [0]}"},
                    {"m_data", "{[nsg,k,i] -> [0]}"},
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
                     {"m_data", "{[nsg,k,i,j] -> [0]}"},
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
    EXPECT_EQ("1","1");

}

TEST(SSATest, SSARenameExample) {

    Computation * comp = new Computation();
    comp->addDataSpace("x", "int");

    comp->addStmt(new Stmt (
            "x(i) = 1;",
            "{[i]:0 <=i<N}",
            "{[i]->[0,i,0]}",
            {},
            {{"x", "{[i]->[i]}"}}
    ));


    comp->addStmt(new Stmt (
            "x(i) = 2;",
            "{[i]:0 <=i<N}",
            "{[i]->[1,i,0]}",
            {},
            {{"x", "{[i]->[i]}"}}
    ));

    comp->addStmt(new Stmt (
            "",
            "{[i]:0 <=i<N}",
            "{[i]->[1,i,1]}",
            {{"x", "{[i]->[i]}"}},
            {}
    ));

    comp->finalize();
    //std::cout <<"------------------------------------------"<<std::endl;
    std:: cout << comp->toDotString();

}


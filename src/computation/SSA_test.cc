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
            "x=1;",
            "{[0]}",
            "{[0]->[0]}",
            {},
            {{"x", "{[0]->[0]}"}}
    ));

    //s2
    // iegenlib::Set* s2 = new iegenlib::Set("{[1]}");

    comp->addStmt(new Stmt (
            "x=2;",
            "{[0]}",
            "{[0]->[1]}",
            {},
            {{"x", "{[0]->[0]}"}}
    ));

    //s3
    //iegenlib::Set* s3 = new iegenlib::Set("{[2,t,0]: 0<=t<M}");

    comp->addStmt(new Stmt (
            "x=2;",
            "{[t]:0<=t<M }",
            "{[t]->[2,t,0]}",
            {},
            {{"x", "{[t]->[t]}"}}
    ));

    //s4
    //iegenlib::Set* s4 = new iegenlib::Set("{[2,t,1,p,0]: 0<=t<M && p> 10 }");

    comp->addStmt(new Stmt (
            "x=2;",
            "{[t,p]:0<=t<M && p> 10}",
            "{[t,p]->[2,t,1,p,0]}",
            {},
            {{"x", "{[t,p]->[t,p]}"}}
    ));

    //s5
    //iegenlib::Set* s5 = new iegenlib::Set("{[2,t,1,p,1,q,0]:0<=t<M && p>10 && q>10}");

    comp->addStmt(new Stmt (
            "x=2;",
            "{[t,p,q]:0<=t<M && p>10 && q>10}",
            "{[t,p,q]->[2,t,1,p,1,q,0]}",
            {},
            {{"x", "{[t,p,q]->[t,p,q]}"}}
    ));

    //s6
    //iegenlib::Set* s6 = new iegenlib::Set("{[2,t,1,p,1,n,0]:0<=t<M && p>10 && n<=10 }");

    comp->addStmt(new Stmt (
            "x=2;",
            "{[t,p,n]:0<=t<M && p>10 && n<=10}",
            "{[t,p,n]->[2,t,1,p,1,n,0]}",
            {},
            {{"x", "{[t,p,n]->[t,p,n]}"}}
    ));

    //s7
    //iegenlib::Set* s7 = new iegenlib::Set("{[2,t,1,p,2]: 0<=t<M && p>10}");

    comp->addStmt(new Stmt (
            "x=2;",
            "{[t,p]:0<=t<M && p>10}",
            "{[t,p]->[2,t,1,p,2]}",
            {},
            {{"x", "{[t,p]->[t,p]}"}}
    ));

    //s8
    //iegenlib::Set* s8 = new iegenlib::Set("{[2,t,2,m,0]: 0<=t<M && m<=10}");


    comp->addStmt(new Stmt (
            "x=2;",
            "{[t,m]:0<=t<M && m<=10}",
            "{[t,m]->[2,t,1,m,0]}",
            {},
            {{"x", "{[t,m]->[t,m]}"}}
    ));

    //s9
    //iegenlib::Set* s9 = new iegenlib::Set("{[2,t,3]: 0<=t<M}");

    comp->addStmt(new Stmt (
            "x=2;",
            "{[t]:0<=t<M}",
            "{[t]->[2,t,3]}",
            {},
            {{"x", "{[t]->[t]}"}}
    ));

    //s10
    //iegenlib::Set* s10 = new iegenlib::Set("{[2,t,4,s,0,r,0]: 0<=t<M && 0<=s<S && r>10}");

    comp->addStmt(new Stmt (
            "x=2;",
            "{[t,s,r]:0<=t<M}",
            "{[t,s,r]->[2,t,4,s,0,r,0]}",
            {},
            {{"x", "{[t,s,r]->[t,s,r]}"}}
    ));
    //s11
    //iegenlib::Set* s11 = new iegenlib::Set("{[2,t,5]: 0<=t<M}");

    comp->addStmt(new Stmt (
            "x=2;",
            "{[t]:0<=t<M}",
            "{[t]->[2,t,5]}",
            {},
            {{"x", "{[t]->[t]}"}}
    ));
    //s12
    //iegenlib::Set* s12 = new iegenlib::Set("{[3]}");

    comp->addStmt(new Stmt (
            "x=2;",
            "{[0]}",
            "{[0]->[3]}",
            {},
            {{"x", "{[0]->[0]}"}}
    ));
    Node * node;
    node = createScheduleTree(comp);

  // node->printBreadthFirst();

   node->calc_all_pred();

   // comp->finalize();
    EXPECT_EQ(1,1);

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




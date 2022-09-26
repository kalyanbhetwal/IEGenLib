/*!
 * \file SSA.cpp
 *
 * \brief Implementation of the SSA class
 *
 * The SSA class is the class that gives us the dominance graph.
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
#include <vector>
#include <tuple>
#include <stack>
#include "set_relation/set_relation.h"
#include <iostream>
#include "Computation.h"
#include <utility>

using namespace SSA;
using namespace iegenlib;
DominanceTree::DominanceTree() {}

DominanceTree::~DominanceTree() {
}

int DominanceTree::push_Back(std::pair<int, iegenlib::Set*> data) {
    Node* node = new Node();
    node->parent= -1;
    node->children = {};
    node->data = data;
    node-> dominanceFrontier = {};
    nodes.push_back(*node);
    return(nodes.size()-1);
}


void DominanceTree::add_edge(int parent, int child) {
    nodes[parent].children.push_back(child);
    nodes[child].parent = parent;
}
// j is predecessor of i
void DominanceTree::add_predecessors(int i, int j){
    nodes[i].predecessors.push_back(j);
}
int DominanceTree::getVectorSize() {
    return nodes.size();
}
Set* DominanceTree::getElem(int i){
    return(nodes[i].data.second);

}
bool DominanceTree::isParent(int parent, int child) {
   for(auto v:  nodes[parent].children){
       if(v == child){
           return true;
       }
   }
    return  false;
}
std::vector<int> DominanceTree::getPredecessor(int i){
    return(nodes[i].predecessors);
}
//calculate Dominance Frontier
void DominanceTree::DFCal() {
    for(int i=0; i<this->nodes.size();i++){
        if( this->nodes[i].predecessors.size() > 1) {
           // std:: cout<< "pred for node " <<i <<'\n';
           //std::cout << "--------------------"<<'\n';
            for (int pred:  this->nodes[i].predecessors) {
                //std:: cout<< "pred "<< pred << " for node " <<i <<'\n';
                int runner = pred;
                std::vector<int> DF_runner =  this->nodes[runner].dominanceFrontier;
                while( runner !=  this->nodes[i].parent){
                    //check if element i isn't in the vector
                    if(std::find(DF_runner.begin(), DF_runner.end(), i) == DF_runner.end()) {
                        DF_runner.push_back(i);
                        //std:: cout << "DF  " << i<<'\n';
                        //looks like additional DF nodes are being added
                        // might have to remove the statement
                        this->nodes[runner].dominanceFrontier.push_back(i);
                    }
                    runner =   this->nodes[runner].parent;
                }
                std::vector<int> test =  this->nodes[i].dominanceFrontier;
                //std::cout << "runner " << runner <<'\n';
                if(std::find(test.begin(), test.end(), i) == test.end()){
                    test.push_back(runner);
                }

            }
        }
        //this->nodes[runner].dominanceFrontier.push_back(i);
    }
    return;
}

//test if two dominator tree are equivalent
bool DominanceTree::equivalent(DominanceTree dt) {

 for(int i=0;i<dt.nodes.size();i++) {
     if (dt.nodes[i].data.second->prettyPrintString() != this->nodes[i].data.second->prettyPrintString()) return false;

     std::vector<int> v1 = dt.nodes[i].children;
     std::vector<int> v2 = this->nodes[i].children;
     std::sort(v1.begin(), v1.end());
     std::sort(v2.begin(), v2.end());
     if(v1.size()) continue;

     if (v1 != v2) {
         return false;
     }
 }
    return true;

}
bool DominanceTree::predecessorEquivalent(DominanceTree dt){
for(int i=0;i<dt.getVectorSize();i++) {
    for(int j = 0; j<dt.nodes[i].predecessors.size();j++){
        if(this->nodes[i].predecessors.size()>0 && this->nodes[i].predecessors.size()== dt.nodes[i].predecessors.size()) {
            if (dt.nodes[i].predecessors[j] != this->nodes[i].predecessors[j]) return false;
            if (dt.nodes[i].predecessors[j] != this->nodes[i].predecessors[j]) return false;
        }
        else{ return  false;}
}
}
    return true;
}


bool SSA::isReverseDominator(iegenlib::Set * s1, iegenlib::Set * s2){
    if(s2->getArity()< s1->getArity()){
        return true;
    }
    while( s2->getArity() > s1->getArity()) {

        TupleDecl tl = s2->getTupleDecl();
        s2 = s2->projectOut(tl.getSize()-1);
    }
    if (s2->isSubset(s1)) return false;
    return true;

}

// function returns list of prefixes
std::vector<Set*> DominanceTree::getPrefixes(Set*s) {

    std::vector<Set*>v;
    Set* res = new Set(*s);
    TupleDecl tl = s->getTupleDecl();
    if( tl.size()==1 ||  tl.size()==2){
        return v;
    }
    for(int i= tl.size()-1; i>0 ;i--){
        if( tl.elemIsConst(i)) {
            //std::cout<< "the const elem are " << i <<'\n';
            res = res->projectOut(i);
            v.push_back(res);
        }
        else{
            res = res->projectOut(i);
        }
    }
    //std::cout << "print "<< s->prettyPrintString() <<'\n';
    return v;
}


DominanceTree* SSA::findPredecessors(DominanceTree* dt) {

    std::vector<std::pair<Set*, int>> stack;
    for (int i = dt->getVectorSize() -1; i >= 0; i--) {
        for (int j = i-1; j >= 0; j--) {
            bool flag = true;
            for(int k: dt->getPredecessor(i)){
                flag = SSA::isReverseDominator( dt->getElem(k),dt->getElem(j) );
            }
            if(flag){
                dt->add_predecessors(i, j);
            }
            // //
            std::vector<Set*>  v = dt->getPrefixes(dt->getElem(j) );
            bool no_match =true;
            std::vector<std::pair<Set*, int>> tmp_stack;

            for(int k=0;k<v.size();k++){
                bool flag  = false;
                int l;
                for( l=0;l<stack.size();l++){
                    if(v[k]==stack[l].first){
                        flag = true;
                        no_match = false;
                        break;
                    }
                }
                if(flag) {
                    int m = 0;
                    while (m < l) {
                        //std::cout << "hello....hello."<<'\n';
                        bool st = SSA::isReverseDominator( dt->getElem(stack[stack.size()-1].second),dt->getElem(j) );
                        if(st){ dt->add_predecessors(stack[stack.size()-1].second, j);}
                       stack.pop_back();
                        m++;
                    }
                    break;
                }else{

                tmp_stack.push_back({v[k], i});}
            }
            if(no_match){stack={};}
            stack.insert(stack.end(), tmp_stack.begin(), tmp_stack.end());
            if(dt->isParent(j,i)) break;
        }
    }
    return dt;
}


bool SSA::isDominator(iegenlib::Set * parent, iegenlib::Set * child){
    // project out constant both parent and child
    Set* childP, *parentP;

    childP = child->projectOutConst(child);
    parentP = parent->projectOutConst(parent);
    

    if( parentP->getArity() > childP->getArity()){
        return  false;
    }
    // project out until parents and child parity matches
    while( parentP->getArity() < childP->getArity()) {
        TupleDecl tl = childP->getTupleDecl();
        childP = childP->projectOut(tl.getSize()-1);
    }
    // child subset of parent;
    return (childP->isSubset(parentP));
}

// create dominator tree
DominanceTree* SSA::createDominanceTree(std::vector<std::pair<int, iegenlib::Set *>>executionS) {
    DominanceTree* pDominanceTree = new DominanceTree();

    //perform the lexicographical sort

    std::sort(executionS.begin(), executionS.end(), []
            (const std::pair<int, iegenlib::Set *> &a, const std::pair<int, iegenlib::Set *> &b) {
        Set* s;
        bool status = s->Set::LexiSort(a.second, b.second);
        return(status);
        delete s;
    });

    // collect nodes into node list
    for (auto v: executionS) {
        pDominanceTree ->push_Back(v);
    }
    // set up the relations( parent and child)
    //
    //
    for (int i = executionS.size() - 1; i >= 0; i--) {
        for (int j = i-1; j >= 0; j--) {
            bool isDominator1 = SSA::isDominator(executionS[j].second, executionS[i].second);
            if(isDominator1){
                pDominanceTree -> add_edge(j, i);
                break;
            }
        }
    }
    return pDominanceTree;
}

void DominanceTree::insertPhiNode(std::vector<std::map<string, std::vector<int>>> globals, Computation* comp ){

    for(int i=0;i<globals.size();i++){
        std:: cout << "hell0 ="<<'\n';
        std::map<string, std::vector<int>> ::iterator it = globals[i].begin();

        std::vector<int> workList = it->second;
        for(int j=0 ;j<workList.size();j++){
            std::cout << "variable's definition node "<< workList[j]<<'\n';
            std::vector<int> DF = this->nodes[workList[j]].dominanceFrontier;
            for(int k=0;k<DF.size();k++){
                //code to insert in phi nodes;

                int df =  DF[k];
                Stmt * s_org = comp->getStmt(nodes[df].data.first);

                //std:: cout<< '---' <<s_org->prettyPrintString()<<'\n';
                std:: cout<< s_org->prettyPrintString()<<'\n';
               // std:: cout << s_org->getIterationSpace()->prettyPrintString()<<'\n';
//
//                Stmt * s2 = new Stmt (
//                        "phi;",
//                        s_org->getIterationSpace(),
//                        s_org->getExecutionSchedule(),
//                        {{it->first, s_org->getIterationSpace()},{it->first, s_org->getIterationSpace()}},
//                        {{it->first,  s_org->getIterationSpace()}}
//                );
//                comp->addStmt(s2);

                std::cout << "insert phi nodes in DF "<< this->nodes[workList[j]].dominanceFrontier[k]<<'\n';
                if (std::find(workList.begin(), workList.end(), DF[k]) == workList.end()) {
                    workList.push_back(DF[k]);
                }
            }
        }
    }
}

void SSA::generateSSA(Computation *  comp) {

    std::vector<Stmt*> stmts  ;
    for( int a=0;a<comp->getNumStmts();a++){
        stmts.push_back(comp->getStmt(a));
    }
    std::vector<std::pair<int, iegenlib::Set*>> executionS;
    std::vector<std::map<string, std::vector<int>>> globals;

    for(int i=0; i<stmts.size(); i++){
            iegenlib::Set* s1 = stmts[i]->getExecutionSchedule()->Apply( stmts[i]->getIterationSpace());
            //std::cout << s1->prettyPrintString()<<'\n';

            std::map<string, std::vector<int>> vec ;
            int numWrites = stmts[i]->getNumWrites();

            for (unsigned int j = 0; j < numWrites; ++j){
                vec[stmts[i]->getWriteDataSpace(j)].push_back(i);
                //std::cout << "variables    " << stmts[i]->getWriteDataSpace(j)<<" ,  "<< i<<'\n';
            }
            globals.push_back(vec);
            executionS.push_back({i, s1});
    }


//    for( auto v: globals){
//        for (auto it = v.cbegin(); it != v.cend(); ++it) {
//            std::cout << "{" << (*it).first << ": "  ;
//            for(int l= 0;l < (*it).second.size();l++) {
//                std::cout << (*it).second[l] <<'\n';
//            }
//
//        }
//    }

    DominanceTree* dt = createDominanceTree(executionS);
    DominanceTree* dt1 = findPredecessors(dt);
    dt1->DFCal();
    dt1->insertPhiNode(globals, comp);
    std::cout << "print random";

}
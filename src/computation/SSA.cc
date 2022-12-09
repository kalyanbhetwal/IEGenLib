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
#include <string>
#include <regex>

using namespace SSA;
using namespace iegenlib;

std::map<Stmt*, std::vector<Stmt*>> SSA::Member::predecessor{};
std::map<Stmt*, std::vector<Stmt*>> SSA::Node::DF{};
std::map<string, std::vector<Stmt*>> SSA::Node::globals{};
std::map<Stmt*, std::vector<Stmt*>> SSA::Member::possiblePaths;

//std::map<Stmt*, std::vector<Stmt*>> SSA::predecessor;


std::vector<Set*> SSA::getPrefixes(Set*s) {

    std::vector<Set*>v;
    Set* res = new Set(*s);

    TupleDecl tl = s->getTupleDecl();
    v.push_back(res);

    if( tl.size()==1 ){
        return(v);
    }

    for(int i= tl.size()-1; i>0 ;i--) {
        res = res->projectOut(i);
        v.push_back(res);
        //std:: cout <<" the extracted set " << res->prettyPrintString()<<std::endl;
    }
    return v;
}

void SSA::Node::printDF() {
    for (auto m: SSA::Node::DF) {
        std::cout << "DF for node " << m.first->getExecutionSchedule()->prettyPrintString() << std::endl;
        for (int i = 0; i < m.second.size(); i++) {
            std::cout << "  is " << m.second[i]->getExecutionSchedule()->prettyPrintString() << std::endl;
        }
        std::cout << "-------===------------" << std::endl;
    }
}
void SSA::Node::printPredDom(){

    for(auto m: SSA::Member::predecessor){
        std::cout<< "the pred dom list for node  " << m.first->getExecutionSchedule()->prettyPrintString() <<std::endl;
        for (int i = 0; i < m.second.size(); i++) {
            std::cout << "  is " << m.second[i]->getExecutionSchedule()->prettyPrintString() << std::endl;
        }
        std::cout << "-------===------------"<<std::endl;
    }

}


 SSA::Node* SSA::createScheduleTree(iegenlib:: Computation* comp){

    std::vector<Stmt*> stmts  ;

    Node * rootNode = new Node();
    rootNode->setOrdered( true);
    rootNode->setCommonArity(1);
    rootNode->setParent(NULL,NULL);

    // remove this for loop
    for( int a=0;a<comp->getNumStmts();a++){
        stmts.push_back(comp->getStmt(a));
    }

    for(int i=0; i<stmts.size(); i++){
        //std:: cout << "the stmt is " << stmts[i]->prettyPrintString() <<std::endl;
        iegenlib::Set* s1 = stmts[i]->getExecutionSchedule()->Apply( stmts[i]->getIterationSpace());

        string s  = s1->getString();
        s.erase(remove(s.begin(), s.end(), '$'), s.end());

        Set * s2 = new Set(s);

        //std:: cout << "execution schedule "<< s2->prettyPrintString() << std::endl;
        int numWrites = stmts[i]->getNumWrites();
        for (int j = 0; j < numWrites; ++j){
            Node::globals[stmts[i]->getWriteDataSpace(j)].push_back(stmts[i]);
        }

        //std::cout << s1->prettyPrintString()<<'\n';
        std::vector<Set*>v;
        v = getPrefixes(s2);
        SSA::Node * current = rootNode;

        for(int j= v.size()-1;j>=0;j--){
            //std:: cout << "prefixes " << (*v[j]).prettyPrintString()<<'\n';
            SSA::Member * m;
            if ( j ==0){
                m = new SSA::Member(v[j], stmts[i]);
            }
            else{
                m = new SSA::Member(v[j], NULL);
            }
            current =  current->insert(m);
        }
    }
    return rootNode;
}

SSA::Node::Node(){
     members = {};
}

SSA::Member::Member(Set * s, Stmt * st) {
    schedule = s;
    child = new SSA::Node();
    stmt = st;
}

Set *Member::getSchedule()  {
    return schedule;
}

void Member::setSchedule(Set *schedule) {
    Member::schedule = schedule;
}

Stmt *Member::getStmt()  {
    if(stmt==NULL){
        return NULL;
    }
    return stmt;
}

void Member::setStmt(Stmt *stmt) {
    Member::stmt = stmt;
}

Node *Member::getChild()  {
    return child;
}

void SSA::Member::setChild(SSA::Node *child) {
    child = child;
}

 std::pair<SSA::Node *, SSA::Member*> &SSA::Node::getParent() {
    return parent;
}

void Node::setParent(Node * n, Member* s) {
    parent = std::make_pair(n,s);
}
void SSA::rename(Computation * comp){
    int counter = 0;
    for (int a = 0; a < comp->getNumStmts(); a++) {
        Stmt *s;
        s = comp->getStmt(a);

        for (int j = 0; j < s->getNumWrites(); j++) {
            std::string write = s->getWriteDataSpace(j);

            write.erase(write.begin());
            write.erase(write.end()-1);
            string newWrite = write + "_"+ std::to_string(counter);

            s->replaceWrite(write,  newWrite);
        }
        counter++;
    }
}

void SSA::renameSSA(Computation* comp){
    Computation* phiComp = new Computation();
    std::map<string, std::vector<Stmt*>>::iterator it;
    std::map<string,  std::map<Stmt *, std::vector<Stmt *>>> readLoc;
    std::map<Stmt* , Stmt*> phi_to_stmt;
    std::map<Stmt* , Stmt*> stmt_to_phi;
    for (it = SSA::Node::globals.begin(); it != SSA::Node::globals.end(); it++) {
        string newName = it->first;
        newName.erase(newName.begin());
        newName.erase(newName.end() - 1);
        std::map<Stmt *, std::vector<Stmt *>> phiLoc;
        for (int v = 0; v < it->second.size(); v++) {
            if (Node::DF.find(it->second[v]) != Node::DF.end()) {
                std::vector<Stmt *> insert_phi_at = Node::DF.at(it->second[v]);
                for (auto stmt: insert_phi_at) {

                    if (std::find(phiLoc[stmt].begin(), phiLoc[stmt].end(), it->second[v]) == phiLoc[stmt].end()) {
                        phiLoc[stmt].push_back(it->second[v]);
                    }

                }
            }
        }

        // if there is direct connection between two nodes or statements
        // check if the immediate node of a phi's read
        for(auto  missingRead=phiLoc.begin();missingRead!=phiLoc.end();missingRead++){
            Stmt* mp = SSA::Member::predecessor[missingRead->first].back();

            if(std::find(missingRead->second.begin(), missingRead->second.end(),mp ) == missingRead->second.end()){
                //std::cout << "adding missing read"<<std::endl;
                phiLoc[missingRead->first].push_back(mp);
            }

        }

        std::map<Stmt *, std::vector<Stmt *>>::iterator phis;
        for (phis = phiLoc.begin(); phis != phiLoc.end(); phis++) {

            int count = 0;
            for(auto i:phis->second){
                for (int j = 0; j < i->getNumWrites(); j++){
                    if(i->getWriteDataSpace(j)==it->first){
                        count= count + 1;
                        break;
                    }
                }
            }
            if(count<2)continue;

            Stmt *phi = new Stmt(
                    "phi",
                    phis->first->getIterationSpace()->getString(),
                    phis->first->getExecutionSchedule()->getString(),
                    {{newName, "{[0]->[0]}"}},
                    {{newName, "{[0]->[0]}"}}
            );
            phi->setPhiNode(true);
            phiComp->addStmt(phi);
            //TODO:replace with multimap
            phi_to_stmt[phi] = phis->first;
            stmt_to_phi[phis->first] = phi;

        }
        readLoc[it->first] = phiLoc;

    }

    //insert definition phis at certain locations
    std::map<Stmt *, Stmt *> merge_to_stmt;
    std::map<Stmt *, Stmt *> stmt_to_merge;
    std::map<string,  std::map<Stmt *, Stmt *>> globalsMap;
    for(int k=0; k<comp->getNumStmts(); k++){
        Stmt * st = comp->getStmt(k);
        if(!st->isPhiNode()){
            for (int j = 0; j < st->getNumWrites(); j++) {
                string newName = st->getWriteDataSpace(j);
                newName.erase(newName.begin());
                newName.erase(newName.end() - 1);
                //std::cout << "the num write space " << st->getWriteDataSpace(j) << std::endl;

               string es = st->getExecutionSchedule()->getString();
               string new_es;
               std::regex rgx("(.*),(.*)](.*)", std::regex::extended);
               std::smatch matches;

                if (std::regex_search(es, matches, rgx)) {
                    new_es = matches[1].str()+',';
                    new_es = new_es + std::to_string((stoi(matches[2])+1)) ;
                    new_es = new_es + ']'+ matches[3].str();
                }
                if(new_es.empty()){
                    new_es = es;
                }
                //std::cout << " the updated execution schedule  " << new_es << std::endl;
                Stmt *phi = new Stmt(
                        "phi",
                        st->getIterationSpace()->getString(),
                        new_es,
                        {{newName, "{[0]->[0]}"}},
                        {{newName, "{[0]->[0]}"}}
                );

                phi->setPhiNode(true);
                phi->setDefPhi(true);
                phiComp->addStmt(phi);
                merge_to_stmt[phi] = st;
                stmt_to_merge[st] = phi;
                globalsMap[st->getWriteDataSpace(j)].insert(std::make_pair(st,phi ));
            }

        }
    }

    //SSA::Member::predecessor={};
    Node* n = createScheduleTree(phiComp);
    n->calc_all_pred();
    //n->printPredDom();
    n->calc_all_backward_paths();

    for (int b = 0; b < phiComp->getNumStmts(); b++) {
        //std::cout<<"the num of statement" <<b<<std::endl;
        Stmt *s1;
        s1 = phiComp->getStmt(b);
        comp->addStmt(s1);
    }

    // rename all the phi nodes variables and statements
    rename(comp);

    for (int b = 0; b < comp->getNumStmts(); b++) {
       // std::cout<<"#no stmt " <<b <<std::endl;
        Stmt *s1;
        s1 = comp->getStmt(b);

        for (int j = 0; j < s1->getNumReads(); j++) {
            std::string read = s1->getReadDataSpace(j);
            std::string test = read;
            test.erase(test.end()-1);

            test = test + '_';

            if(s1->isDefPhi()){
                std::vector<Stmt*> s;
                s = SSA::Member::possiblePaths[s1];
                //std::cout<< "var "<<s1->getWriteDataSpace(0)<<std::endl;
                //std:: cout << "the stmt "<< s1->getExecutionSchedule()->prettyPrintString() <<std::endl;
                s.push_back(merge_to_stmt[s1]);
                bool flag = false;
                for(auto x: s){
                    //std::cout << "reads " << x->getExecutionSchedule()->prettyPrintString() <<std::endl;
                    //if(!x->isPhiNode()) continue;
                    for (int j = 0; j < x->getNumWrites(); j++){
                        if(x->getWriteDataSpace(j).find(test)!= std::string::npos){
                            flag = true;
                            s1->addRead(x->getWriteDataSpace(j),"{[0]->[0]}");
                            break;
                        }
                    }
                }
               //std::cout <<"-----------------------"<<std::endl;
                if(flag){s1->removeReadDataSpace(0);
                break;}
                //std::cout <<"-----------------------"<<std::endl;
            }
            if( s1->isPhiNode() && !s1->isDefPhi()){
                std::vector<Stmt*> slist;
                Stmt * actualS = phi_to_stmt[s1];
                slist = SSA::Member::predecessor[actualS];
                //std::cout <<"-----------------------"<<std::endl;
                //std:: cout << "the stmt "<< s1->getExecutionSchedule()->prettyPrintString() <<std::endl;
                for(auto v: slist){
                    if(globalsMap.find(read)!= globalsMap.end()){
                        v = globalsMap[read][v];
                       // std::cout <<"match"<<std::endl;
                    }
                    if(!v)continue; ///TODO:: check why is this happening
                    if(!v->isPhiNode()) continue;
                   // std::cout << "reads " << v->getExecutionSchedule()->prettyPrintString() <<std::endl;

                    for (int j = 0; j < v->getNumWrites(); j++){
                        if(v->getWriteDataSpace(j).find(test)!= std::string::npos){
                            s1->addRead(v->getWriteDataSpace(j),"{[0]->[0]}");
                            break;
                        }
                    }
                }
                s1->removeReadDataSpace(0);
                break;
            }
            // trying to edit reads on complementary node of phi node
            else if(readLoc[read].find(s1)!=readLoc[read].end()     && stmt_to_phi.find(s1)!=stmt_to_phi.end() ){
                Stmt * s_phi = stmt_to_phi[s1];
                //std::cout << s_phi->getExecutionSchedule()->prettyPrintString()<<std::endl;
                if(!s_phi)continue; //TODO::check why is this happening
                for (int l = 0; l < s1->getNumReads(); l++){
                    if(s1->getReadDataSpace(l)==read){
                       // std::cout <<"match"<<std::endl;
                        s1->replaceReadDataSpace( s1->getReadDataSpace(l), s_phi->getWriteDataSpace(0));
                        //std::cout<< "the dataspace " << s1->getReadDataSpace(l)<<std::endl;
                        break;
                    }
                }
                //s1->removeReadDataSpace(0);
            }
            else{
                if(s1->isPhiNode())continue;
                std::vector<Stmt*> pred= SSA::Member::predecessor[s1];
                //std::cout << "the pred size is "<< pred.size()<<std::endl;
                Stmt *s_pred;
                if(pred.size()==0)continue;
                if(pred.size()>1) {
                    /// if multiple predecessors are writing to same data space then we should throw and exception
                //    std::cout << "read "<< read<< std::endl;
                   // std::cout << "stmt " << s1->getExecutionSchedule()->prettyPrintString() << std::endl;

                      //  std::cout <<"pred "<< pred[pred.size()-1]->getExecutionSchedule()->prettyPrintString() << std::endl;

                        if (globalsMap.find(read) != globalsMap.end()) {
                            if(globalsMap[read].size()>1){
                              //  std::cout << "------multimap -----"<<std::endl;
                                s_pred = globalsMap[read][pred[pred.size()-1]];
                            }else s_pred = globalsMap[read].begin()->second;


                            //s_pred = stmt_to_merge[sorg];
                            //std::cout << "matched  "<< sorg->getExecutionSchedule()->prettyPrintString() << std::endl;
//                            if (globalsMap[read].find(pred[d]) != globalsMap[read].end()) {
//                                s_pred = globalsMap[read][pred[d]];
//                                //std::cout << "matched  "<< pred[d]->getExecutionSchedule()->prettyPrintString() << std::endl;
//                                //std::cout<<"__________________________________"<<std::endl;
//                                break;
//                            }
                        }

                }else {
                    if (globalsMap.find(read) != globalsMap.end()) {
                        if(globalsMap[read].size()>1){
                            //std::cout << "multimap "<<std::endl;
                            s_pred = globalsMap[read][pred[0]];
                        }else s_pred = globalsMap[read].begin()->second;
                    }
                    //s_pred =stmt_to_merge[pred[0]];

                }
                    if(!s_pred)continue; //TODO:: check why is this happening
                    int k;
                    bool match = false;
                    for ( k = 0; k < s_pred->getNumWrites(); k++){
                        if(s_pred->getWriteDataSpace(k).find(test)!= std::string::npos){
                            match = true;
                            break;
                        }
                    }
                    if( s_pred->getNumWrites()>0 and match) {
                        for (int l = 0; l < s1->getNumReads(); l++) {
                           // std:: cout << s1->getReadDataSpace(l) << "  tttyyyttt   "<< test <<std::endl;
                            if (s1->getReadDataSpace(l) == read) {
                                //std:: cout << s1->getReadDataSpace(l) << "  t  "<< s_pred->getWriteDataSpace(k) << "  the test is " << test<<std::endl;
                                s1->replaceReadDataSpace(s1->getReadDataSpace(l), s_pred->getWriteDataSpace(k));
                                break;
                            }
                        }
                    }

                }


            }
        }

    }

void SSA::generateSSA(iegenlib::Computation *comp) {
    Node * node = createScheduleTree(comp);
    node->calc_all_pred();

    node->printPredDom();

    node-> computeDF();

    //node->printDF();

    SSA::renameSSA(comp);



}

void SSA::Node::computeDF() {
    std::map<Stmt*, std::vector<Stmt*>>::iterator it;
    //for all statements in computation
    for (it = Member::predecessor.begin(); it != Member::predecessor.end(); it++)
    {
        Stmt* runner;
        //for all pred of that statement
       // std::cout << "it "<< it->first->getExecutionSchedule()->prettyPrintString()<<std::endl;
        if(it->second.size()> 1) {
            for (int j = 0; j < it->second.size(); j++) {
                runner =(Stmt*)it->second[j];
                while (runner->getExecutionSchedule()->toString() != it->second[it->second.size() -1]->getExecutionSchedule()->toString()) {
                    //std:: cout << "p   "<< it->first->getExecutionSchedule()->prettyPrintString()<<std::endl;
                   // std:: cout << "r  "<< runner->getExecutionSchedule()->prettyPrintString()<<std::endl;
                   if(runner!=it->first) {
                       Node::DF[runner].push_back(it->first);
                   }
                    runner = Member::predecessor[runner][Member::predecessor[runner].size()-1];
                }
                //std::cout <<"-----------------------------------------------"<<std::endl;
            }

        }
    }
}

SSA::Node* SSA::Node::insert(SSA::Member * m){
    if(m->getSchedule()->getArity() != common_arity){
        return NULL;
    }
    if(ordered){
        for(auto current=members.begin(); current!=members.end(); ++current ){
            //current = members[i];
            if((*(*current)->getSchedule())== (*m->getSchedule())){
                return (*current)->getChild();
            }

            if(!((*current)->getSchedule())->LexiLess(m->getSchedule())){
                members.emplace(current, m);
                m->getChild()->setParent(this, m);
                m->getChild()->setCommonArity(getCommonArity()+1);
                return m->getChild();
            }

        }
        members.push_back(m);
        m->getChild()->setCommonArity(getCommonArity()+1);
        m->getChild()->setParent(this,m );
        return m->getChild();
    }
    for(auto current=members.begin(); current!=members.end();current++ ){
        if((*(*current)->getSchedule())== (*m->getSchedule())){
            return (*current)->getChild();
        }
    }
    members.push_back(m);
    m->getChild()->setCommonArity(getCommonArity()+1);
    m->getChild()->setParent(this,m ); // parent == null
    return m->getChild();
}

bool SSA::Node::isOrdered()  {
    return ordered;
}

void SSA::Node::setOrdered(bool ordered) {
    Node::ordered = ordered;
}

int SSA::Node::getCommonArity()  {
    return common_arity;
}

void SSA::Node::setCommonArity(int commonArity) {
    common_arity = commonArity;
}

void SSA::Node::printBreadthFirst() {
    for(auto it=members.begin(); it!=members.end();it++){
        (*it)->printBreadthFirst();
    }
        std::cout << "------------------"<<'\n';
}


void SSA::Node::calc_all_pred() {

    for(auto it=members.begin(); it!=members.end();it++){
        (*it)->calc_all_pred(this);
    }
    //std::cout << "------------------"<<'\n';


}

void SSA::Node::calc_all_backward_paths() {

    for(auto it=members.begin(); it!=members.end();it++){
        (*it)->calc_all_backward_path(this);
    }
    //std::cout << "------------------"<<'\n';


}

void Node::setMembers( std::vector<Member *> &members) {
    Node::members = members;
}


void SSA::Member::printBreadthFirst() {
    std::cout << schedule->prettyPrintString()<<'\n';
    child->printBreadthFirst();
}

void SSA::Member::calc_all_pred(Node * n){

    if(stmt!=NULL) {
        int j;
        for (j = 0; j < n->getMembers().size(); j++) {
            if (this == n->getMembers()[j]) {
                break;
            }
        }
        std::vector<Stmt *> stmtList;

        stmtList = pred_and_dom(n, j - 1);

        std::vector<Stmt*> rduplicates{};

        //rduplicates  = predecessor[stmt];

        for (int i = 0; i < stmtList.size(); i++) {
            if(stmtList[i]== stmt){
                continue;
            }
            if(std::find(rduplicates.begin(), rduplicates.end(),stmtList[i] ) == rduplicates.end()){
                rduplicates.push_back(stmtList[i]);
            }
        }
        predecessor[stmt] = rduplicates;

    }
    child->calc_all_pred();
}

void SSA::Member::calc_all_backward_path(Node * n){

    if(stmt!=NULL) {
        int j;
        for (j = 0; j < n->getMembers().size(); j++) {
            if (this == n->getMembers()[j]) {
                break;
            }
        }
        std::vector<Stmt *> stmtList;

        stmtList = pred_path(n, j - 1);

        std::vector<Stmt*> rduplicates;

        rduplicates  = possiblePaths[stmt];

        for (int i = 0; i < stmtList.size(); i++) {
            if(stmtList[i]== stmt){
                continue;
            }
            if(std::find(rduplicates.begin(), rduplicates.end(),stmtList[i] ) == rduplicates.end()){
                rduplicates.push_back(stmtList[i]);
            }
        }
        possiblePaths[stmt] = rduplicates;

    }
    child->calc_all_backward_paths();
}

std::vector<Stmt*> SSA::Member::pred_path(Node* n, int idx) {
    std::vector < Stmt * > listOfStatements{};
    int i;
    for (i = idx; i >= 0; i--) {

        //this case is for when we hit a dominator
        if (n->getMembers()[i]->getStmt() != NULL) {
            listOfStatements.push_back(n->getMembers()[i]->getStmt());
	    //return listOfStatements;
        }
        //this case is for when we are adding predecessors that aren't dominators
        for (auto c: n->getMembers()[i]->getChild()->getMembers()) {
            std::vector < Stmt * > s;
            s = pred_and_dom(c->getChild(), c->getChild()->getMembers().size() - 1);
            listOfStatements.insert(listOfStatements.end(), s.begin(), s.end());
        }
    }
    // if we are here we did not find a dominator within the loop.
    // Now we have to do three things.
    // 1. we have to find possible predecessors by looking at parent unordered node
    // 2. we have to find the dominator by looking at grandparent node
    // 3. we have to find possible intraloop predecessors by looking from
    // the end of the loop backward until we reach idx
    // I am going to do step 3 first, then 1, and then 2
    //
    // This is step 3 above.
    if(n->getCommonArity()!=1 && n->getMembers().size()>1) {
        for (i = n->getMembers().size() - 1; i != idx; i--) {
            //this case is for when we hit a dominator
            if (n->getMembers()[i]->getStmt() != NULL) {
                listOfStatements.push_back(n->getMembers()[i]->getStmt());
                break;
            }
            //this case is for when we are adding predecessors that aren't dominators
            for (auto c: n->getMembers()[i]->getChild()->getMembers()) {
                std::vector<Stmt *> s;
                s = pred_and_dom(c->getChild(), c->getChild()->getMembers().size() - 1);
                listOfStatements.insert(listOfStatements.end(), s.begin(), s.end());
            }
        }
    }
    // this is for the root node
   // std::cout <<" root node "<< n->getParent().first << std::endl;

    if (n->getParent().first == NULL) {
        return listOfStatements;
    }
    // stepping up to find the location of the dominator in the member vector
    Node* p = n->getParent().first;
    for (auto c:p->getMembers()){
        if(c->getChild()!= n) {
            std::vector<Stmt *> s;
            s = pred_and_dom(c->getChild(), c->getChild()->getMembers().size() - 1);
            listOfStatements.insert(listOfStatements.end(), s.begin(), s.end());
        }
    }
   // std::cout <<" get parent "<< p->getParent().second->getSchedule()->prettyPrintString() << std::endl;

    Node * gp = p->getParent().first;
    Member * gpm = p->getParent().second;
    if(gp != NULL){
        std::vector<Stmt*> s;
        int j;
        for(j=0;j<gp->getMembers().size();j++ ){
            if(gpm==gp->getMembers()[j] ){
                break;
            }
        }
        s = pred_and_dom(gp,j-1);
        listOfStatements.insert(listOfStatements.end(), s.begin(), s.end());
    }
    return listOfStatements;
}


std::vector<Stmt*> SSA::Member::pred_and_dom(Node* n, int idx) {

   // std:: cout << "The ordered value is "<< n->isOrdered()<<std::endl;
    std::vector < Stmt * > listOfStatements{};
    int i;
    for (i = idx; i >= 0; i--) {

        //this case is for when we hit a dominator
        if (n->getMembers()[i]->getStmt() != NULL) {
            listOfStatements.push_back(n->getMembers()[i]->getStmt());
            return listOfStatements;
        }
        //this case is for when we are adding predecessors that aren't dominators
        for (auto c: n->getMembers()[i]->getChild()->getMembers()) {
            std::vector < Stmt * > s;
            s = pred_and_dom(c->getChild(), c->getChild()->getMembers().size() - 1);
            listOfStatements.insert(listOfStatements.end(), s.begin(), s.end());
        }
    }
    // if we are here we did not find a dominator within the loop.
    // Now we have to do three things.
    // 1. we have to find possible predecessors by looking at parent unordered node
    // 2. we have to find the dominator by looking at grandparent node
    // 3. we have to find possible intraloop predecessors by looking from
    // the end of the loop backward until we reach idx
    // I am going to do step 3 first, then 1, and then 2
    //
    // This is step 3 above.
    if(n->getCommonArity()!=1 && n->getMembers().size()>1) {
        for (i = n->getMembers().size() - 1; i != idx; i--) {
            //this case is for when we hit a dominator
            if (n->getMembers()[i]->getStmt() != NULL) {
                listOfStatements.push_back(n->getMembers()[i]->getStmt());
                break;
            }
            //this case is for when we are adding predecessors that aren't dominators
            for (auto c: n->getMembers()[i]->getChild()->getMembers()) {
                std::vector<Stmt *> s;
                s = pred_and_dom(c->getChild(), c->getChild()->getMembers().size() - 1);
                listOfStatements.insert(listOfStatements.end(), s.begin(), s.end());
            }
        }
    }
    // this is for the root node
   // std::cout <<" root node "<< n->getParent().first << std::endl;

    if (n->getParent().first == NULL) {
        return listOfStatements;
    }
    // stepping up to find the location of the dominator in the member vector
    Node* p = n->getParent().first;
    for (auto c:p->getMembers()){
        if(c->getChild()!= n) {
            std::vector<Stmt *> s;
            s = pred_and_dom(c->getChild(), c->getChild()->getMembers().size() - 1);
            listOfStatements.insert(listOfStatements.end(), s.begin(), s.end());
        }
    }
   // std::cout <<" get parent "<< p->getParent().second->getSchedule()->prettyPrintString() << std::endl;

    Node * gp = p->getParent().first;
    Member * gpm = p->getParent().second;
    if(gp != NULL){
        std::vector<Stmt*> s;
        int j;
        for(j=0;j<gp->getMembers().size();j++ ){
            if(gpm==gp->getMembers()[j] ){
                break;
            }
        }
        s = pred_and_dom(gp,j-1);
        listOfStatements.insert(listOfStatements.end(), s.begin(), s.end());
    }
    return listOfStatements;
}

std::vector<Member*> SSA::Node::getMembers(){
    if(members.empty())  return std::vector<Member*>();
    return members;
}

/*!
 * \file simplifyForPartialParallelDriver.cc
 *
 * This file is a driver for using simplifyForPartialParallel 
 * and pertaining functions in IEgenLib.
 * These functions can be used to carry out simplification of 
 * data dependence relation in non-affine relations containing
 * uninterpreted function calls.
 *
 * To compile, after building IEgenLib with:

     ./configure    
     make

   Run the following command with YOUR OWN ADDRESSES:

     g++ -o EXECUTABLE simplifyDriver.cc -I IEGENLIB_HOME/src \
         IEGENLIB_HOME/build/src/libiegenlib.a -lisl -std=c++11

   You can also run following command after running "make install" for IEgenLIB:

       g++ -o EXECUTABLE simplifyDriver.cc \
           -I INSTALLATION_FOLDER/include/iegenlib cpp_api_example.cc \
              INSTALLATION_FOLDER/lib/libiegenlib.a -lisl -std=c++11
   
   IEGENLIB_HOME indicates where you have your copy of IEgenLIB.
   For example if you are compiling this file in its original location that is 
   IEGENLIB_HOME/src/drivers, FROM IEGENLIB_HOME run following to compile:

     g++ -o src/drivers/simplifyDriver_DL src/drivers/simplifyDriver_DL.cc -I src build/src/libiegenlib.a -lisl -std=c++11

 * Now to run the driver, you should put your dependence relations inside
 * JSON files and give them as inputs to the driver, one or more files
 * at a time. JSON file format is demonstrated by example sample.json 
 * in the same directory as this driver as well as several *.json files 
 * in IEGENLIB_HOME/data directory. The json input files include comment 
 * fields that have the application code, and says where in the code each
 * dependence relation is getting extracted from. So you can run following
 * after compiling the driver to get the simplified relations:
 
   build/bin/simplifyDriver sample.json

 * Note: the driver can read one or as many input json files as given.

 *     
 *
 * \date Date Started: 4/16/2017
 *
 * \authors Michelle Strout, Mahdi Soltan Mohammadi
 *
 * Copyright (c) 2016, University of Arizona <br>
 * All rights reserved. <br>
 * See COPYING for details. <br>
 * 
 *
 * HOW TO USE INTERFACE TO SIMPLIFICATION ALGORITHM for
 * simplifying constraints relation containing Uninterpreted Function Calls.
 * the overall process is explained in the 7 steps bellow, the simplify function
 * implementes these steps using IEGenLib's Functions:
 * 
 * (1) You need to define uninterpreted function calls (UFCs) that appear
 *        in constraints for iegenlib environment. These informations are 
 *        the first set of information that simplify function reads from 
 *        input json files and puts them in the IEGenLib enviroment.
 *        Note that, you need to do this only once for 
 *        relations that have same UFSs during one run.
 *
 * (2) You need to put constraints in iegenlib Relation (or Set).
 * 
 * (3) Determing unsatisfiability
 * 
 * (4) Create a std::set that includes which tuple variables (loop iterators)
 *     we are not going to project out. 
 *
 * (5) Apply a heuristic to remove expensive constraints that
 *     is keeping us from projecting out an iterator:
         for instance: 
                        col(j) < n  would keep us from projecting 'j'
       We only remove constraints up to a maximum number determined by user.
 * 
 * (6) Use simplifyForPartialParallel function that is main interface for the
 *     algorithm. If relation is not satisfiable the function would return NULL.
 *     Otherwise it would return the simplified result as iegenlib Relation.
 * 
 * (7) Print out result (if not NULL) using toISLString() function.  
 
 *  We have demonstrated these steps in the simplify function
 *  This function reads information from a JSON file (inputFile), and applies
 *  the simplification algorithm to the sets found in the file. 
 */


#include <iostream>
#include "iegenlib.h"
#include "parser/jsoncons/json.hpp"

using jsoncons::json;
using namespace iegenlib;
using namespace std;


// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//! If you wish not to see the output, change the value: verbose = false;
bool verbose=true;

void simplify(string inputFile);

// Utility function
void printRelation(string msg, Relation *rel);
void EXPECT_EQ(string a, string b);
void EXPECT_EQ(Relation *a, Relation *b);
int str2int(string str);


//----------------------- MAIN ---------------
int main(int argc, char **argv)
{

  if (argc == 1){
    cout<<"\n\nYou need to specify the input JSON files (one or more) "
          "that contain dependence relations:"
          "\n./simplifyDriver file1.json file2.json\n\n";
  } else if (argc >= 2){
    // Parsing command line arguments and reading given files.
    for(int arg = 1; arg < argc ; arg++){
      simplify(string(argv[arg]));
    }
  }

  return 0;
}


// Reads information from a JSON file (inputFile), and applies
// the simplification algorithm to the sets found in the file. 
void simplify(string inputFile)
{

  iegenlib::setCurrEnv();
  std::set<int> parallelTvs;
  // (0) Read the data from inputFile
  ifstream in(inputFile);
  json data;
  in >> data;

  string dataLogOutputfile = inputFile.erase(inputFile.length()-5,5);
  string resultFile = dataLogOutputfile;
  dataLogOutputfile.append(".dl");
  resultFile.append(".result");;
  char buffer[200];
  sprintf(buffer, "rm -f %s", resultFile.c_str() );
  system(buffer);

  int p = 0;
  //cout<<data[p][0]["Name"].as<string>()<<"\n\n";
  for (size_t i = 0; i < data[p].size(); ++i){// Conjunctions found for one DR in the file

    sprintf(buffer, "cp data/genRules.dl %s", dataLogOutputfile.c_str() );
    system(buffer);

    // (1) Putting constraints in an iegenlib::Relation
    // Reading original set.
    Relation* rel = new Relation(data[p][i]["Relation"].as<string>());

    // (2) Introduce the uninterpreted function symbols to environment, and 
    // indicate their domain, range, whether they are bijective, or monotonic.
    if( i == 0 ){  // Read these data only once. 
                   // They are stored in the first conjunction.

      json ufcs = data[p][i];
      // Read UFCs' data for code No. p, from ith relation
      addUFCs(ufcs);

      // Add defined domain information to environment
      json uqCons = data[p][0]["User Defined"];
      adduniQuantConstraints(uqCons);

      // Specify loops that are going to be parallelized,
      json npJ = data[p][i]["Do Not Project Out"];
      notProjectIters( rel, parallelTvs, npJ);
    }

    Relation* relWithDR = rel->boundDomainRange();

//if(i==0) printRelation(string("Orig? = "), rel);

    //Relation *datalog_rel = rel->checkUnsat(dataLogOutputfile);
    Relation *datalog_rel = relWithDR->checkUnsat(dataLogOutputfile);

    ofstream rOut(resultFile.c_str(), std::ofstream::out | std::ofstream::app);
    rOut<<"\n\n####### Relation "<<i<<":"<<rel->toISLString()<<"\n\n";
    rOut.close();

    // Reading datalog rules that are drectly put in the Json files
    std::ofstream dlOut(dataLogOutputfile, std::ofstream::out | std::ofstream::app);//
    json uqCons = data[p][0]["User Defined"];
    for (size_t j = 0; j < uqCons.size(); ++j){

      // Read user defined Datalog rulesconstraints based on universally quantified
      if( uqCons[j]["Type"].as<string>() == "UserDefDataLogRule"){

        dlOut<<"\n"<<uqCons[j]["Rule:"].as<string>();
      }
    }
    dlOut.close();

    sprintf(buffer, "souffle %s >> out", dataLogOutputfile.c_str());//resultFile.c_str());
    system(buffer);
//    printRelation(string("Orig = "), rel);
//    printRelation(string("Purf = "), datalog_rel);
  }

// } // End of p loop
}








void printRelation(string msg, Relation *rel){

    if ( rel ) {
        cout<<"\n\n"<<msg<<rel->toISLString()<<"\n\n";
    } else {
        cout<<"\n"<<msg<<"Not Satisfiable"<<"\n";
    }
}

void EXPECT_EQ(string a, string b){

    if( a != b ){
        cout<<"\n\nExpected: "<<a;
        cout<<"\n\nActual:"<< b <<"\n\n";
    }
}

void EXPECT_EQ(Relation *a, Relation *b){

    if( a != b ){
        cout<<"\n\nIncorrect results: Expected or Actual is NULL.\n\n";
    }
}

int str2int(string str){
  int i;
  sscanf (str.c_str(),"%d",&i);
  return i;
}

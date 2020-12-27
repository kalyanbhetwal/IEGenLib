/*!
 * \file Computation.h
 *
 * \brief Declarations for the Computation and Stmt classes.
 *
 * The Computation class is the SPF representation of a logical computation.
 * It contains a Stmt class for each statement, which in turn contains
 * information about that statement as mathematical objects.
 * Originally part of spf-ie.
 *
 * \date Started: 10/09/20
 *
 * \authors Anna Rift
 *
 * Copyright (c) 2020, University of Arizona <br>
 * Copyright (c) 2020, Boise State University <br>
 * All rights reserved. <br>
 * See ../../COPYING for details. <br>
 */

#ifndef COMPUTATION_H_
#define COMPUTATION_H_

#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "set_relation/Visitor.h"
#include "set_relation/environment.h"
#include "set_relation/set_relation.h"

namespace iegenlib {

class Stmt;

/*!
 * \class VisitorChangeUFsForOmega
 *
 * \brief Visitor that makes modifications to UF call terms so that they are
 * acceptable for what Omega supports.
 *
 * Intended for use on sets/relations. Also gathers information needed to pass
 * along to Omega codegen, such as #define macros.
 */
class VisitorChangeUFsForOmega : public Visitor {
   private:
    //! string stream for building up necessary UF call macros
    std::ostringstream macros;
    //! equality constraints that must be added to the current conjunction to
    //! make nested UF substitutions valid
    std::vector<Exp*> UFSubstitutionConstraints;
    //! symbolic constants to use in place of UF calls that become 0-args
    std::map<UFCallTerm*, VarTerm*> zeroArgUFReplacements;
    //! next number to use in creating unique function names
    int nextFuncReplacementNumber = 0;
    //! next number to use in creating replacement variable names
    int nextVarReplacementNumber = 0;

   public:
    //! Construct a new VisitorChangeUFsForOmega
    VisitorChangeUFsForOmega() {}

    //! Get the UF call macros required for the code corresponding to the
    //! set/relation to function correctly, as a string
    std::string getMacrosString() { return macros.str(); }

    void postVisitConjunction(Conjunction* conj) {
        // add constraints on replacement variables to make UF subs valid
        for (const auto& constraint : UFSubstitutionConstraints) {
            conj->addEquality(constraint);
        }
        UFSubstitutionConstraints.clear();
    }

    void postVisitExp(Exp* exp) {
        // replace 0-arg UF calls we found while traversing with symbolic
        // constants
        for (const auto& originalTerm : exp->getTermList()) {
            UFCallTerm* originalTermAsUF =
                dynamic_cast<UFCallTerm*>(originalTerm);
            auto it = zeroArgUFReplacements.begin();
            while (it != zeroArgUFReplacements.end()) {
                // match term with one that must be replaced
                if (it->first == originalTermAsUF) {
                    // perform replacement
                    // subtract original UF term
                    UFCallTerm* subtractionTerm =
                        static_cast<UFCallTerm*>(originalTermAsUF->clone());
                    subtractionTerm->setCoefficient(
                        -1 * subtractionTerm->coefficient());
                    exp->addTerm(subtractionTerm);
                    // add new symbolic constant term
                    exp->addTerm(it->second->clone());
                    // remove replacement rule after it has been applied, to
                    // avoid extra searching later
                    it = zeroArgUFReplacements.erase(it);
                } else {
                    it++;
                }
            }
        }
    }

    void postVisitUFCallTerm(UFCallTerm* callTerm) {
        // set up macro outputs
        std::ostringstream os_replaceFrom;
        std::ostringstream os_replaceTo;

        // set new function name
        std::string replacementName =
            callTerm->name() + "_" +
            std::to_string(nextFuncReplacementNumber++);
        os_replaceFrom << replacementName;
        os_replaceTo << callTerm->name() << "(";
        callTerm->setName(replacementName);

        // process every parameter
        bool pastFirstParam = false;
        int paramNumber = 0;
        bool haveAddedToOutput;
        bool haveAddedToInput;
        Exp* paramExp;
        // maintain a list of parameters that will remain in the call
        std::vector<Term*> termsToSave;
        unsigned int i;
        for (i = 0; i < callTerm->numArgs(); ++i) {
            // loop through all terms, adding them into the 'to' and 'from'
            // appropriately
            haveAddedToInput = false;
            haveAddedToOutput = false;
            if (pastFirstParam) {
                os_replaceTo << ",";
            }
            paramExp = callTerm->getParamExp(i);
            for (const auto& term : paramExp->getTermList()) {
                if (term->isConst()) {
                    // add the term to the function call, without making an
                    // input param for it
                    os_replaceTo << (haveAddedToOutput ? "+" : "") << "("
                                 << term->toString() << ")";
                } else {
                    // add the term to both the input and output function call
                    if (!haveAddedToInput && !pastFirstParam) {
                        os_replaceFrom << "(";
                    }
                    os_replaceFrom
                        << ((pastFirstParam || haveAddedToInput) ? "," : "")
                        << "p" << paramNumber;
                    os_replaceTo << (haveAddedToOutput ? "+" : "") << "p"
                                 << paramNumber;
                    termsToSave.push_back(term->clone());
                    paramNumber++;
                    haveAddedToInput = true;
                }
                haveAddedToOutput = true;
            }
            pastFirstParam = true;
        }

        if (termsToSave.size() != 0) {
            // rewrite argument list, one arg per term in original call args
            callTerm->resetNumArgs(termsToSave.size());
            i = 0;
            for (const auto& savedTerm : termsToSave) {
                Exp* newParamExp = new Exp();
                if (dynamic_cast<UFCallTerm*>(savedTerm)) {
                    // use a replacement variable, which will be constrained in
                    // this conjunction to be equal to the UF
                    VarTerm* replacementVar = new VarTerm(
                        savedTerm->coefficient(),
                        "rvar_" + std::to_string(nextVarReplacementNumber++));
                    newParamExp->addTerm(replacementVar);

                    // create a constraint in the current conjunction to make
                    // the replacement valid, for example: if we have a call
                    // A(B(i)), it will become A(rvar_0), and we will add the
                    // constraint B(i) - rvar_0 = 0
                    Exp* replacementExp = new Exp();
                    VarTerm* replacementVarForConstraint =
                        static_cast<VarTerm*>(replacementVar->clone());
                    replacementVarForConstraint->setCoefficient(-1);
                    savedTerm->setCoefficient(1);
                    replacementExp->addTerm(savedTerm);
                    replacementExp->addTerm(replacementVarForConstraint);

                    // add the constraint to a list to be added later, to avoid
                    // double-processing
                    UFSubstitutionConstraints.push_back(replacementExp);
                } else {
                    newParamExp->addTerm(savedTerm);
                }
                callTerm->setParamExp(i, newParamExp);

                i++;
            }
        } else {
            // replace 0-arg UF calls with symbolic constant
            VarTerm* replacementSymbol =
                new VarTerm(callTerm->coefficient(), callTerm->name());
            zeroArgUFReplacements.emplace(callTerm, replacementSymbol);
        }

        // complete outputs for this UF call
        if (haveAddedToInput) {
            os_replaceFrom << ")";
        }
        os_replaceTo << ")";
        macros << "#define " << os_replaceFrom.str() << " "
               << os_replaceTo.str() << "\n";
    }
};

/*!
 * \class Computation
 *
 * \brief SPF representation of a computation (group of statements such as a
 * function).
 */
class Computation {
   public:
    //! Construct an empty Computation
    Computation(){};

    //! Copy constructor
    Computation(const Computation& other);

    //! Assignment operator (copy)
    Computation& operator=(const Computation& other);

    //! Equality operator
    bool operator==(const Computation& other) const;

    //! Add a statement to this Computation and get a pointer to it.
    //! Statements are numbered sequentially from 0 as they are inserted
    Stmt* addStmt(const Stmt& stmt);
    //! Get a statement by index
    Stmt* getStmt(unsigned int index);

    //! Add a data space to this Computation
    void addDataSpace(std::string dataSpaceName);
    //! Get data spaces
    std::unordered_set<std::string> getDataSpaces() const;

    //! Get the number of statements in this Computation
    unsigned int getNumStmts() const;

    //! Print out all the information represented in this Computation for
    //! debug purposes
    void printInfo() const;

    //! Get whether this Computation's required information is filled out,
    //! including that of the Stmts it contains
    bool isComplete() const;

    //! Clear all data from this Computation
    void clear();

    //! Environment used by this Computation
    Environment env;

    //! Method generates c code.
    std::string codeGen();

   private:
    //! Information on all statements in the Computation
    std::vector<Stmt> stmts;
    //! Data spaces accessed in the computation
    std::unordered_set<std::string> dataSpaces;
};

/*!
 * \class Stmt
 *
 * \brief Information attached to a statement represented as mathematical
 * objects.
 */
class Stmt {
   public:
    //! Construct an empty Stmt
    Stmt(){};

    //! Construct a complete Stmt, given strings that will be used to
    //! construct each set/relation.
    Stmt(std::string stmtSourceCode, std::string iterationSpaceStr,
         std::string executionScheduleStr,
         std::vector<std::pair<std::string, std::string>> dataReadsStrs,
         std::vector<std::pair<std::string, std::string>> dataWritesStrs);

    //! Copy constructor
    Stmt(const Stmt& other);

    //! Assignment operator (copy)
    Stmt& operator=(const Stmt& other);

    //! Equality operator
    //! Checks equality, NOT mathematical equivalence
    bool operator==(const Stmt& other) const;

    //! Get whether or not all necessary information for this Stmt is set
    bool isComplete() const;

    //! Get the source code of the statement
    std::string getStmtSourceCode() const;
    //! Set the source code of the statement
    void setStmtSourceCode(std::string newStmtSourceCode);

    //! Get the iteration space Set
    Set* getIterationSpace() const;
    //! Set the iteration space, constructing it from the given string
    void setIterationSpace(std::string newIterationSpaceStr);

    //! Get the execution schedule Relation
    Relation* getExecutionSchedule() const;
    //! Set the execution schedule, constructing it from the given string
    void setExecutionSchedule(std::string newExecutionScheduleStr);

    //! Add a data read
    void addRead(std::string dataSpace, std::string relationStr);
    //! Get the number of data reads for this statement
    unsigned int getNumReads() const;
    //! Get a data read's data space by index
    std::string getReadDataSpace(unsigned int index) const;
    //! Get a data read's relation by index
    Relation* getReadRelation(unsigned int index) const;

    //! Add a data write
    void addWrite(std::string dataSpace, std::string relationStr);
    //! Get the number of data writes for this statement
    unsigned int getNumWrites() const;
    //! Get a data write's data space by index
    std::string getWriteDataSpace(unsigned int index) const;
    //! Get a data write's relation by index
    Relation* getWriteRelation(unsigned int index) const;

   private:
    //! Source code of the statement, for debugging purposes
    std::string stmtSourceCode;
    //! Iteration space of a statement
    std::unique_ptr<Set> iterationSpace;
    //! Execution schedule of a single statement
    std::unique_ptr<Relation> executionSchedule;
    //! Read dependences of a statement, pairing data space name to relation
    std::vector<std::pair<std::string, std::unique_ptr<Relation>>> dataReads;
    //! Write dependences of a statement, pairing data space name to relation
    std::vector<std::pair<std::string, std::unique_ptr<Relation>>> dataWrites;
};

}  // namespace iegenlib

#endif
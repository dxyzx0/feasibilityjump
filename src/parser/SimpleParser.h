//
// Created by psw on 5/19/24.
//

#ifndef FEASIBILITYJUMP_HEURISTIC_SIMPLEPARSER_H_
#define FEASIBILITYJUMP_HEURISTIC_SIMPLEPARSER_H_

#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <vector>
#include <cassert>
#include <algorithm>
#include "SimpleParser.h"
#include "DefaultCallback.h"
#include "PboCallback.h"
using namespace std;

/**
 * this class stores products of literals (as a tree) in order to
 * associate unique identifiers to these product (for linearization)
 */
template< typename Callback >
class ProductStore
{
 private:
	// we represent each node of a n-ary tree by a vector<ProductNode>
	struct ProductNode
	{
		int lit; // ID of the literal
		int productId; // identifier associated to the product of the
		// literals found from the root up to this node
		vector< ProductNode >* next; // list of next literals in a product

		ProductNode(int l)
		{
			lit = l;
			productId = 0;
			next = NULL;
		}

		// if we define a destructor to free <next>, we'll have to define
		// a copy constructor and use reference counting. It's not worth it.
	};

	vector< ProductNode > root; // root of the n-ary tree
	int nextSymbol{}; // next available variable

	/**
	 * define an order on ProductNode based on the literal (used to
	 * speed up look up)
	 */
	class ProductNodeLessLit
	{
	 public:
		bool operator()(const ProductNode& a, const ProductNode& b)
		{
			return a.lit < b.lit;
		}
	};

 public:
	/**
	 * give the first extra variable that can be used
	 */
	void setFirstExtraSymbol(int id)
	{
		nextSymbol = id;
	}

	/**
	 * get the identifier associated to a product term (update the list
	 * if necessary)
	 */
	int getProductVariable(vector< int >& list)
	{
		vector< ProductNode >* p = &root;
		typename vector< ProductNode >::iterator pos;

		// list must be sorted
		sort(list.begin(), list.end());

		// is this a known product ?
		for (size_t i = 0; i < list.size(); ++i)
		{
			assert(p != NULL);

			// look for list[i] in *p
			pos = lower_bound(p->begin(), p->end(), list[i], ProductNodeLessLit());
			if (pos == p->end() || (*pos).lit != list[i])
				pos = p->insert(pos, ProductNode(list[i])); // insert at the right place

			if (i != list.size() - 1 && (*pos).next == NULL)
				(*pos).next = new vector< ProductNode >;

			p = (*pos).next;
		}

		if ((*pos).productId == 0)
			(*pos).productId = nextSymbol++;

		return (*pos).productId;
	}

	/**
	 * add the constraints which define all product terms
	 *
	 */
	void defineProductVariable(Callback& cb)
	{
		vector< int > list;

		defineProductVariableRec(cb, root, list);
	}

	/**
	 * free all allocated product data
	 *
	 */
	void freeProductVariable()
	{
		freeProductVariableRec(root);
	}

 private:
	/**
	 * add the constraints which define all product terms
	 *
	 */
	void defineProductVariableRec(Callback& cb,
		vector< ProductNode >& nodes, vector< int >& list)
	{
		for (size_t i = 0; i < nodes.size(); ++i)
		{
			list.push_back(nodes[i].lit);
			if (nodes[i].productId)
				cb.linearizeProduct(nodes[i].productId, list);

			if (nodes[i].next)
				defineProductVariableRec(cb, *nodes[i].next, list);

			list.pop_back();
		}
	}

	/**
	 * free all allocated product data
	 *
	 */
	void freeProductVariableRec(vector< ProductNode >& nodes)
	{
		for (int i = 0; i < nodes.size(); ++i)
		{
			if (nodes[i].next)
			{
				freeProductVariableRec(*nodes[i].next);
				delete nodes[i].next;
			}
		}

		nodes.clear();
	}

};

template< typename Callback >
class SimpleParser
{
 public:
	Callback cb{};

 private:
	ifstream in; // the stream we're reading from
	int nbVars{}, nbConstr{}; // MetaData: #Variables and #Constraints in file.

	int nbProduct{}, sizeProduct{}; // MetaData for non linear format
	int nEqual{}; // MetaData for number of equal constraints
	int intSize{};
	ProductStore< Callback > store;
	bool autoLinearize{}; // should the parser linearize constraints ?

	/**
	 * get the next character from the stream
	 */
	char get()
	{
		return in.get();
	}

	/**
	 * put back a character into the stream (only one chr can be put back)
	 */
	void putback(char c)
	{
		in.putback(c);
	}

	/**
	 * return true iff we've reached EOF
	 */
	bool eof()
	{
		return !in.good();
	}

	/**
	 * skip white spaces
	 */
	void skipSpaces()
	{
		char c;

		while (isspace(c = get()));

		putback(c);
	}

	/**
	 * read an identifier from stream and append it to the list "list"
	 * @param list: the current list of identifiers that were read
	 * @return true iff an identifier was correctly read
	 */
	bool readIdentifier(vector< int >& list)
	{
		char c;
		bool negated = false;

		skipSpaces();

		// first char (must be 'x')
		c = get();
		if (eof())
			return false;

		if (c == '~')
		{
			negated = true;
			c = get();
		}

		if (c != 'x')
		{
			putback(c);
			return false;
		}

		int varID = 0;

		// next chars (must be digits)
		while (true)
		{
			c = get();
			if (eof())
				break;

			if (isdigit(c))
				varID = varID * 10 + c - '0';
			else
			{
				putback(c);
				break;
			}
		}

		//Small check on the coefficient ID to make sure everything is ok
		if (varID > nbVars)
			throw runtime_error("variable identifier larger than "
								"#variables in metadata.");

		if (negated)
			varID = -varID;

		list.push_back(varID);

		return true;
	}

	/**
	 * read a relational operator from stream and store it in s
	 * @param s: the variable to hold the relational operator we read
	 * @return true iff a relational operator was correctly read
	 */
	bool readRelOp(string& s)
	{
		char c;

		skipSpaces();

		c = get();
		if (eof())
			return false;

		if (c == '=')
		{
			s = "=";
			return true;
		}

		if (c == '>' && get() == '=')
		{
			s = ">=";
			return true;
		}

		return false;
	}

	/**
	 * read the first comment line to get the number of variables and
	 * the number of constraints in the file
	 *
	 * calls metaData with the data that was read
	 */
	void readMetaData()
	{
		char c;
		string s;

		// get the number of variables and constraints
		c = get();
		if (c != '*')
			throw runtime_error("First line of input file should be a comment");

		in >> s;
		if (eof() || s != "#variable=")
			throw runtime_error("First line should contain #variable= as first keyword");

		in >> nbVars;
		store.setFirstExtraSymbol(nbVars + 1);

		in >> s;
		if (eof() || s != "#constraint=")
			throw runtime_error("First line should contain #constraint= as second keyword");

		in >> nbConstr;

		skipSpaces();

		c = get();
		putback(c);

		if (c == '#')
		{
			// assume non linear format
			in >> s;
			if (eof() || (s != "#product=" && s != "#equal="))
				throw runtime_error("First line should contain #product= or #equal= as third keyword");
			else
			{
				if (s == "#equal=")
					in >> nEqual;
				else
					in >> nbProduct;
			}

			in >> s;
			if (eof() || (s != "sizeproduct=" && s != "intsize="))
				throw runtime_error("First line should contain sizeproduct= or intsize= as fourth keyword");
			else
			{
				if (s == "sizeproduct=")
					in >> sizeProduct;
				else
					in >> intSize;
			}
		}

		// skip the rest of the line
		getline(in, s);

		// callback to transmit the data
		if (nbProduct && autoLinearize)
		{
#ifdef ONLYCLAUSES
			cb.metaData(nbVars+nbProduct,nbConstr+nbProduct+sizeProduct);
#else
			cb.metaData(nbVars + nbProduct, nbConstr + 2 * nbProduct);
#endif
		}
		else
			cb.metaData(nbVars, nbConstr);

		// final check on the instance
		if (nbProduct > 0)
			throw runtime_error("nonlinear format not supported");

		// check if the intSize is less than 63
		if (intSize >= 64) {
			throw runtime_error("intsize should be less than 64");
		}
	}

	/**
	 * skip the comments at the beginning of the file
	 */
	void skipComments()
	{
		string s;
		char c;

		// skip further comments

		while (!eof() && (c = get()) == '*')
		{
			getline(in, s);
		}

		putback(c);
	}

	/**
	 * read a term into coeff and list
	 * @param coeff: the coefficient of the variable
	 * @param list: the list of literals identifiers in the product
	 */
	void readTerm(IntegerType& coeff, vector< int >& list)
	{
		list.clear();

		in >> coeff;

		skipSpaces();

		while (readIdentifier(list));

		if (list.size() == 0)
			throw runtime_error("identifier expected");
	}

	/**
	 * read the objective line (if any)
	 *
	 * calls beginObjective, objectiveTerm and endObjective
	 */
	void readObjective()
	{
		char c;
		string s;

		IntegerType coeff;
		vector< int > list;

		// read objective line (if any)

		skipSpaces();
		c = get();
		if (c != 'm')
		{
			// no objective line
			putback(c);
			return;
		}

		if (get() == 'i' && get() == 'n' && get() == ':')
		{
			cb.beginObjective(); // callback

			while (!eof())
			{
				readTerm(coeff, list);
				if (list.size() == 1 && list[0] > 0)
					cb.objectiveTerm(coeff, list[0]);
				else
					handleProduct(true, coeff, list);

				skipSpaces();
				c = get();
				if (c == ';')
					break; // end of objective
				else if (c == '-' || c == '+' || isdigit(c))
					putback(c);
				else
					throw runtime_error("unexpected character in objective function");
			}

			cb.endObjective();
		}
		else
			throw runtime_error("input format error: 'min:' expected");
	}

	/**
	 * read a constraint
	 *
	 * calls beginConstraint, constraintTerm and endConstraint
	 */
	void readConstraint()
	{
		string s;
		char c;

		IntegerType coeff;
		vector< int > list;

		cb.beginConstraint();

		while (!eof())
		{
			readTerm(coeff, list);
			if (list.size() == 1 && list[0] > 0)
				cb.constraintTerm(coeff, list[0]);
			else
				handleProduct(false, coeff, list);

			skipSpaces();
			c = get();
			if (c == '>' || c == '=')
			{
				// relational operator found
				putback(c);
				break;
			}
			else if (c == '-' || c == '+' || isdigit(c))
				putback(c);
			else
				throw runtime_error("unexpected character in constraint");
		}

		if (eof())
			throw runtime_error("unexpected EOF before end of constraint");

		if (readRelOp(s))
			cb.constraintRelOp(s);
		else
			throw runtime_error("unexpected relational operator in constraint");

		in >> coeff;
		cb.constraintRightTerm(coeff);

		skipSpaces();
		c = get();
		if (eof() || c != ';')
			throw runtime_error("semicolon expected at end of constraint");

		cb.endConstraint();
	}

	/**
	 * passes a product term to the solver (first linearizes the product
	 * if this is wanted)
	 */
	void handleProduct(bool inObjective, IntegerType coeff, vector< int >& list)
	{
		if (autoLinearize)
		{
			// get symbol corresponding to this product
			int var = store.getProductVariable(list);

			if (inObjective)
				cb.objectiveTerm(coeff, var);
			else
				cb.constraintTerm(coeff, var);
		}
		else
		{
			if (inObjective)
				cb.objectiveProduct(coeff, list);
			else
				cb.constraintProduct(coeff, list);
		}
	}

 public:
	/**
	 * constructor which only opens the file
	 */
	explicit SimpleParser(const char* filename)
	{
		in.open(filename, ios_base::in);

		if (!in.good())
			throw runtime_error("error opening input file");
	}

	~SimpleParser()
	{
		store.freeProductVariable();
	}

	/**
	 * ask the parser to linearize the constraints using a simple
	 * default method
	 */
	void setAutoLinearize(bool lin = true)
	{
		autoLinearize = lin;
	}

	/**
	 * parses the file and uses the callbacks to send the data
	 * back to the program
	 */
	void parse()
	{
		char c;

		readMetaData();

		skipComments();

		readObjective();

		// read constraints
		int nbConstraintsRead = 0;
		while (!eof())
		{
			skipSpaces();
			if (eof())
				break;

			putback(c = get());
			if (c == '*')
				skipComments();

			if (eof())
				break;

			readConstraint();
			nbConstraintsRead++;
		}

		//Small check on the number of constraints
		if (nbConstraintsRead != nbConstr)
			throw runtime_error("number of constraints read is different from metadata.");

		if (autoLinearize)
		{
			store.defineProductVariable(cb);
		}
	}
}; // class SimpleParser

#endif //FEASIBILITYJUMP_HEURISTIC_SIMPLEPARSER_H_

/*=============================================================================
 * parser for pseudo-Boolean instances
 * 
 * Copyright (c) 2005-2007 Olivier ROUSSEL and Vasco MANQUINHO
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *=============================================================================
 */

/* version 2.9.4 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#ifdef useGMP
#include <gmp.h>

//#define debug

typedef 
  mpz_t IntegerType;
#else
#warning this IntegerType may not be suitable for some input file. Consider using GMP
typedef
long IntegerType;
#endif

/**
 * stop the program when an error is found
 *
 * @param msg: the error message
 */
void runtimeError(char* msg)
{
	printf("Error: %s\n", msg);
	exit(1);
}

/*
 * This section defines the callback that are used by the
 * parser to send data back to the program.
 *
 * These functions must be redefined by the user.
 */

/**
 * callback called when we get the number of variables and the
 * expected number of constraints
 * @param nbvar: the number of variables
 * @param nbconstr: the number of contraints
 */
void metaData(int nbvar, int nbconstr)
{
	printf("[nbvar=%d]\n", nbvar);
	printf("[nbconstr=%d]\n", nbconstr);
}

/**
 * callback called before we read the objective function
 */
void beginObjective()
{
	printf("objective:  ");
}

/**
 * callback called after we've read the objective function
 */
void endObjective()
{
	printf("\n");
}

/**
 * callback called when we read a term of the objective function
 *
 * @param coeff: the coefficient of the term
 * @param idVar: the numerical identifier of the variable
 */
void objectiveTerm(IntegerType coeff, int idVar)
{
#ifdef useGMP
	gmp_printf("[%+Zd x%d] ",coeff,idVar);
#else
	printf("[%+d x%d] ", coeff, idVar);
#endif
}

/**
 * callback called when we read a term of the objective function
 * which is a product of literals
 *
 * @param coeff: the coefficient of the term
 * @param list: list of literals which appear in the product
 * @param n: number of literals in the list
 */
void objectiveProduct(IntegerType coeff, int* list, int n)
{
	int i;

#ifdef useGMP
	gmp_printf("[%+Zd ",coeff);
#else
	printf("[%+d ", coeff);
#endif

	for (i = 0; i < n; ++i)
	{
		if (list[i] < 0)
			printf("~x%d ", -list[i]);
		else
			printf("x%d ", list[i]);
	}

	printf("] ");
}

/**
 * callback called before we read a constraint
 */
void beginConstraint()
{
	printf("constraint: ");
}

/**
 * callback called after we've read a constraint
 */
void endConstraint()
{
	printf("\n");
}

/**
 * callback called when we read a term of the constraint
 *
 * @param coeff: the coefficient of the term
 * @param idVar: the numerical identifier of the variable
 */
void constraintTerm(IntegerType coeff, int idVar)
{
#ifdef useGMP
	gmp_printf("[%+Zd x%d] ",coeff,idVar);
#else
	printf("[%+d x%d] ", coeff, idVar);
#endif
}

/**
 * callback called when we read a term of the constraint
 * which is a product of literals
 *
 * @param coeff: the coefficient of the term
 * @param list: list of literals which appear in the product
 * @param n: number of literals in the list
 */
void constraintProduct(IntegerType coeff, int* list, int n)
{
	int i;

#ifdef useGMP
	gmp_printf("[%+Zd ",coeff);
#else
	printf("[%+d ", coeff);
#endif

	for (i = 0; i < n; ++i)
	{
		if (list[i] < 0)
			printf("~x%d ", -list[i]);
		else
			printf("x%d ", list[i]);
	}

	printf("] ");
}

/**
 * callback called when we read the relational operator of a constraint
 *
 * @param relop: the relational operator (>= or =)
 */
void constraintRelOp(char* relop)
{
	printf("[%s] ", relop);
}

/**
 * callback called when we read the right term of a constraint (also
 * known as the degree)
 *
 * @param val: the degree of the constraint
 */
void constraintRightTerm(IntegerType val)
{
#ifdef useGMP
	gmp_printf("[%+Zd] ",val);
#else
	printf("[%+d] ", val);
#endif
}

/**
 * add the necessary constraints to define newSymbol as equivalent
 * to the product (conjunction) of literals in product.
 */
void linearizeProduct(int newSymbol, int* product, int n)
{
	IntegerType r, one, minusOne;
	int i;

#ifdef debug
	printf("linearize(%d,{",newSymbol);
	for(i=0;i<n;++i)
	  printf("%d ",product[i]);
	printf("})\n");
#endif

#ifdef useGMP
	mpz_init(r);
	mpz_init_set_si(one,1);
	mpz_init_set_si(minusOne,-1);
#else
	one = 1;
	minusOne = -1;
#endif

	// product => newSymbol (this is a clause)
	// not x0 or not x1 or ... or not xn or newSymbol
#ifdef useGMP
	mpz_set_si(r,1);
#else
	r = 1;
#endif
	beginConstraint();
	constraintTerm(one, newSymbol);
	for (i = 0; i < n; ++i)
		if (product[i] > 0)
		{
			constraintTerm(minusOne, product[i]);
#ifdef useGMP
			mpz_sub_ui(r,r,1);
#else
			r -= 1;
#endif
		}
		else
			constraintTerm(one, -product[i]);
	constraintRelOp(">=");
	constraintRightTerm(r);
	endConstraint();

#ifdef ONLYCLAUSES
	// newSymbol => product translated as
	// not newSymbol of xi (for all i)
	for(i=0;i<n;++i)
	{
#ifdef useGMP
	  mpz_set_si(r,0);
#else
	  r=0;
#endif
	  beginConstraint();
	  constraintTerm(minusOne,newSymbol);
	  if (product[i]>0)
		constraintTerm(one,product[i]);
	  else
	  {
		constraintTerm(minusOne,-product[i]);
#ifdef useGMP
		mpz_sub_ui(r,r,1);
#else
		r-=1;
#endif
	  }
	  constraintRelOp(">=");
	  constraintRightTerm(r);
	  endConstraint();
	}
#else
	// newSymbol => product translated as
	// x0+x1+x3...+xn-n*newSymbol>=0
	beginConstraint();
#ifdef useGMP
	mpz_set_si(r,-n);
#else
	r = -n;
#endif
	constraintTerm(r, newSymbol);
#ifdef useGMP
	mpz_set_si(r,0);
#else
	r = 0;
#endif
	for (i = 0; i < n; ++i)
		if (product[i] > 0)
			constraintTerm(one, product[i]);
		else
		{
			constraintTerm(minusOne, -product[i]);
#ifdef useGMP
			mpz_sub_ui(r,r,1);
#else
			r -= 1;
#endif
		}
	constraintRelOp(">=");
	constraintRightTerm(r);
	endConstraint();
#endif

#ifdef useGMP
	mpz_clear(r);
	mpz_clear(one);
	mpz_clear(minusOne);
#endif
}


/*
 * private section
 *
 */

/*
 * a tiny expandable array to represent a list of integers
 */

typedef
struct
{
	int* tab; /* current array containing the list */
	int capacity; /* current capacity of the array */
	int size; /* current number of elements in the array */
} Vector_Int;

void Vector_Int_constructor(Vector_Int* this)
{
	this->size = 0;
	this->capacity = 16;
	this->tab = calloc(this->capacity, sizeof(int));
	if (!this->tab)
		runtimeError("failed to allocate memory for a list");
}

void Vector_Int_clear(Vector_Int* this)
{
	this->size = 0;
}

void Vector_Int_push_back(Vector_Int* this, int x)
{
	if (this->size + 1 >= this->capacity)
	{
		this->capacity += 16;
		this->tab = realloc(this->tab, this->capacity * sizeof(int));
		if (!this->tab)
			runtimeError("failed to expand memory for a Vector_Int");
	}

	this->tab[this->size] = x;
	this->size++;
}

void Vector_Int_pop_back(Vector_Int* this)
{
	this->size--;
}

void Vector_Int_print(Vector_Int* this)
{
	int i;

	printf("{");

	if (this->size)
		printf("%d", this->tab[0]);

	for (i = 1; i < this->size; ++i)
		printf(",%d", this->tab[i]);

	printf("}");
}

void Vector_Int_destructor(Vector_Int* this)
{
	free(this->tab);
}

// for qsort()
int compareInt(const void* a, const void* b)
{
	return *(int*)a - *(int*)b;
}

/**
 * the pseudo-class ProductStore stores products of literals (as a
 * tree) in order to associate unique identifiers to these product
 * (for linearization)
 */

typedef
struct
{
	int lit; // ID of the literal
	int productId; // identifier associated to the product of the
	// literals found from the root up to this node
	struct _Vector_ProductNode* next; // list of next literals in a product
} ProductNode;

void ProductNode_constructor(ProductNode* this, int l);

// a node of the tree
typedef
struct _Vector_ProductNode
{
	ProductNode* tab; // current array
	int size; // current number of elements in tab
	int capacity; // current capacity of the array
} Vector_ProductNode;

void Vector_ProductNode_constructor(Vector_ProductNode* this)
{
	this->size = 0;
	this->capacity = 16;
	this->tab = calloc(this->capacity, sizeof(ProductNode));
}

void Vector_ProductNode_clear(Vector_ProductNode* this)
{
	this->size = 0;
}

Vector_ProductNode* new_Vector_ProductNode()
{
	Vector_ProductNode* p = malloc(sizeof(Vector_ProductNode));
	Vector_ProductNode_constructor(p);
	return p;
}

void ProductNode_constructor(ProductNode* this, int l)
{
	this->lit = l;
	this->productId = 0;
	this->next = NULL;
}

void Vector_ProductNode_insert(Vector_ProductNode* this,
	int pos, int lit)
{
	int i;

	if (this->size + 1 >= this->capacity)
	{
		this->capacity += 16;
		this->tab = realloc(this->tab, this->capacity * sizeof(ProductNode));
		if (!this->tab)
			runtimeError("failed to expand memory for a Vector_ProductNode");
	}

	for (i = this->size; i > pos; --i)
		this->tab[i] = this->tab[i - 1];

	ProductNode_constructor(&this->tab[pos], lit);
	++this->size;
}

/**
 * return the position of elem in l (binary search)
 */
int Vector_ProductNode_binary_search(Vector_ProductNode l, int elem)
{
	// binary search in [beg..end[
	int beg = 0, end = l.size, mid;

	while (beg < end)
	{
		mid = (beg + end) / 2;
		if (l.tab[mid].lit == elem)
			return mid;
		else if (l.tab[mid].lit > elem)
			end = mid;
		else
			beg = mid + 1;
	}

	// not found
	return beg;
}

typedef
struct
{
	Vector_ProductNode root; // root of the n-ary tree
	int nextSymbol; // next available variable
} ProductStore;

void ProductStore_constructor(ProductStore* this)
{
	Vector_ProductNode_constructor(&this->root);
	this->nextSymbol = 0;
}

/**
 * give the first extra variable that can be used
 */
void ProductStore_setFirstExtraSymbol(ProductStore* this, int id)
{
	this->nextSymbol = id;
}

/**
 * get the identifier associated to a product term (update the list
 * if necessary)
 */
int ProductStore_getProductVariable(ProductStore* this, Vector_Int* list)
{
	Vector_ProductNode* pn = &this->root;
	ProductNode* pNode = NULL;
	int i, pos;

#ifdef debug
	printf("getProductVariable for ");
	Vector_Int_print(list);
	printf("\n");
#endif

	// list must be sorted
	qsort(list->tab, list->size, sizeof(int), compareInt);

	// is this a known product ?
	for (i = 0; i < list->size; ++i)
	{
		assert(pn != NULL);

		// look for list[i] in *pn
		pos = Vector_ProductNode_binary_search(*pn, list->tab[i]);
		if (pos >= pn->size || pn->tab[pos].lit != list->tab[i])
			Vector_ProductNode_insert(pn, pos, list->tab[i]); // insert at the right place

		pNode = &pn->tab[pos];
		if (i != list->size - 1 && pNode->next == NULL)
			pNode->next = new_Vector_ProductNode();

		pn = pNode->next;
	}

	if (pNode->productId == 0)
		pNode->productId = this->nextSymbol++;

	return pNode->productId;
}

/**
 * add the constraints which define all product terms
 *
 */
void ProductStore_defineProductVariableRec(Vector_ProductNode* nodes,
	Vector_Int* list)
{
	int i;
	for (i = 0; i < nodes->size; ++i)
	{
		Vector_Int_push_back(list, nodes->tab[i].lit);
		if (nodes->tab[i].productId)
			linearizeProduct(nodes->tab[i].productId, list->tab, list->size);

		if (nodes->tab[i].next)
			ProductStore_defineProductVariableRec(nodes->tab[i].next, list);

		Vector_Int_pop_back(list);
	}
}

/**
 * free all allocated product data
 *
 */
void ProductStore_freeProductVariableRec(Vector_ProductNode* nodes)
{
	int i;
	for (i = 0; i < nodes->size; ++i)
	{
		if (nodes->tab[i].next)
		{
			ProductStore_freeProductVariableRec(nodes->tab[i].next);
			free(nodes->tab[i].next);
		}
	}

	Vector_ProductNode_clear(nodes);
}

/**
 * add the constraints which define all product terms
 *
 */
void ProductStore_defineProductVariable(ProductStore* this)
{
	Vector_Int list;
	Vector_Int_constructor(&list);

	ProductStore_defineProductVariableRec(&this->root, &list);

	Vector_Int_destructor(&list);
}

/**
 * free all allocated product data
 *
 */
void ProductStore_freeProductVariable(ProductStore* this)
{
	ProductStore_freeProductVariableRec(&this->root);
}

FILE* in; // the file we're reading from
int nbVars, nbConstr; // MetaData: #Variables and #Constraints in file.

int nbProduct, sizeProduct; // MetaData for non linear format
char autoLinearize = 0; // should the parser linearize constraints ?
ProductStore store;

/**
 * get the next character from the stream
 */
char get()
{
	return fgetc(in);
}

/**
 * put back a character into the stream (only one chr can be put back)
 */
void putback(char c)
{
	ungetc(c, in);
}

/**
 * return 1 iff we've reached EOF
 */
int eof()
{
	return feof(in);
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
 * @return 1 iff an identifier was correctly read
 */
int readIdentifier(Vector_Int* list)
{
	char c;
	int varID;
	int negated = 0;

	skipSpaces();

	// first char (must be 'x')
	c = get();
	if (eof())
		return 0;

	if (c == '~')
	{
		negated = 1;
		c = get();
	}

	if (c != 'x')
	{
		putback(c);
		return 0;
	}

	varID = 0;

	//next chars (must be digits)
	while (1)
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
		runtimeError("Variable identifier larger than #variables in metadata.");

	if (negated)
		varID = -varID;

	Vector_Int_push_back(list, varID);

	return 1;
}

/**
 * read a relational operator from stream and store it in s
 * @param s: the variable to hold the relational operator we read. Must be at lest 3 characters long.
 * @return 1 iff a relational operator was correctly read
 */
int readRelOp(char* s)
{
	char c;

	skipSpaces();

	c = get();
	if (eof())
		return 0;

	if (c == '=')
	{
		strcpy(s, "=");
		return 1;
	}

	if (c == '>' && get() == '=')
	{
		strcpy(s, ">=");
		return 1;
	}

	return 0;
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
	char s[1024];

	// get the number of variables and constraints
	c = get();
	if (c != '*')
		runtimeError("First line of input file should be a comment");

	fscanf(in, "%20s", s);
	if (eof() || strcmp(s, "#variable=") != 0)
		runtimeError("First line should contain #variable= as first keyword");

	fscanf(in, "%d", &nbVars);

	ProductStore_setFirstExtraSymbol(&store, nbVars + 1);

	fscanf(in, "%20s", s);
	if (eof() || strcmp(s, "#constraint=") != 0)
		runtimeError("First line should contain #constraint= as second keyword");

	fscanf(in, "%d", &nbConstr);

	skipSpaces();

	c = get();
	putback(c);

	if (c == '#')
	{
		// assume non linear format
		fscanf(in, "%20s", s);
		if (eof() || strcmp(s, "#product=") != 0)
			runtimeError("First line should contain #product= as third keyword");

		fscanf(in, "%d", &nbProduct);

		fscanf(in, "%20s", s);
		if (eof() || strcmp(s, "sizeproduct=") != 0)
			runtimeError("First line should contain sizeproduct= as fourth keyword");

		fscanf(in, "%d", &sizeProduct);
	}

	// skip the rest of the line
	fgets(s, sizeof(s), in);

	// callback to transmit the data
	if (nbProduct && autoLinearize)
	{
#ifdef ONLYCLAUSES
		metaData(nbVars+nbProduct,nbConstr+nbProduct+sizeProduct);
#else
		metaData(nbVars + nbProduct, nbConstr + 2 * nbProduct);
#endif
	}
	else
		metaData(nbVars, nbConstr);
}

/**
 * skip the comments at the beginning of the file
 */
void skipComments()
{
	char s[1024];
	char c;

	// skip further comments

	while (!eof() && (c = get()) == '*')
	{
		fgets(s, sizeof(s), in);
	}

	putback(c);
}

/**
 * read a term into coeff and var
 * @param coeff: the coefficient of the variable
 * @param list: the list of literals identifiers in the product
 */
void readTerm(IntegerType* coeff, Vector_Int* list)
{
	char c;

	Vector_Int_clear(list);

#ifdef useGMP
	gmp_fscanf(in,"%Zd",*coeff);
#else
	fscanf(in, "%d", coeff);
#endif

	skipSpaces();

	while (readIdentifier(list));

	if (list->size == 0)
		runtimeError("identifier expected");
}

/**
 * passes a product term to the solver (first linearizes the product
 * if this is wanted)
 */
void handleProduct(char inObjective, IntegerType coeff, Vector_Int* list)
{
	if (autoLinearize)
	{
		// get symbol corresponding to this product
		int var = ProductStore_getProductVariable(&store, list);

		if (inObjective)
			objectiveTerm(coeff, var);
		else
			constraintTerm(coeff, var);
	}
	else
	{
		if (inObjective)
			objectiveProduct(coeff, list->tab, list->size);
		else
			constraintProduct(coeff, list->tab, list->size);
	}
}

/**
 * read the objective line (if any)
 *
 * calls beginObjective, objectiveTerm and endObjective
 */
void readObjective()
{
	char c;
	char s[1024];
	IntegerType coeff;
	Vector_Int list;

	Vector_Int_constructor(&list);

#ifdef useGMP
	mpz_init(coeff);
#endif

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
		beginObjective(); // callback

		while (!eof())
		{
			readTerm(&coeff, &list);
			if (list.size == 1 && list.tab[0] > 0)
				objectiveTerm(coeff, list.tab[0]);
			else
				handleProduct(1, coeff, &list);

			skipSpaces();
			c = get();
			if (c == ';')
				break; // end of objective
			else if (c == '-' || c == '+' || isdigit(c))
				putback(c);
			else
				runtimeError("unexpected character in objective function");
		}

		endObjective();
	}
	else
		runtimeError("input format error: 'min:' expected");

#ifdef useGMP
	mpz_clear(coeff);
#endif
}

/**
 * read a constraint
 *
 * calls beginConstraint, constraintTerm and endConstraint
 */
void readConstraint()
{
	char s[1024];
	char c;

	IntegerType coeff;
	Vector_Int list;

	Vector_Int_constructor(&list);

#ifdef useGMP
	mpz_init(coeff);
#endif

	beginConstraint();

	while (!eof())
	{
		readTerm(&coeff, &list);
		if (list.size == 1 && list.tab[0] > 0)
			constraintTerm(coeff, list.tab[0]);
		else
			handleProduct(0, coeff, &list);

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
			runtimeError("unexpected character in constraint");
	}

	if (eof())
		runtimeError("unexpected EOF before end of constraint");

	if (readRelOp(s))
		constraintRelOp(s);
	else
		runtimeError("unexpected relational operator in constraint");

#ifdef useGMP
	gmp_fscanf(in,"%Zd",coeff);
#else
	fscanf(in, "%d", &coeff);
#endif
	constraintRightTerm(coeff);

	skipSpaces();
	c = get();
	if (eof() || c != ';')
		runtimeError("semicolon expected at end of constraint");

	endConstraint();

#ifdef useGMP
	mpz_clear(coeff);
#endif
}


/*
 * public section
 *
 */

/**
 * parses the file and uses the callbacks to send to send the data
 * back to the program
 *
 * @param filename: the file to parse
 */
void parse(char* filename)
{
	char c;
	in = fopen(filename, "rt");

	if (in == NULL)
		runtimeError("error opening input file");

	ProductStore_constructor(&store);

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
	{
		runtimeError("Number of constraints read is different from metadata.");
	}

	if (autoLinearize)
	{
		ProductStore_defineProductVariable(&store);
	}

}

int main(int argc, char* argv[])
{
	if (argc != 2)
		printf("usage: SimpleParser <filename>\n");
	else
	{
		autoLinearize = 1;

		parse(argv[1]);
	}

	return 0;
}

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

// version 2.9.0

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.math.BigInteger;
import java.util.Vector;
import java.util.Iterator;
import java.lang.Integer;

class SimpleParser {
    /*
     * This protected section defines the callback that are used by the parser
     * to send data back to the program.
     * 
     * These functions must be redefined in a subclass.
     */

    /**
     * callback called when we get the number of variables and the expected
     * number of constraints
     * 
     * @param nbvar:
     *            the number of variables
     * @param nbconstr:
     *            the number of contraints
     */
    protected void metaData(int nbvar, int nbconstr) {
        System.out.println("[nbvar=" + nbvar + "]");
        System.out.println("[nbconstr=" + nbconstr + "]");
    }

    /**
     * callback called before we read the objective function
     */
    protected void beginObjective() {
        System.out.print("objective:  ");
    }

    /**
     * callback called after we've read the objective function
     */
    protected void endObjective() {
        System.out.println();
    }

    /**
     * callback called when we read a term of the objective function
     * 
     * @param coeff:
     *            the coefficient of the term
     * @param idVar: 
     *            the numerical identifier of the variable
     */
    protected void objectiveTerm(BigInteger coeff, int idVar) {
        System.out.print("[" + (coeff.signum() == 1 ? "+" : "") + coeff + " x"
                + idVar + "] ");
    }

   /**
     * callback called when we read a term of the objective function
     * which is a product of literals
     * 
     * @param coeff:
     *            the coefficient of the term
     * @param list: 
     *            list of literals which appear in the product
     */
    protected void objectiveProduct(BigInteger coeff, Vector list) {
        System.out.print("[" + (coeff.signum() == 1 ? "+" : "") + coeff + " ");

	Iterator i=list.iterator();

	while (i.hasNext())
	{
	    int id=((Integer)i.next()).intValue();
	    
	    if (id<0)
		System.out.print("~x"+(-id)+" ");
	    else
		System.out.print("x"+id+" ");
	}

        System.out.print("] ");
    }

    /**
     * callback called before we read a constraint
     */
    protected void beginConstraint() {
        System.out.print("constraint: ");
    }

    /**
     * callback called after we've read a constraint
     */
    protected void endConstraint() {
        System.out.println();
    }

    /**
     * callback called when we read a term of a constraint
     * 
     * @param coeff:
     *            the coefficient of the term
     * @param idVar:
     *            the numerical identifier of the variable
     */
    protected void constraintTerm(BigInteger coeff, int idVar) {
        System.out.print("[" + (coeff.signum() == 1 ? "+" : "") + coeff + " x"
                + idVar + "] ");
    }

    /**
     * callback called when we read a term of a constraint
     * which is a product of literals
     * 
     * @param coeff:
     *            the coefficient of the term
     * @param list: 
     *            list of literals which appear in the product
     */
    protected void constraintProduct(BigInteger coeff, Vector list) {
        System.out.print("[" + (coeff.signum() == 1 ? "+" : "") + coeff + " ");

	Iterator i=list.iterator();

	while (i.hasNext())
	{
	    int id=((Integer)i.next()).intValue();
	    
	    if (id<0)
		System.out.print("~x"+(-id)+" ");
	    else
		System.out.print("x"+id+" ");
	}

        System.out.print("] ");
    }

    /**
     * callback called when we read the relational operator of a constraint
     * 
     * @param relop:
     *            the relational oerator (>= or =)
     */
    protected void constraintRelOp(String relop) {
        System.out.print("[" + relop + "] ");
    }

    /**
     * callback called when we read the right term of a constraint (also known
     * as the degree)
     * 
     * @param val:
     *            the degree of the constraint
     */
    protected void constraintRightTerm(BigInteger val) {
        System.out.print("[" + (val.signum() == 1 ? "+" : "") + val + "]");
    }

    BufferedReader in; // the stream we're reading from

    char savedChar; // a character read from the file but not yet consumed

    boolean charAvailable = false; // true iff savedChar contains a character

    boolean eofReached = false; // true iff we've reached EOF
    
    int nbVars, nbConstr; // MetaData: #Variables and #Constraints in file.
    
    public class ParseException extends Exception {
        /**
         * Constructor for ParseException.
         * 
         * @param message
         *            the error message
         */
        public ParseException(String message) {
            super(message);
        }
    }

    /**
     * get the next character from the stream
     * 
     * @throws IOException
     */
    private char get() throws IOException {
        int c;

        if (charAvailable) {
            charAvailable = false;
            return savedChar;
        }

        c = in.read();
        if (c == -1)
            eofReached = true;

        return (char) c;
    }

    /**
     * put back a character into the stream (only one chr can be put back)
     */
    private void putback(char c) {
        savedChar = c;
        charAvailable = true;
    }

    /**
     * return true iff we've reached EOF
     */
    private boolean eof() {
        return eofReached;
    }

    /**
     * skip white spaces
     * 
     * @throws IOException
     */
    private void skipSpaces() throws IOException {
        char c;

        while (Character.isWhitespace(c = get()))
            ;

        putback(c);
    }

    /**
     * read a word from file
     * 
     * @return the word we read
     * @throws IOException
     */
    public String readWord() throws IOException {
        StringBuffer s = new StringBuffer();
        char c;

        skipSpaces();

        while (!Character.isWhitespace(c = get()) && !eof())
            s.append(c);

        return s.toString();
    }

    /**
     * read a integer from file
     * 
     * @param s:
     *            a StringBuffer to store the integer that was read
     * @throws IOException
     */
    public void readInteger(StringBuffer s) throws IOException {
        char c;

        skipSpaces();
        s.setLength(0);

        c = get();
        if (c == '-' || Character.isDigit(c))
            s.append(c);
        // note: BigInteger don't like a '+' before the number, we just skip it

        while (Character.isDigit(c = get()) && !eof())
            s.append(c);

        putback(c);
    }

    /**
     * read an identifier from stream and store it in s
     * 
     * @return the identifier we read or null
     * @throws IOException, ParseException
     */
    private boolean readIdentifier(Vector list) throws IOException, ParseException {
	StringBuffer s=new StringBuffer();
        char c;
	boolean negated=false;
	
        skipSpaces();

        // first char (must be a 'x')
        c = get();
        if (eof())
            return false;
	
	if (c == '~')
	{
	    negated=true;
	    c=get();
	}

        if (c != 'x') {
            putback(c);
            return false;
        }

        s.append(c);

        // next chars (must be digits)
        while (true) {
            c = get();
            if (eof())
                break;

            if (Character.isDigit(c))
                s.append(c);
            else {
                putback(c);
                break;
            }
        }
	
	//Small check on the coefficient ID to make sure everything is ok
	int varID = Integer.parseInt(s.substring(1));
	if (varID > nbVars) {
	    throw new ParseException("Variable identifier larger than #variables in metadata.");
	}
	
	if (negated)
	    varID=-varID;

	list.add(new Integer(varID));
	
        return true;
    }

    /**
     * read a relational operator from stream and store it in s
     * 
     * @return the relational operator we read or null
     * @throws IOException
     */
    private String readRelOp() throws IOException {
        char c;

        skipSpaces();

        c = get();
        if (eof())
            return null;

        if (c == '=')
            return "=";

        if (c == '>' && get() == '=')
            return ">=";

        return null;
    }

    /**
     * read the first comment line to get the number of variables and the number
     * of constraints in the file
     * 
     * calls metaData with the data that was read
     * 
     * @throws IOException
     * @throws ParseException
     */
    private void readMetaData() throws IOException, ParseException {
        char c;
        String s;

        // get the number of variables and constraints
        c = get();
        if (c != '*')
            throw new ParseException(
                    "First line of input file should be a comment");

        s = readWord();
        if (eof() || !s.equals("#variable="))
            throw new ParseException(
                    "First line should contain #variable= as first keyword");

        nbVars = Integer.parseInt(readWord());

        s = readWord();
        if (eof() || !s.equals("#constraint="))
            throw new ParseException(
                    "First line should contain #constraint= as second keyword");

        nbConstr = Integer.parseInt(readWord());

        // skip the rest of the line
        in.readLine();

        // callback to transmit the data
        metaData(nbVars, nbConstr);
    }

    /**
     * skip the comments at the beginning of the file
     * 
     * @throws IOException
     */
    private void skipComments() throws IOException {
        char c = ' ';

        // skip further comments

        while (!eof() && (c = get()) == '*') {
            in.readLine();
        }

        putback(c);
    }

    /**
     * read a term into coeff and var
     * 
     * @param coeff:
     *            the coefficient of the variable
     * @param list:
     *            the list of literals identifiers in the product
     * @throws IOException
     * @throws ParseException
     */
    private void readTerm(StringBuffer coeff, Vector list)
            throws IOException, ParseException {
        char c;

	list.clear();

        readInteger(coeff);

        skipSpaces();

        while(readIdentifier(list));

	if (list.size()==0)
            throw new ParseException("identifier expected");
    }

    /**
     * read the objective line (if any)
     * 
     * calls beginObjective, objectiveTerm and endObjective
     * 
     * @throws IOException
     * @throws ParseException
     */
    private void readObjective() throws IOException, ParseException {
        char c;
        Vector list = new Vector();
        StringBuffer coeff = new StringBuffer();

        // read objective line (if any)

        skipSpaces();
        c = get();
        if (c != 'm') {
            // no objective line
            putback(c);
            return;
        }

        if (get() == 'i' && get() == 'n' && get() == ':') {
            beginObjective(); // callback

            while (!eof()) {
                readTerm(coeff, list);
		if (list.size()==1)
		{
		    int val=((Integer)list.firstElement()).intValue();
		    if (val>0)
			objectiveTerm(new BigInteger(coeff.toString()), 
				      val); // callback
		    else
			objectiveProduct(new BigInteger(coeff.toString()), 
					 list); // callback
		}
		else
		    objectiveProduct(new BigInteger(coeff.toString()), 
				     list); // callback

                skipSpaces();
                c = get();
                if (c == ';')
                    break; // end of objective
                else if (c == '-' || c == '+' || Character.isDigit(c))
                    putback(c);
                else
                    throw new ParseException(
                            "unexpected character in objective function");
            }

            endObjective();
        } else
            throw new ParseException("input format error: 'min:' expected");
    }

    /**
     * read a constraint
     * 
     * calls beginConstraint, constraintTerm and endConstraint
     * 
     * @throws ParseException
     * @throws IOException
     */
    private void readConstraint() throws IOException, ParseException {
        Vector list = new Vector();
        StringBuffer coeff = new StringBuffer();
        char c;

        beginConstraint();

        while (!eof()) {
	    readTerm(coeff, list);
	    if (list.size()==1) {
		int val=((Integer)list.firstElement()).intValue();
		if (val>0)
		    constraintTerm(new BigInteger(coeff.toString()), 
				   val); // callback
		else
		    constraintProduct(new BigInteger(coeff.toString()), 
				      list); // callback
	    }
	    else
		constraintProduct(new BigInteger(coeff.toString()), 
				  list); // callback

            skipSpaces();
            c = get();
            if (c == '>' || c == '=') {
                // relational operator found
                putback(c);
                break;
            } else if (c == '-' || c == '+' || Character.isDigit(c))
                putback(c);
            else
                throw new ParseException("unexpected character in constraint");
        }

        if (eof())
            throw new ParseException("unexpected EOF before end of constraint");

        String relop;
        if ((relop = readRelOp()) != null)
            constraintRelOp(relop);
        else
            throw new ParseException(
                    "unexpected relational operator in constraint");

        readInteger(coeff);
        constraintRightTerm(new BigInteger(coeff.toString()));

        skipSpaces();
        c = get();
        if (eof() || c != ';')
            throw new ParseException("semicolon expected at end of constraint");

        endConstraint();
    }

    /**
     * constructor which only opens the file
     * 
     * @throws FileNotFoundException
     */
    public SimpleParser(String filename) throws FileNotFoundException {
        in = new BufferedReader(new InputStreamReader(new FileInputStream(
                filename)));
    }

    /**
     * parses the file and uses the callbacks to send to send the data back to
     * the program
     * 
     * @throws IOException
     * @throws ParseException
     */
    public void parse() throws IOException, ParseException {
	char c;

        readMetaData();

        skipComments();

        readObjective();

        // read constraints
	int nbConstraintsRead = 0;
        while (!eof()) {
            skipSpaces();
            if (eof())
                break;

	    c=get();
	    putback(c);
	    if(c=='*')
		skipComments();

            if (eof())
                break;

            readConstraint();
	    nbConstraintsRead++;
        }
	
	//Small check on the number of constraints
	if (nbConstraintsRead != nbConstr) {
	    throw new ParseException("Number of constraints read is different from metadata.");
	}
    }

    public static void main(String[] args) {
        try {
            if (args.length != 1)
                System.out.println("usage: SimpleParser <filename>");
            else {
                SimpleParser parser = new SimpleParser(args[0]);

                parser.parse();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}; // class SimpleParser


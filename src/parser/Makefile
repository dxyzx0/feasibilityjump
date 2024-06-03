TESTFILES=testlin.pb testnlc.pb mult.pb

all: runC runCC runJava

runCC: testCC
	@echo -e '\n#\n# Parser for the C++ language\n#\n'
	@for f in $(TESTFILES) ; do echo -e "\n#\n# $$f\n#\n" ; ./testCC $$f ; done

runC: testC
	@echo -e '\n#\n# Parser for the C language\n#\n'
	@for f in $(TESTFILES) ; do echo -e "\n#\n# $$f\n#\n" ; ./testC $$f; done

runJava: SimpleParser.class
	@echo -e '\n#\n# Parser for the Java language\n#\n'
	@for f in $(TESTFILES) ; do  echo -e "\n#\n# $$f\n#\n" ; java SimpleParser $$f ; done


testCC : SimpleParser.cc
	g++ -g -DuseGMP -o testCC SimpleParser.cc -lgmpxx -lgmp

testC : SimpleParser.c
	gcc -g -DuseGMP -o testC SimpleParser.c -lgmp

SimpleParser.class : SimpleParser.java
	javac SimpleParser.java

clean:
	rm -f testC testCC SimpleParser*.class

realclean: clean
	rm -f *~


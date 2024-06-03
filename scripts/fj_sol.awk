BEGIN {
    FS = " "
    OFS = ","
    min = 2e100
    nSol = 0
}

$2 == "(FJSOL)" {
    nSol++
    # print prev
    # print $0
    if (!min || $5 < min) {
        min = $5
        min_line = $0
        min_line_prev = prev
    }
}

{
    prev = $0
}

END {
    basename = FILENAME
    sub(".*/", "", basename)
    print "File processed: " basename

    # print "Line above with minimum fourth field: " min_line_prev
    sol = min_line_prev
    sub(".*:", "", sol)
    print basename, sol >> "bestSol.csv"
    
    # print "Line with minimum fourth field: " min_line
    # print "Min Obj: " bestSol[4]
    # print "Min Time (s): " bestSol[3]
    # print "Number of solutions: " nSol
    # if (FNR == NR) {
    #     print "filename,nSol,objBest,minTimeBest" > "results.csv"
    # }
    split(min_line, bestSol, " ")
    print basename, nSol, bestSol[5], bestSol[4] >> "results.csv"
}

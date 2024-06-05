BEGIN {
    FS = " "
    OFS = ","
    min = 2e100
    nSol = 0
}

($2 == "(FJSOL)") && ($4 == "time=") {
    nSol++
    # print prev
    # print $0
    if (!min || $5 < min) {  # this can overflow if $5 is too large
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
    sub(".*solution:", "", sol)
    print basename, sol >> "bestSol.csv"
    
    # print "Line with minimum fourth field: " min_line
    # print "Min Obj: " bestSol[7]
    # print "Min Time (s): " bestSol[5]
    # print "Number of solutions: " nSol
    # if (FNR == NR) {
    #     print "filename,nSol,objBest,minTimeBest" > "results.csv"
    # }
    split(min_line, bestSol, " ")
    print basename, nSol, bestSol[7], bestSol[5] >> "results.csv"
}

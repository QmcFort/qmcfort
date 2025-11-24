import os
import sys
from datetime import datetime

def read_report(fname):
    try:
        with open(fname) as file:
            report = [line.split() for line in file]
        del report[:2]
        return report
    except:
        print("File " + fname + " is not found!")
        #raise RuntimeError from None

def finish_report(report, fname):
    tend = datetime.now().strftime('%d/%m/%Y %H:%M:%S')
    tend = datetime.strptime(tend, '%d/%m/%Y %H:%M:%S')
    with open("timestamp", "r") as file:
        tstart = file.readline().strip()
        tstart = datetime.strptime(tstart, "%d/%m/%Y %H:%M:%S")
    os.remove("timestamp")


    l = "################################################################################"
    n = 0
    nsuccess = 0
    
    for item in report:
        n +=1
        if (item[2] == "True"):
            nsuccess += 1
            
    f = open(fname, "a")
    f.write(l + "\n")
    f.write("\n")
    f.write("Total number of tests    = " + str(n) + "\n")
    f.write("Number of passed tests   = " + str(nsuccess) + "\n")
    f.write("Number of failed tests   = " + str(n-nsuccess) + "\n")
    f.write("\n")
    f.write("\n")
    f.write("Tests started  at : " + str(tstart) + "\n")
    f.write("Tests finished at : " + str(tend)   + "\n")
    f.write("Time difference   : " + str(tend-tstart) + "\n")
    f.write(l + "\n")
    f.write(l + "\n")
    f.close()
    
    print("Integration test summary")
    print("========================")
    print("")
    print("    Total number of tests    = " + str(n))
    print("    Number of passed tests   = " + str(nsuccess))
    print("    Number of failed tests   = " + str(n-nsuccess))
    print("")
    print("")
    print("    Tests started  at : " + str(tstart))
    print("    Tests finished at : " + str(tend))
    print("    Time difference   : " + str(tend-tstart))
    print(l)
    print(l)
    return None

fname = sys.argv[1]
report = read_report(fname)
finish_report(report, fname)

import sys
import os

if len(sys.argv) >= 2:
    testname = sys.argv[1]
else:
    testname = "noname_test"

if len(sys.argv) == 3:
    tol = float(sys.argv[2])
else:
    tol = 1.0e-05

def read_values(fname):
    try:
        with open(fname) as file:
            lines = [line.rstrip() for line in file]
        values = []
        for item in lines:
            values.append(float(item))
        return values
    except:
        print("File " + fname + " is not found!")
        return []

def write_report(testname, success, maxdiff):
    fname = "../../report.log"
    str0  = "Counter"
    str0_ = "======="
    str1  = "Testname"
    str1_ = "========"
    str2  = "Success"
    str2_ = "======="
    str3  = "Maxdiff"
    str3_ = "======="
    l = "################################################################################"
    f = open(fname, "a")
    if (os.path.getsize(fname) == 0):
        f.write('%-10s %-40s %10s %20s \n' % (str0, str1, str2, str3))
        f.write('%-10s %-40s %10s %20s \n' % (str0_, str1_, str2_, str3_))
    index = max(sum(1 for line in open(fname))-1, 1)
    f.write('%-10d %-40s %10s %20s \n' % (index, testname, str(success), str("{:12.8f}".format(float(maxdiff)))))
    f.close()

    print(l)
    #print(l)
    print("Integration test " + testname)
    print("")
    if success:
        print("    Test passed!")
    else:
        print("    Test failed!")
    print("    Maximal difference is " + str("{:12.8f}".format(float(maxdiff))))
    print("")
    #print(l)
    print(l)
    return None

def test_arrays(a, b, tol=1.0e-06):
    success = True
    maxdiff = 0.0
    if (len(a)==0 or len(b)==0):
        success = False
        maxdiff = 1.0
    elif (len(a) != len(b)):
        success = False
        maxdiff = 1.0
    else:
        for i in range(len(a)):
            success = success and abs(a[i]-b[i])<tol
            maxdiff = max(maxdiff, abs(a[i]-b[i]))
    return success, maxdiff

a = read_values("values.txt")
b = read_values("ref_values.txt")

success, maxdiff = test_arrays(a, b, tol=tol)

write_report(testname, success, maxdiff)


import sys

origds = sys.argv[1]    # dataset (data + labels)
newds_x = origds+".X"   # data
newds_y = origds+".Y"   # labels

fp_origds  = open(origds, 'r')
fp_newds_x = open(newds_x, 'w')
fp_newds_y = open(newds_y, 'w')

counter = 1
empty_example = 0  # number of examples with all feature values being zero
for line in fp_origds.readlines():    
    tokens = line.strip().split(" ",1)
    if len(tokens) <= 1:
        print "line %d : %s" %(counter,line)
        empty_example = empty_example + 1
        continue
    fp_newds_y.write(tokens[0]+"\n")
    fp_newds_x.write(tokens[1]+"\n")
# this remove svmlight ids
#    for item in tokens[1:]:
#        if item.startswith("qid"):
#            continue
#        elif item.startswith("sid"):
#            continue
#        elif item.startswith("cost"):
#            continue
#        else:
#            fp_newds_x.write(item+" ")
#    fp_newds_x.write("\n")

    counter = counter + 1

print "empty example : ", empty_example
        
        
        



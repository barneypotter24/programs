#!/usr/bin/python
import sys
import timeit
from directed_hypergraph import *

def main(argv):
    # Specify the name of the file containing the nodes and edges
    title = argv[0]
    print("...\n...\n...")
#    print("...\n...\n...")
#    print("...\n...\n...")
#    dc = input("Use the D&C formulation? (y/n) ");
    nodeFile = "../hypergraphs/"+title+"-nodes.txt"
    edgeFile = "../hypergraphs/"+title+"-edges.txt"
    out = title+".lp"
    outDC = title+"_DC.lp"

    H = DirectedHypergraph()

    node_delimeter = ","
    column_delimeter = ";"

    # Read the edge and node files to create a weighted hypergraph.
    H.read(edgeFile, node_delimeter, column_delimeter)
    H.weight_nodes(nodeFile, node_delimeter, column_delimeter)
    weightNodes(H,title)

    # A = timeit.Timer(lambda: build_lp(H,out))
    # B = timeit.Timer(lambda: build_lp_dc(H,outDC))
    #
    # listA = A.repeat(repeat=100,number=1)
    #
    # listB = B.repeat(repeat=100,number=1)
    #
    # print("Average run time of build_lp:")
    # print(sum(listA)/len(listA))
    #
    # print("Average run time of build_lp_dc:")
    # print(sum(listB)/len(listB))

    build_lp(H,out)
    build_lp_dc(H,outDC)

def weightNodes(H,prefix):

    name2prize = nameToPrize(prefix)
    uniprot2name = uniprotToName(prefix)
    reactome2uniprot = reactomeIdToUniprot(prefix)

    count = 0

    for node in H.node_iterator():

        name = str(node)

        if not name in reactome2uniprot:
            continue

        UP = reactome2uniprot[name]

        if not UP in uniprot2name:
            continue

        up2n = uniprot2name[UP]

        if not up2n in name2prize:
            continue

        prize = float(name2prize[up2n])
        prize = prize * 10
        prize = str(prize)

        atts = H.get_node_attributes(node)
        pen = atts["penalty"]

        if name[0] == "P":

            attributes = {"prize" : prize, "penalty" : pen}
            H.add_node(name, attributes)
            count += 1

        elif name[0] == "C":

            #stuff
            proteins = proteinsInComplex(filename,name)
            prizes = []

            for prot in proteins:
                up = reactome2uniprot[prot]
                pz = float(name2prize[up])
                pz = pz * 10
                pz = str(pz)
                prizes.append(pz)

            minPrize = min(prizes)

            attributes = {"prize" : minPrize, "penalty" : pen}
            H.add_node(name, attributes)
            count += 1

    print("Nodes weighted:",count)

def nameToPrize(prefix):

    D = {}

    sep = ","

    filename = "outputFiles/"+prefix+"-name2prize.csv"

    in_file = open(filename, 'r')

    # Skip the header line
    in_file.readline()

    line_number = 2
    for line in in_file.readlines():
        line = line.strip()
        # Skip empty lines
        if not line:
            continue

        words = line.split(sep)
        if not (len(words) == 3):
            raise \
                IOError("Line {} ".format(line_number) +
                        "contains {} ".format(len(words)) +
                        "columns -- must contain only 3.")

        D[words[0]] = words[2]

        line_number += 1

    in_file.close()

    return D

def proteinsInComplex(prefix, nodeName):

    proteins = []

    sep = "\t"
    delim = ";"

    filename = "../../reactomeData/testOutput/"+prefix+"-complexes.txt"

    in_file = open(filename, 'r')

    # Skip the header line
    in_file.readline()

    line_number = 2
    for line in in_file.readlines():
        line = line.strip()
        # Skip empty lines
        if not line:
            continue

        words = line.split(sep)
        prots = words[3].split(delim)

        if words[0] == nodeName:

            proteins = prots


        line_number += 1

    in_file.close()

    return proteins

def reactomeIdToUniprot(prefix):

    D = {}

    filename = "../../reactomeData/testOutput/"+prefix+"-elements.txt"

    sep = "\t"
    delim = ";"
    delim2 = ":"

    in_file = open(filename, 'r')

    # Skip the header line
    in_file.readline()

    line_number = 2
    for line in in_file.readlines():
        line = line.strip()
        # Skip empty lines
        if not line:
            continue

        words = line.split(sep)
        if words[0][0] == "P":
            aliases = words[3].split(delim)
            #print(aliases)
            if len(aliases) > 1:
                UP = aliases[1].split(delim2)
            #print(UP)

            D[words[0]] = UP[1]

        line_number += 1

    in_file.close()

    return D

def uniprotToName(prefix):

    D = {}

    sep = "\t"

    filename = "../hypergraphs/"+prefix+"-conversion.txt"

    in_file = open(filename, 'r')

    # Skip the header line
    in_file.readline()

    line_number = 2
    for line in in_file.readlines():
        line = line.strip('\n')
        # Skip empty lines
        if not line:
            continue

        words = line.split(sep)
        D[words[1]] = words[0]

        line_number += 1

    in_file.close()

    #print("up2n: ",D)

    return D

def build_lp(Hypergraph,outputFile):

    """ Defines a function that will build a .lp file after being Given
    a hypergraph. The hypergraph must be weighted with
        - Edge weights (Added during H.read())
        - Node prized (Added during H.weight_nodes())
        - Node penalties (Added during H.weight_nodes())

    """

    # Open a new text file (default 'text.lp').
    # Edit this name to make a new file based on which hypergraph
    # is being used for the example.
    lp_file = open("outputFiles/"+outputFile, "w")

    # Write the first line of the .lp file.
    # This can be "Maximize\n" or "Minimize\n".
    lp_file.write("Maximize\n")
    # Begin the line for the objective function.
    lp_file.write(" obj: ")

    # Create three lists of the variables present in the linear program.
    xVariables = []
    # xVarVals = []
    dVariables = []
    # dVarVals = []
    eVariables = []
    # eVarVals = []

    # First, write sum_{v \in V}{g_{v}'a_{v}}

    isFirst = True
    for n in Hypergraph.node_iterator():
        if isFirst == True:

#            print(Hypergraph.get_node_attribute(n,"prize"))
            lp_file.write(str(Hypergraph.get_node_attribute(n,"prize")))
            lp_file.write(" ")
            lp_file.write(str(n)+"x")
            isFirst = False
        else:
            lp_file.write(" + ")

#            print(Hypergraph.get_node_attribute(n,"prize"))
            lp_file.write(str(Hypergraph.get_node_attribute(n,"prize")))
            lp_file.write(" ")
            lp_file.write(str(n)+"x")

        xVariables.append(str(n)+"x")

    # Next, write sum_{v \in V}{h_{v}'d_{v}}
    for n in Hypergraph.node_iterator():
        lp_file.write(" - ")
        if str(Hypergraph.get_node_attribute(n,"penalty")) == "inf":
            lp_file.write(str(0.0))
        else:
            lp_file.write(str(Hypergraph.get_node_attribute(n,"penalty")))
        lp_file.write(" ")
        lp_file.write(str(n)+"d")

        dVariables.append(str(n)+"d")

    # Last, write sum_{e \in E}{c_{e}'a_{e}}
    for e in Hypergraph.hyperedge_id_iterator():
        lp_file.write(" - ")
        lp_file.write(str(Hypergraph.get_hyperedge_weight(e)))
        lp_file.write(" ")
        lp_file.write(str(e))

        eVariables.append(str(e))

    xVariables.sort()
    dVariables.sort()
    eVariables.sort()

    # print("xVariables:",xVariables)
    # print("dVariables:",dVariables)
    # print("eVariables:",eVariables)

    # Done writing the objective function
    lp_file.write("\n")
    print("Objective function written successfully.")
    print("...\n...\n...")


    # Begin the linear constraints.
    lp_file.write("Subject To\n")


    # 4.3 and 4.4
    for e in Hypergraph.hyperedge_id_iterator():

        lp_file.write(" c44_"+str(e)+": ")

        isFirst = True

        for tail in Hypergraph.get_hyperedge_tail(e):

            if isFirst == True:
                lp_file.write(str(tail)+"x")
                isFirst = False
            else:
                lp_file.write(" + ")
                lp_file.write(str(tail)+"x")

        lp_file.write(" - "+str(len(Hypergraph.get_hyperedge_tail(e)))
                                        +" "+e+" >= 0\n")

        lp_file.write(" c43_"+str(e)+": ")

        isFirst = True

        for head in Hypergraph.get_hyperedge_head(e):

            if isFirst == True:
                lp_file.write(str(head)+"x")
                isFirst = False
            else:
                lp_file.write(" + ")
                lp_file.write(str(head)+"x")

        lp_file.write(" - "+str(len(Hypergraph.get_hyperedge_head(e)))
                                        +" "+e+" >= 0\n")

    # 4.5
    for i in range(len(xVariables)):
        lp_file.write(" c45_"+xVariables[i][0]+": ")
        lp_file.write(str(dVariables[i])+" - "+str(xVariables[i])+" <= 0\n")

    # 4.6
    for n in Hypergraph.node_iterator():
        lp_file.write(" c46_"+str(n)+": ")
        lp_file.write(str(n)+"d - "+str(n)+"x")

        for edge in Hypergraph.get_backward_star(n):
            lp_file.write(" + "+str(edge))

        lp_file.write(" >= 0\n")

    # Done writing the linear constraints.
    print("Linear constraints written successfully.")
    print("...\n...\n...")


    # Begin writing bounds.
    lp_file.write("Bounds\n")

    for n in Hypergraph.node_iterator():
        if str(Hypergraph.get_node_attribute(n,"penalty")) == "inf":
            lp_file.write(" "+str(n)+"x = 1\n")
            lp_file.write(" "+str(n)+"d = 0\n")

    #Done writing the bounds.
    print("Bounds written successfully.")
    print("...\n...\n...")


    # Write which binary variables are being optimized.
    lp_file.write("Binary\n")

    for x in xVariables:
        lp_file.write(" "+str(x)+"\n")

    for d in dVariables:
        lp_file.write(" "+str(d)+"\n")

    for e in eVariables:
        lp_file.write(" "+str(e)+"\n")

    # Done adding binary variables.

    # End the .lp file
    lp_file.write("End\n")
    print("Done.")
    lp_file.close()

def goesTo(Hypergraph,node,L,D,firstNode):

    for e in Hypergraph.get_forward_star(node):
        for h in Hypergraph.get_hyperedge_head(e):
            if h not in L:

                #L.append(h)
                L[h]=True
                name = str(firstNode)+'To'+str(h)
                D[name] = True
                goesTo(Hypergraph,h,L,D,firstNode)

def connectivityFinder(Hypergraph):

    # Make a dictionary that contains whether or not any 2 nodes are connected

    D = {}
    connected = set()

    count = 0

    for node in Hypergraph.node_iterator():
        for otherNode in Hypergraph.node_iterator():

            if count % 100000 == 0:
                print("count:",count)

            name = str(node)+'To'+str(otherNode)

            if node == otherNode:
                D[name] = True
            else:
                D[name] = False

            count += 1

    for n in Hypergraph.node_iterator():

        #connected = []
        connected = {}
        goesTo(Hypergraph,n,connected,D,n)

    return D

def build_lp_dc(Hypergraph,outputFile):

    """ Defines a function that will build a .lp file after being Given
    a hypergraph. The hypergraph must be weighted with
        - Edge weights (Added during H.read())
        - Node prized (Added during H.weight_nodes())
        - Node penalties (Added during H.weight_nodes())

    """

    D = connectivityFinder(Hypergraph)
    # print(D)

    ts = 0
    for key in D:
        if D[key] == True:
            ts += 1
    #print("ts: ",ts)


    # Open a new text file (default 'text.lp').
    # Edit this name to make a new file based on which hypergraph
    # is being used for the example.
    lp_file = open("outputFiles/"+outputFile, "w")

    # Write the first line of the .lp file.
    # This can be "Maximize\n" or "Minimize\n".
    lp_file.write("Maximize\n")
    # Begin the line for the objective function.
    lp_file.write(" obj: ")

    # Create three lists of the variables present in the linear program.
    xVariables = []
    # xVarVals = []
    dVariables = []
    # dVarVals = []
    eVariables = []
    # eVarVals = []
    cVariables = []


    # First, write sum_{v \in V}{g_{v}'a_{v}}

    isFirst = True
    for n in Hypergraph.node_iterator():
        if isFirst == True:
            lp_file.write(str(Hypergraph.get_node_attribute(n,"prize")))
            lp_file.write(" ")
            lp_file.write(str(n)+"x")
            isFirst = False
        else:
            lp_file.write(" + ")
            lp_file.write(str(Hypergraph.get_node_attribute(n,"prize")))
            lp_file.write(" ")
            lp_file.write(str(n)+"x")

        xVariables.append(str(n)+"x")

    # Next, write sum_{v \in V}{h_{v}'d_{v}}
    for n in Hypergraph.node_iterator():
        lp_file.write(" - ")
        if str(Hypergraph.get_node_attribute(n,"penalty")) == "inf":
            lp_file.write(str(0.0))
        else:
            lp_file.write(str(Hypergraph.get_node_attribute(n,"penalty")))
        lp_file.write(" ")
        lp_file.write(str(n)+"d")

        dVariables.append(str(n)+"d")

    # Dummy variables for chads
    for n in Hypergraph.node_iterator():
        lp_file.write(" + ")
        lp_file.write(str(0.0))
        lp_file.write(str(n)+"c")

        cVariables.append(str(n)+"c")

    # Dummy variables for connections
    for key in D:
        lp_file.write(" + 0.0 "+key)

    # Last, write sum_{e \in E}{c_{e}'a_{e}}
    for e in Hypergraph.hyperedge_id_iterator():
        lp_file.write(" - ")
        lp_file.write(str(Hypergraph.get_hyperedge_weight(e)))
        lp_file.write(" ")
        lp_file.write(str(e))

        eVariables.append(str(e))


    xVariables.sort()
    dVariables.sort()
    eVariables.sort()
    cVariables.sort()
    #print("proof that something is happening: ",cVariables)

    # print("xVariables:",xVariables)
    # print("dVariables:",dVariables)
    # print("eVariables:",eVariables)

    # Done writing the objective function
    lp_file.write("\n")
    print("Objective function written successfully.")
    print("...\n...\n...")


    # Begin the linear constraints.
    lp_file.write("Subject To\n")


    # 4.3 and 4.4
    for e in Hypergraph.hyperedge_id_iterator():

        lp_file.write(" c44_"+str(e)+": ")

        isFirst = True

        for tail in Hypergraph.get_hyperedge_tail(e):

            if isFirst == True:
                lp_file.write(str(tail)+"x")
                isFirst = False
            else:
                lp_file.write(" + ")
                lp_file.write(str(tail)+"x")

        lp_file.write(" - "+str(len(Hypergraph.get_hyperedge_tail(e)))
                                        +" "+e+" >= 0\n")

        lp_file.write(" c43_"+str(e)+": ")

        isFirst = True

        for head in Hypergraph.get_hyperedge_head(e):

            if isFirst == True:
                lp_file.write(str(head)+"x")
                isFirst = False
            else:
                lp_file.write(" + ")
                lp_file.write(str(head)+"x")

        lp_file.write(" - "+str(len(Hypergraph.get_hyperedge_head(e)))
                                        +" "+e+" >= 0\n")

    # 4.5
    for i in range(len(xVariables)):
        lp_file.write(" c45_"+xVariables[i][0]+": ")
        lp_file.write(str(dVariables[i])+" - "+str(xVariables[i])+" <= 0\n")

    # 4.6
    for n in Hypergraph.node_iterator():
        lp_file.write(" c46_"+str(n)+": ")
        lp_file.write(str(n)+"d - "+str(n)+"x")

        for edge in Hypergraph.get_backward_star(n):
            lp_file.write(" + "+str(edge))

        lp_file.write(" >= 0\n")

    # Danglers and Chads issue

    # D&C1
    for i in range(len(xVariables)):
        lp_file.write(" dc1_"+xVariables[i][0]+": ")
        lp_file.write(str(cVariables[i])+" - "+str(xVariables[i])+" <= 0\n")

    # D&C2
    for n in Hypergraph.node_iterator():
        lp_file.write(" dc2_"+str(n)+": ")
        lp_file.write(str(n)+"c - "+str(n)+"x")

        for edge in Hypergraph.get_forward_star(n):
            lp_file.write(" + "+str(edge))

        lp_file.write(" >= 0\n")

    # D&C3

    for n in Hypergraph.node_iterator():
        for m in Hypergraph.node_iterator():
            name = str(n)+'To'+str(m)

            lp_file.write(" dc3_"+name+": ")
            lp_file.write(str(n)+"c + "+str(m)+"d + "+name+" <= "+str(2)+"\n")


    # Done writing the linear constraints.
    print("Linear constraints written successfully.")
    print("...\n...\n...")


    # Begin writing bounds.
    lp_file.write("Bounds\n")

    for n in Hypergraph.node_iterator():
        if str(Hypergraph.get_node_attribute(n,"penalty")) == "inf":
            lp_file.write(" "+str(n)+"x = 1\n")
            lp_file.write(" "+str(n)+"d = 0\n")

    # C&D bounds

    for key in D:
        if D[key] == True:
            lp_file.write(" "+key+" = 1\n")
        else:
            lp_file.write(" "+key+" = 0\n")

    #Done writing the bounds.
    print("Bounds written successfully.")
    print("...\n...\n...")


    # Write which binary variables are being optimized.
    lp_file.write("Binary\n")

    for x in xVariables:
        lp_file.write(" "+str(x)+"\n")

    for d in dVariables:
        lp_file.write(" "+str(d)+"\n")

    for e in eVariables:
        lp_file.write(" "+str(e)+"\n")

    for key in D:
        lp_file.write(" "+key+"\n")

    # Done adding binary variables.

    # End the .lp file
    lp_file.write("End\n")
    print("Done.")
    lp_file.close()

def hypershrub(Hypergraph,outputFile):

    """ Defines a function that will build a .lp file after being Given
    a hypergraph. The hypergraph must be weighted with
        - Edge weights (Added during H.read())
        - Node prized (Added during H.weight_nodes())
        - Node penalties (Added during H.weight_nodes())

    """

    #D = connectivityFinder(Hypergraph)
    # print(D)

    #ts = 0
    #for key in D:
    #    if D[key] == True:
    #        ts += 1
    #print("ts: ",ts)


    # Open a new text file (default 'text.lp').
    # Edit this name to make a new file based on which hypergraph
    # is being used for the example.
    lp_file = open("outputFiles/"+outputFile, "w")

    ################################################################
    # Objective Function #
    ######################

    ################################

    # Write the first line of the .lp file.
    # This can be "Maximize\n" or "Minimize\n".
    lp_file.write("Maximize\n")
    # Begin the line for the objective function.
    lp_file.write(" obj: ")

    # Create three lists of the variables present in the linear program.
    xVariables = []
    # xVarVals = []
    dVariables = []
    # dVarVals = []
    eVariables = []
    # eVarVals = []
    # cVariables = []

    ################################
    # First, write sum_{v \in V}{p_{v}\alpha_{v}}

    isFirst = True
    for n in Hypergraph.node_iterator():
        if isFirst == True:
            lp_file.write(str(Hypergraph.get_node_attribute(n,"prize")))
            lp_file.write(" ")
            lp_file.write("alpha_"+str(n))
            isFirst = False
        else:
            lp_file.write(" + ")
            lp_file.write(str(Hypergraph.get_node_attribute(n,"prize")))
            lp_file.write(" ")
            lp_file.write("alpha_"+str(n))

        xVariables.append("alpha_"+str(n))

    ################################
    # Next, write sum_{v \in V}{r_{v}\rho_{v}}
    for n in Hypergraph.node_iterator():
        lp_file.write(" - ")
        if str(Hypergraph.get_node_attribute(n,"penalty")) == "inf":
            lp_file.write(str(0.0))
        else:
            lp_file.write(str(Hypergraph.get_node_attribute(n,"penalty")))
        lp_file.write(" ")
        lp_file.write("rho_"+str(n))

        dVariables.append("rho_"+str(n))

    ################################
    # Write sum_{e \in E}{c_{e}\alpha_{e}}
    for e in Hypergraph.hyperedge_id_iterator():
        lp_file.write(" - ")
        lp_file.write(str(Hypergraph.get_hyperedge_weight(e)))
        lp_file.write(" ")
        lp_file.write("alpha_"+str(e))

        eVariables.append("alpha_"+str(e))

    ################################
    # Dummys #
    ##########

    ################################
    # Dummy variable for sigma:

    sigmas = []
    for n in Hypergraph.node_iterator():
        lp_file.write(" + 0.0 sigma_"+str(n))

        sigmas.append("sigma_"+str(n))

    ################################
    # Dummy variable for pi:

    pis = []
    for n in Hypergraph.node_iterator():
        for e in Hypergraph.get_backward_star(n):
            lp_file.write(" + 0.0 pi_"+str(e)+","+str(v))

            pis.append("pi_"+str(e)+","+str(v))

    ################################
    # Dummy variable for beta:

    betas = []
    for n in Hypergraph.node_iterator():
        lp_file.write(" + 0.0 beta_"+str(n))

        betas.append("beta_"+str(n))

    ################################
    # Sorting #
    ###########

    xVariables.sort()
    dVariables.sort()
    eVariables.sort()
    sigmas.sort()
    pis.sort()
    betas.sort()
    #cVariables.sort()
    #print("proof that something is happening: ",cVariables)

    # print("xVariables:",xVariables)
    # print("dVariables:",dVariables)
    # print("eVariables:",eVariables)

    ################################

    # Done writing the objective function
    lp_file.write("\n")
    print("Objective function written successfully.")
    print("...\n...\n...")

    ################################################################
    # Linear Constraints #
    ######################

    # Begin the linear constraints.
    lp_file.write("Subject To\n")

    ################################
    # (12)

    for v in Hypergraph.node_iterator():

        lp_file.write(" 12_"+str(v)+": ")

        lp_file.write(" sigma_"+str(v)+" = ")

        first = True
        for e in Hypergraph.get_backward_star(v):
            if first == True:
                lp_file.write("pi_"+str(e)+","+str(v)+"\n")
                first = False
            else:
                lp_file.write(" + pi_"+str(e)+","+str(v)+"\n")

    ################################
    # (13)

    #
    # In "bounds"
    #

    ################################
    # (14)

    for v in Hypergraph.node_iterator():
        lp_file.write(" 14_"+str(v)+": ")
        lp_file.write("rho_"+str(v)+" <= alpha_"+str(v)+"\n")

    ################################
    # (15)

    for v in Hypergraph.node_iterator():
        lp_file.write(" 15_"+str(v)+": ")
        lp_file.write("sigma_"+str(v)+" <= 1 - rho_"+str(v)+"\n")

    ################################
    # (16)

    for v in Hypergraph.node_iterator():
        lp_file.write(" 16_"+str(v)+": ")
        lp_file.write("sigma_"+str(v)+" <= alpha_"+str(v)+"\n")

    ################################
    # (17)

    for e in Hypergraph.hyperedge_id_iterator():
        for v in Hypergraph.get_hyperedge_tail(e):
            lp_file.write(" 17_"+str(e)+","+str(v)+": ")
            lp_file.write("alpha_"+str(e)+" <= alpha_"+str(v)+"\n")

    ################################
    # (18)

    for e in Hypergraph.hyperedge_id_iterator():
        lp_file.write(" 18_"+str(e)+": ")
        lp_file.write("alpha_"+str(e)+" <= ")

        first = True
        for v in Hypergraph.get_hyperedge_head(e):
            if first == True:
                lp_file.write("pi_"+str(e)+","+str(v))
                first = False
            else:
                lp_file.write(" + pi_"+str(e)+","+str(v))

        lp_file.write("\n")

    ################################
    # (19)

    for e in Hypergraph.hyperedge_id_iterator():
        lp_file.write(" 19_"+str(e)+": ")
        lp_file.write("alpha_"+str(e)+" >= ")

        first = True
        for v in Hypergraph.get_hyperedge_head(e):
            if first == True:
                lp_file.write("pi_"+str(e)+","+str(v))
                first = False
            else:
                lp_file.write(" + pi_"+str(e)+","+str(v))

        lp_file.write("\n")

    #########################
    #                       #
    #  Talk to Jim or Anna  #
    #                       #
    #########################

    ################################
    # (20)

    #
    # Doesn't require a constraint.
    #

    ################################
    # (21)

    #
    # In "bounds"
    #

    ################################
    # (22)

    epsilon = 1/len(xVariables)

    for e in Hypergraph.hyperedge_id_iterator():
        for u in Hypergraph.get_hyperedge_tail(e):
            for v in Hypergraph.get_hyperedge_tail(e):
                lp_file.write(" 22_"+str(e)+","+str(u)+","+str(v)+": ")
                lp_file.write("1 - pi_"+str(e)+","+str(v)+" + beta_"+str(v)+" >= beta_"+str(u)+" + "+str(epsilon)+"\n")

    ################################

    # Done writing the linear constraints.
    print("Linear constraints written successfully.")
    print("...\n...\n...")

    ################################################################
    # Bounds #
    ##########

    # Begin writing bounds.
    lp_file.write("Bounds\n")

    ################################
    # Any node designated as a target is automatically included in the solution
    # Nodes are designated as targets by having their rho value set to inf

    for n in Hypergraph.node_iterator():
        if str(Hypergraph.get_node_attribute(n,"penalty")) == "inf":
            lp_file.write(" "+str(n)+"x = 1\n")
            lp_file.write(" "+str(n)+"d = 0\n")

    ################################
    # Sigma constraints

    for v in Hypergraph.node_iterator():
        lp_file.write(" 0 <= sigma_"+str(v)+" <= 1\n")

    ################################
    # Beta constraints

    ffor v in Hypergraph.node_iterator():
        lp_file.write(" 0 <= beta_"+str(v)+" <= 1\n")

    #Done writing the bounds.
    print("Bounds written successfully.")
    print("...\n...\n...")

    ################################################################
    # Binary #
    ##########

    # Write which binary variables are being optimized.
    lp_file.write("Binary\n")

    for x in xVariables:
        lp_file.write(" "+str(x)+"\n")

    for d in dVariables:
        lp_file.write(" "+str(d)+"\n")

    for e in eVariables:
        lp_file.write(" "+str(e)+"\n")

    for s in sigmas:
        lp_file.write(" "+str(s)+"\n")

    for p in pis:
        lp_file.write(" "+str(p)+"\n")

    # Done adding binary variables.

    ################################################################
    # Rational #
    ############

    # Write which rational variables are being optimized.
    # lp_file.write("Rational\n")
    #
    # for b in betas:
    #     lp_file.write(" "+str(b)+"\n")

    ################################################################
    # End #
    #######

    # End the .lp file
    lp_file.write("End\n")
    print("Done.")
    lp_file.close()

    ################################################################

if __name__ == "__main__":
    main(sys.argv[1:])

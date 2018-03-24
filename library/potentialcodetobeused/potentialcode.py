def structureRetriever(structure_id,filename):
    # this function returns to the caller a structure object from the pdb file name by parsing the pdb file and returns a structure object
    parser = PDBParser()
    structure = parser.get_structure(structure_id,filename)
    return structure

def structureParser(structure_id,filename):
    # this function returns to the caller a structure object from the pdb file name by parsing the pdb file
    parser = PDBParser()
    structure = parser.get_structure(structure_id,filename)
    # The Structure object follows the so-called SMCRA (Structure/Model/Chain/Residue/Atom)
    for model in structure:
        print model
        for chain in model:
            print chain
            for residue in chain:
                print residue
                for atom in residue:
                    print atom

def neighborsingleSearcher(structure,distance):
    atom_list = Selection.unfold_entities(structure, 'A')
    ns = NeighborSearch(atom_list)
    # print dir(atom_list[0]) #print all possible attributes or methods of the object
    # print atom_list[0].__dict__ #print all info of this object
    print len(atom_list)
    center = atom_list[0].get_coord()
    neighbors = ns.search(center,distance)
    residue_list = Selection.unfold_entities(neighbors, 'R')
    print residue_list
    print "There are %s of nearest residues", len(residue_list)
    for i in residue_list:
        print dir(i)
        print i.get_resname()
        print i.child_list
        # the next thing we can do is get the individual atoms in the residue, if it's a hydrogen atom then it's hydrogen bonding
        # we can start classifyign the contacts here if this is then incrementer counter to. If this is then increment the counter here.

def neighborSearcher(structure,distance):
    atom_list = Selection.unfold_entities(structure, 'A')
    ns = NeighborSearch(atom_list)
    # print dir(atom_list[0]) #print all possible attributes or methods of the object
    # print atom_list[0].__dict__ #print all info of this object
    indexing = [0,0,0,0]
    # for atom in atom_list:
    for atom in atom_list[0:100]:
        center = atom.get_coord()
        neighbors = ns.search(center,distance)
        residue_list = Selection.unfold_entities(neighbors, 'R')
        print atom
        print residue_list
        print "There are %s nearest residues." % len(residue_list)
        for residue in residue_list:
            # index = 0
            # indexing = [0,0,0,0]
            print residue.get_resname()
            # if residue.get_resname().strip() == "  G":
            # we have to clean this string up using .strip
            for atom in residue:
                print atom.get_name()
                if residue.get_resname().strip() == "G" and atom.get_name() == "H1":
                #print atom
                  print "matching G"
                #   index += 1
                  indexing[0] +=1
                elif residue.get_resname() == "C" and atom.get_name() == "H1":
                  print "matching C"
                  indexing[1] +=1
                elif residue.get_resname() == "A" and atom.get_name() == "H1":
                  print "matching A"
                  indexing[2] +=1
                elif residue.get_resname() == "U" and atom.get_name() == "H1":
                  print "matching U"
                  indexing[3] +=1
            # print "There are %s G to H1, %s C to H1, %s A to H1, and %s U to H1 interactions." %(indexing[0],indexing[1],indexing[2],indexing[3])
    print "There are %s G to H1, %s C to H1, %s A to H1, and %s U to H1 interactions." %(indexing[0],indexing[1],indexing[2],indexing[3])
            # print residue.child_list
        # the next thing we can do is get the individual atoms in the residue, if it's a hydrogen atom then it's hydrogen bonding
        # we can start classifyign the contacts here if this is then incrementer counter to. If this is then increment the counter here.
# run hbplus in hydrogenbond mode


def aminoacidrnamatcher(aminoacid,nucleotidebase):
    aminoacidlist=["ARG","ALA","ARG","GLY","CYS,""ILE","LYS","MET","PHE,""PRO","SER","THR","TYR","VAL"]
    nucleotidebase=["U","G","C","A"]
    for aa in aminoacidlist:
        for nb in nucleotidebase:
            print aa
            print nb

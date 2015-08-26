from Bio import pairwise2
from Bio import SeqIO
from py2neo import Path, neo4j
from py2neo import authenticate, Graph, Node, Relationship

# set up authentication parameters
authenticate("localhost:7474", "neo4j", "neo4j")

# connect to authenticated graph database
sgraph = Graph()

seqs = []
scores = []


#tx = sgraph.cypher.begin()


list_of_nodes = dict()
handle = open("PF16840_full.txt", "rU")

for record in SeqIO.parse(handle, "fasta") :
	seqs.append(record)
handle.close()

'''
This one can be optimized by calculating triangle 
of the matrix and then copy this to the second triangle
'''

for i in seqs:
	for j in seqs:
		if i != j:
			alignment = pairwise2.align.globalxx(i.seq.tostring(), j.seq.tostring())
			scores.append((i.id, j.id, alignment[0][4]))

'''
Let's create sample database with sequences and scores
'''

for seq in seqs:
	list_of_nodes[seq.id] = sgraph.create({"name": seq.id})[0]


for record in scores:
	rel = sgraph.create((list_of_nodes[record[0]], record[2], list_of_nodes[record[1]]))


print list_of_nodes["U3CDF9_CALJA/1-64"]

'''for pair in scores:
	node1 = Node("Protein", name= pair[0])
	node2 = Node("Protein", name= pair[1])
	distance = Relationship(node1, pair[2], node2)
	graph.create(distance)
'''

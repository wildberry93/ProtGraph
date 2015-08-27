from Bio import pairwise2
from Bio import SeqIO
from py2neo import Path, neo4j
from py2neo import authenticate, Graph, Node, Relationship

# set up authentication parameters
authenticate("localhost:7474", "neo4j", "3bitbit")

# connect to authenticated graph database
sgraph = Graph()

seqs = []
scores = []
THRESHOLD = 70

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
	list_of_nodes[seq.id].add_labels("Protein")

'''
Remove "backwards" relations and remove pairs with similarity lower than 70
'''
uniq_scores = {d[:2]:d for d in scores if d[2] > THRESHOLD}
uniq_scores = uniq_scores.values()

print len(uniq_scores)

for record in uniq_scores:
	rel = sgraph.create((list_of_nodes[record[0]], record[2], list_of_nodes[record[1]]))

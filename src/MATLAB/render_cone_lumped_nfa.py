from pythomata import SimpleDFA
from pythomata.impl.simple import SimpleNFA

characters = {'a':1,'b':2,'c':3,'d':4,'e':5,'f':6,'g':7,'h':8,'i':9,'j':10,'k':11,'l':12,'m':13,\
    'n':14,'o':15,'p':16,'q':17,'r':18,'s':19,'t':20,'u':21,'v':22,'w':23,'x':24,'y':25,'z':26}

from sys import argv
''' arg1: NFA / DFA file name '''

file = open(argv[1], 'r')
lines = file.readlines()
file.close()

states = int(lines[0].split()[0])

accepting_states = [ i+1 for i in list(map(int, lines[2].split()[1::])) ]

initial_state = int(lines[3].split()[0])+1

partition = list(range(states))

print(lines[1])
alphabet = len(lines[1])-1
print("size alphabet: ",alphabet)

det_trans_table = [[None]*states for i in range(alphabet)]
for i in lines[4::]:
    head, char, tail = i.split()
    head=int(head)+1
    tail=int(tail)+1
    char=characters[char]-1
    det_trans_table[char][head-1]=tail

num_states = states
states = [str(i) for i in range(1,states + 1)]
initial_state = str(initial_state)
accepting_states = {str(i) for i in accepting_states }
alphabet = list(characters.keys())[:alphabet]


trans_dict = {i:{j:0 for j in alphabet} for i in states}

for j in range(len(det_trans_table)):
    for k in range(len(det_trans_table[j])):
        if det_trans_table[j][k] != None:
            trans_dict[states[k]][alphabet[j]] = str(det_trans_table[j][k])

jpop = [[] for i in states]
for i in states:
    jflag = False
    for j in alphabet:
        if trans_dict[i][j] == 0:
            trans_dict[i].pop(j)
#     if len(trans_dict[i]) 
# for l in range(len(jpop)):
#     for j in jpop[l]:
#         trans_dict[str(l+1)].pop(j)
# ipop = []
# for l in trans_dict.keys():
#     if len(trans_dict[l])==0:
#         ipop.append(l)
# for i in ipop:
#     trans_dict.pop(i)

#print(trans_dict)

alphabet = set(alphabet)
states = set(states)

dfa = SimpleDFA(
    states,
    alphabet,
    initial_state,
    accepting_states,
    trans_dict
)

#dfa = dfa.minimize().trim()
#print("states: ",len(dfa.states))
#dfa = dfa.complete()
#print("DFA initial: ",dfa.initial_state)

digraph = dfa.to_graphviz()

digraph.render(argv[1]+'_image')
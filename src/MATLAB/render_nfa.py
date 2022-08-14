from sys import argv

from finite_automata import NFA#, DFA

#import pythomata
#from pythomata.alphabets import MapAlphabet
from pythomata.impl.simple import SimpleNFA

# read input from NFA description file
file = open(argv[1], 'r')
lines = file.readlines()
file.close()

nfa = NFA()

nfa.construct_nfa_from_file(lines)
# print(nfa.num_states)
print(nfa.states)
# print(nfa.symbols)
# print(nfa.num_accepting_states)
# print(nfa.accepting_states)
# print(nfa.start_state)
# print(nfa.transition_functions)

nfa.symbols = [x for x in nfa.symbols if x != ' ']
trans_dict = dict()
for i in nfa.states:
    trans_dict[str(i)] = {}
    for a in nfa.symbols:
        for x,y,z in nfa.transition_functions:
            if i == x and a == y:
                if str(a) not in trans_dict[str(i)]:
                    trans_dict[str(i)][str(a)]=[str(z)]
                else:
                    trans_dict[str(i)][str(a)]+=[str(z)]

inds = [x for x in trans_dict if not trans_dict[x]]
for x in inds:
    trans_dict.pop(x)

for i in trans_dict:
    for j in trans_dict[i]:
        trans_dict[i][j] = set(trans_dict[i][j])

nfa = SimpleNFA(
            {str(i) for i in nfa.states},
            set(nfa.symbols),
            str(nfa.start_state),
            {str(i) for i in nfa.accepting_states},
            trans_dict,
        )

digraph = nfa.to_graphviz()
digraph.render(argv[1]+"_render_nfa")

# dfa = nfa.determinize().minimize().trim()
# digraph = dfa.to_graphviz()

# digraph.render(argv[1]+"_render_dfa")
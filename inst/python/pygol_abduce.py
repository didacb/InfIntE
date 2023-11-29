
from pygol import *

"""
def dany():
	with open('BK.pl') as f: data = [line.strip() for line in f.read().splitlines() if line]
	const=["yes","no","up","down"]
	constant_details={"species":[1]}
	positive_example_list=read_example("pos_example.f")
	const=read_constants_bk(file=data, constants=const, relation=constant_details)
	P, N = bottom_clause_generation(file=data, constant_set = const,  container = "memory",positive_example=positive_example_list, negative_example=[])
	H=["abundance(C1,C2,S1,up):-presence(C2,S2,yes)&presence1(C1,S2,no)&effect_up(S2,S1)","abundance(C1,C2,S1,down):-presence(C2,S2,yes)&presence1(C1,S2,no)&effect_down(S2,S1)"]
	abduce=['effect_up','effect_down' ]
	#random_positive_example= random.choices(positive_example_list, k=500)
	coverage=pygolm_abduction(P, abduce,  positive_example_list=random_positive_example, constant_set=const, meta_rule=H, metric="predictive_power")

"""

def generate_bottom_clause(file, constant_set, positive_example, negative_example, container="dict"):
	P, N = bottom_clause_generation( negative_example=negative_example, file=file, constant_set =constant_set , positive_example=positive_example, container = container, depth=2)
	return P



def abduction(P, abduce, positive_example_list, constant_set, meta_rule, metric):
	coverage=pygol_abduction(P, abduce,  positive_example_list=positive_example_list, constant_set=constant_set, meta_rule=meta_rule, metric="predictive_power")
	return coverage

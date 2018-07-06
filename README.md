# MFS-greedy
Supporting docs for the "Advances in Greedy Optimization Algorithm of MFS Nodes and Collocation Points" paper

singleImage: Main file. Generates an image of the bunny based on # of nodes, whether approximation or error was chosen to be displayed, can display collocations and sources, some misc options

bunny_shift_5_5_5: (x,y,z) data for the bunny. xls format

bunnyGenApprox: generates approximation function on the surface of the bunny

bunnyGenError: generates error function on the surface of the bunny

errorFn: used by bunnyGenApprox and bunnyGenError to determine best collocation-source pairs

nodePoolGenGrid: generates pool of nodes around the bunny in a grid

nodePoolGenRand: generates pool of nodes around the bunny, uniformly distributed. Intended for future research

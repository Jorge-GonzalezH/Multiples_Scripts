from modeller import *
#Este script trabaja para generar modelos por homologia de la proteina metamorfica MAD2
env = environ()
aln = alignment(env)
for (pdb, chain) in (('4aez', 'B'), ('1go4', 'A'), ('2vfx', 'A'), 
                    ('2v64', 'D'), ('2qyf', 'C')):
    m = model(env, file=pdb, model_segment=('FIRST:'+chain, 'LAST:'+chain))
    aln.append_model(m, atom_files=pdb, align_codes=pdb+chain)
aln.malign()
aln.malign3d()
aln.compare_structures()
aln.id_table(matrix_file='family.mat')
env.dendrogram(matrix_file='family.mat', cluster_cut=-1.0)

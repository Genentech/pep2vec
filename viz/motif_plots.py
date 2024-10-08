import logomaker as lm
import numpy as np

def motif_plot_depletion(axes, seq_array, special_chars=True):
    AA_dict = {'A': 0.06871025,
     'C': 0.02257031,
     'D': 0.04791498,
     'E': 0.0706464,
     'F': 0.03679829,
     'G': 0.06437484,
     'H': 0.02621163,
     'I': 0.04396082,
     'K': 0.05717299,
     'L': 0.09942579,
     'M': 0.02212694,
     'N': 0.03604085,
     'P': 0.06247835,
     'Q': 0.0478216,
     'R': 0.05615158,
     'S': 0.084151,
     'T': 0.05399473,
     'V': 0.0605339,
     'W': 0.01199009,
     'Y': 0.02692467,
     '*': 0.0001,
     '$': 0.0001,
     '[':.0001,
     ']':.0001,
     '#':.0001,
     'U':.0001,
     'X':.0001,}
    
    counts_mat = lm.alignment_to_matrix(seq_array)
    probs = lm.transform_matrix(counts_mat, from_type='counts', to_type='probability')
    probs_start = probs.copy()

    if special_chars == False:
        probs[['*', '$', '[', ']', 'U', 'X','#']] = 0.0
        probs = probs.drop(columns=['*', '$', '[', ']', 'U', 'X','#'])
        probs[:] = (1 / probs.sum(axis=1)).values.reshape(-1,1) * probs.values
    
    for col in probs:
        baseline = AA_dict[col]
        probs[col] = probs[col] * np.log2(probs[col]/baseline)

    if axes is not None:      lm.Logo(probs, ax=axes, color_scheme='chemistry')
    return probs, probs_start

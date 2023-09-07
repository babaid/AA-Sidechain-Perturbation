# AA-Sidechain-Perturbation
This is a proof of concept python project, where the aim is to perturb sidechains of amino acids by rotating them around predefined axes. This can be of use for GNNs.


The goal is to perform random rotations of amino acid sidechains as described in [Pre-training of Graph Neural Network for Modeling Effects of Mutations on Protein-Protein Binding Affinity](https://arxiv.org/abs/2008.12473).

For now I want not to sample the distributions descirbed in the [Smooth Backbone-Dependent Rotamer Library 2010](http://dunbrack.fccc.edu/lab/bbdep2010), but also allow (probably) unrealistic states, where the only condition is that atoms do not clash due to the rotation.

In this stage, the python package is just a prototype, highly unoptimized. 
I will create an optimized version in C++/CUDA so that it can be used to create large datasets on clusters in a short amount of time.


# Requirements
The only requirements are
- rdkit
- numpy
- scipy
- loguru

Start by installing rdkit in a conda/mamba environment, and everything else should be straightforward.
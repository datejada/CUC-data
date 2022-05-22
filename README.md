# CUC-data
This repository has two case studies for Clustered Unit Commitment problems:
1. Modified version of the New England IEEE 39-bus System
2. Modified version of the IEEE 118-bus system

The optimization model description is available in:

G. Morales-España and D. A. Tejada-Arango, _"Modeling the Hidden Flexibility of Clustered Unit Commitment,"_ in IEEE Transactions on Power Systems, vol. 34, no. 4, pp. 3294-3296, July 2019.
doi: 10.1109/TPWRS.2019.2908051

And as a preprint in:

G. Morales-Espana and D. A. Tejada-Arango, _“Modelling the Hidden Flexibility of Clustered Unit Commitment,”_ arXiv e-prints, p. [arXiv:1811.02622](https://arxiv.org/abs/1811.02622), Nov. 2018.

The file __'StarGenLite_Clustered-UC.gms'__ has the GAMS code with the implementation of the optimization models. The input data is defined in GAMS using the keyword user1 in the GAMS' command line. For example: _user1=CUC_IEEE118-busSystem_

# Active-Polymer-Conformation
This project allows studying the conformation of active polymers that are driven by correlated random excitations.
It is based on analytic formulae that relate 

1. the excitation correlation function, and 
2. the mean squared separation between different monomers of this polymer.

Using these relations, one can transform back and forth between excitation correlation matrix and polymer conformation. 
Then, for any given polymer configuration, one can extract an apparent exctitation correlation matrix that leads to such a conformation.

We here consider discrete polymers whose response matrix, in Fourier space, is diagonal.
This allows us to exploit several simplifications, by transforming to Fourier space and back.

This code accompanies the following publication:

**Polymer folding through active processes recreates features of genome organization**, Andriy Goychuk<sup>\*</sup>, Deepti Kannan<sup>\*</sup>, Arup K. Chakrabory<sup>+</sup>, and Mehran Kardar<sup>+</sup>, [bioRxiv (2022)](https://doi.org/10.1101/2022.12.24.521789)

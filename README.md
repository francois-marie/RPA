# RANDOM PHASE APPROXIMATION (RPA) FOR THE TWO-DIMENSIONAL HUBBARD MODEL


This research project focuses on the investigation of magnetic fluctuations in one of the most fundamental models for correlated electron systems (like the cuprates), the two-dimensional Hubbard model. By applying a quantum many-body technique, the random phase approximation (RPA), the magnetic fluctuations of the system are analysed by means of its spin-spin susceptibility and magnetic correlation length for specific parameter regimes.

## Configuration

A configuration is defined thanks to the following parameters:
* n_k, defines the mesh in the k space
* n_q, same for a vector q in k-space
* t,  hopping parameter
* beta,  inverse temperature
* mu, chemical potential
* U, interaction parameter

From this set of parameters we can define a configuration, for example :

```
config = Configuration(40, 40, 1., 15., 0., 2.)
```

which will instanciate the values for the energy, the fermi occupation, the suceptibility in the tight binding approximation and the susceptibility in the RPA.

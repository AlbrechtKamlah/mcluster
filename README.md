# mcluster-LevelC

This code is the version of mcluster presented in the paper (appendix B)

Kamlah et al. (2022a) - https://ui.adsabs.harvard.edu/abs/2022MNRAS.511.4060K/abstract 

and it is forked from the mcluster version ([see here](https://github.com/agostinolev/mcluster)) presented in 

Leveque et al. (2022) - https://ui.adsabs.harvard.edu/abs/2022MNRAS.514.5739L/abstract

It contains the so-called _level C_ stellar evolution, which is described in Kamlah et al. (2022a).
This code was used in setting up initial conditions in, e.g

Kamlah et al. (2022b) - https://ui.adsabs.harvard.edu/abs/2022MNRAS.516.3266K/abstract

## mcluster - a tool to make a star cluster

The original mcluster software ([see here](https://github.com/ahwkuepper/mcluster)) is an open source code, which is
used to either set up initial conditions for _N_-body computations or
to generate artificial star clusters for direct investigation.

Kuepper et al. (2011) - https://ui.adsabs.harvard.edu/abs/2011MNRAS.417.2300K/abstract


## level C stellar evolution 

We call level C stellar evolution in reference to the paper Kamlah et al. (2022a), 
where an implementation of new stellar evolution into Nbody6++GPU ([see here for the Heidelberg/Beijing version of the code](https://github.com/kaiwu-astro/Nbody6PPGPU-beijing)), 
MOCCA and mcluster are discussed in detail. By our definition you also have Levels A, B and an upcoming level D stellar evolution in planning. 
We define all the levels as follows (see appendix A in Kamlah et al. 2022a):

* Level A: stellar evolution settings that mirror in part the settings in the Dragon simulations of GCs (see )
* Level B: stellar evolution settings that have been tested extensively and may be used without concern. A selection of these
should be enabled in the next gravitational million-body simulations.
* Level C: stellar evolution settings that are available in the codes, but those that are not present in level B have not yet
undergone sufficient testing and are therefore deemed experimental as of the writing of this paper.
* Level D: stellar evolution settings that will be added in the next iteration of stellar evolution updates, see also section 5.2 for
details on these.
* ...

In the more distant future, we will sequentially add new levels (the next one would be Level E), where we group further planned stellar
evolution updates on top of the preceding level (in this case level D) in Nbody6++GPU, MOCCA & McLuster together. We hope that
this will greatly help in the documentation and aid the future user of the codes to properly choose SSE & BSE settings in his or her
simulations.


_**This README is incomplete and will be updated over the next weeks**_

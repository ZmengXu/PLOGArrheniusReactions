[![OpenFOAM version](https://img.shields.io/badge/OpenFOAM-8-brightgreen)](https://github.com/OpenFOAM/OpenFOAM-8)
# PLOGArrheniusReactions
A library for OpenFOAM to handel the `PLOG` keywords in reactions. This has been reported but [has not yet been solved](https://bugs.openfoam.org/view.php?id=3523).

> For more information, please refer to the [defination of PLOG](http://engine.princeton.edu/modelreduction/PLOG-documents/PLOG-software_distribution.pdf), and also the `Troe` keyword in the OpenFOAM [source code](https://github.com/OpenFOAM/OpenFOAM-8/blob/master/src/thermophysicalModels/specie/reaction/reactionRate/fallOffFunctions/TroeFallOffFunction/TroeFallOffFunction.H).

In Chemkin, `PLOG` is shown as
```C++
C3H4+H=CH3CCH2 9.2E+38 -8.65 7000. !(This is a dummy reaction rate for Chemkin) 
 PLOG /0.1 9.2E+38 -8.65 7000./ 
 PLOG /1. 9.5E+42 -9.43 11200./ 
 PLOG /10. 1.5E+45 -9.69 15100./ 
 PLOG /100. 1.8E+43 -8.78 16800./ 
 PLOG /1.0E+5 4.4E+09 1.45 2400./ 
```

which is written as
```C++
  type            reversibleArrheniusPLOG;
  reaction        "C3H4+H=CH3CCH2";
  A               9.2E+35;//convert mol to kmol
  beta            -8.65;
  Ta              3522.73274;//convert Ea to Ta, Ta = Ea/R (8.314 J⋅K−1⋅mol−1), cal to J (4.18400 J/cal)
  ArrheniusData
  (
      // convert MPa to Pa
      (1E5  9.2E+35 -8.65 3522.73274   )   // PLOG /p A beta Ta/
      (1E6  9.5E+39 -9.43 5636.372385  )
      (1E7  1.5E+42 -9.69 7599.037769  )
      (1E8  1.8E+40 -8.78 8454.558577  )
      (1E11 4.4E+06  1.45 1207.794082  )
  );
```

# lorentzboost_test
This program randoms some Higgs -> Zg events in the lab frame, then boost the system to compute the Theta angle
The definition of `Theta` is defined as following

Boost **WITHOUT** rotation, `Theta` is the angle between Z boson and z-axis
```
root -l -b -q  boostmultivector.C\(1\)

Rotate the system before boosting, `Theta` is the angle between Z boson and z-axis
```
root -l -b -q  boostmultivector.C\(2\)
```

To boost **WITHOUT** rotation, but theta is the angle between the Z in Higgs' rest frame and the Higgs in lab frame
```
root -l -b -q  boostmultivector.C\(3\)
```

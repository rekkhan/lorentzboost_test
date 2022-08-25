# lorentzboost_test
This program randoms some Higgs -> Zg events in the lab frame, then boost the system to compute the Theta angle
The definition of `Theta` is defined as following

## the event generator
The event generator first generates Z and photon decay from a stationary Higgs boson, the angle between Z and the z-axis is called the "`Decay angle of Z`". This angle is randomly generated.

The whole system will then be boosted a long the Z axis to create an event when a moving Higgs decay into Z and photon. 

The last step randomly rotates the whole system by a random 3D angle to mimic the "real" experiment.

## Usage
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

method 1 and 3 should return the `Decay angle of Z`, as shown in the plots attached.

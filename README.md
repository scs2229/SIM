# SIM models – README

**SIM33** – outputs various structures for crystal number and process contributions, uses Lawson parameterization for fragments from droplet shattering, has the Paukert freezing probability in these droplet shattering tendencies 
- function [Nice, cont, INPt, tauNUC] = SIM_33(etabr, etaagg, etaRSg, etaRSG, etaCOA, etaSH, duration, step, INPstop) 
- function [Nice, cont, INPt, tauNUC] = SIM_33thermo(etabr, etaagg, etaRSg, etaRSG, etaCOA, etaSH, duration, step, INPstop, uz, Tstart, tau) 
- function [Nice, cont, INPt, tauNUC] = SIM_33param(etabr, etaagg, etaRSg, etaRSG, etaCOA, etaSH, duration, step, INPstop, pmax, Fshatter, Fbr, Tminbr, FRS, factor) 

Extra version here which outputs the freezing probability values as well: 
- function [probFRZ, Nice, cont, INPt, tauNUC] = SIM_33prob(etabr, etaagg, etaRSg, etaRSG, etaCOA, etaSH, duration, step, INPstop) 

**SIM34** – outputs enhancements, uses Lawson parameterization for fragments from droplet shattering, has the Paukert freezing probability in these droplet shattering tendencies 
- function enhancement = SIM_34(etabr, etaagg, etaRSg, etaRSG, etaCOA, etaSH, duration, step, INPstop) 
- function enhancement = SIM_34thermo(etabr, etaagg, etaRSg, etaRSG, etaCOA, etaSH, duration, step, INPstop, uz, Tstart, tau) 
- function enhancement = SIM_34param(etabr, etaagg, etaRSg, etaRSG, etaCOA, etaSH, duration, step, INPstop, pmax, Fshatter, Fbr, Tminbr, FRS, factor)

**SIM35** and **SIM36param** – expansion of the SIM33 series, still outputs various structures for crystal number and process contribution and has the Paukert freezing probability in droplet shattering tendencies, coll is Boolean (true for collisional DS, false for non-collisional) and sig indicates the functional form for the fragments coming from droplet shattering ( = 1 for Lawson parameterization, 2 for D^3 dependence within the Lawson parameterization, and 3 for a sigmoid; in this last case Fshatter should be a 3-element array for the 3 parameters of the sigmoid) 
- function [Nice, cont, INPt,tauNUC] = SIM_35param(etabr, etaagg, etaRSg, etaRSG, etaCOA, etaSH, duration, step, INPstop, pmax, Fshatter, Fbr, Tminbr, FRS, factor, coll, sig) 
- function [Nice, cont, INPt, tauNUC] = SIM_35thermo(etabr, etaagg, etaRSg, etaRSG, etaCOA, etaSH, duration, step, INPstop, uz, Tstart, tau, Fshatter, coll, sig) 
- function enhancement = SIM_36param(etabr, etaagg, etaRSg, etaRSG, etaCOA, etaSH, duration, step, INPstop, pmax, Fshatter, Fbr, Tminbr, FRS, factor, coll, sig) 

#!/bin/bash

energy=4
while [ $energy -le 15 ]
do
	cp Analysis.C muon-neutrino/${energy}GeV/
	cp Analysis.C muon-neutrino/anti-${energy}GeV/
	cp Analysis.C electron-neutrino/${energy}GeV/
	cp Analysis.C electron-neutrino/anti-${energy}GeV/
	((energy++))
done


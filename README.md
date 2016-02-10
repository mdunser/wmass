# wmass
run the python wrapper:
=======================

python runWmassAnalyzer.py /store/group/phys_smp/Wmass/perrozzi/ntuples/ntuples_2014_05_23_53X/WMinusPOWHEG/InclWeights/ -b -q 1nd -o minusOutput
python runWmassAnalyzer.py /store/group/phys_smp/Wmass/perrozzi/ntuples/ntuples_2014_05_23_53X/WPlusPOWHEG/InclWeights/  -b -q 1nd -o plusOutput

compile the code in root:
    root
    .L wmassAnalyzer.cc+
makes a wmassAnalyzer_cc.so

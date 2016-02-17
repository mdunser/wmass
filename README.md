# wmass
run the python wrapper:
=======================

python runWmassAnalyzer.py /store/group/phys_smp/Wmass/perrozzi/ntuples/ntuples_2014_05_23_53X/WMinusPOWHEG/InclWeights/ -b -q 1nd -o minusOutput
python runWmassAnalyzer.py /store/group/phys_smp/Wmass/perrozzi/ntuples/ntuples_2014_05_23_53X/WPlusPOWHEG/InclWeights/  -b -q 1nd -o plusOutput

compile the code in root:
    root
    .L wmassAnalyzer.cc+
makes a wmassAnalyzer_cc.so



combine   -M MaxLikelihoodFit --saveNLL datacardShapes.txt -t -1                 --expectSignal=1     -S 1
              do NLL           save it                     use asimov dataset    normalization        with/without systematics (bounds should change)

## or do all masses at once
for i in {0..201}; do combine   -M MaxLikelihoodFit --saveNLL -m $i datacardMass${i}.txt -t -1 --expectSignal=1     -S 1; done


## make asimov toy dataset from nominal mass (ID95)
combine datacardMass95.txt -M GenerateOnly -m 888 -t -1 --expectSignal=1 --saveToys -S 0


## combining 
combine datacard2.txt -M MaxLikelihoodFit --toysFile <fileWithToy>  -t -1 --saveNLL

combine datacardMass${i}.txt -M MaxLikelihoodFit --toysFile higgsCombineTest.GenerateOnly.mH999.123456.root  -t -1  --saveNLL -m 100${i} 

for i in {0..201}; do combine   -M MaxLikelihoodFit --toysFile higgsCombineTest.GenerateOnly.mH999.123456.root --saveNLL -m $i datacardMass${i}.txt -t -1 -S 1; done

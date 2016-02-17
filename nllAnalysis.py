import ROOT, optparse, os, itertools, math
from multiprocessing import Pool


class nllStruct:
    def __init__(self, infile, syst, expSig, sign = 'both'):
        self.ifile = ROOT.TFile(infile)
        self.tree  = self.ifile.Get('limit')
        self.syst  = syst
        self.expS  = expSig
        self.sign  = sign
        self.hname = 'syst%sexpSig%sW%s'%(str(self.syst), str(self.expS).replace('.','p'), self.sign)
        self.histo = ROOT.TH1F(self.hname, self.hname, 1000, -1.*delta, delta)
        self.color = ROOT.kBlue if sign == 'both' else ROOT.kOrange if sign == 'minus' else ROOT.kBlack
        if self.expS != 1:
            self.color = self.color+3
        self.style = (24 if syst else 20) if sign == 'both' else (26 if syst else 22) if sign == 'plus' else (25 if syst else 21)
        self.histo.SetMarkerColor(self.color); self.histo.SetMarkerStyle(self.style); self.histo.SetMarkerSize(0.9);
        self.wname = 'w+w-' if sign=='both' else 'w+' if sign == 'plus' else 'w-'
        self.sname = 'with syst.' if syst else 'w/o syst'
        self.ename = 'expSig. %.3f'%expSig
        self.fname = ' '.join([self.wname,self.sname,self.ename])
        self.fillHisto()

    def fillHisto(self):
        for i,m in enumerate(self.tree):
            #if i%4: continue
            if m.quantileExpected == 0.5:
                nll  = m.nll
                nll0 = m.nll0
                mass = masses[m.mh]
                if 90 < m.mh < 100:
                    print 'filling mass %.4f with nll %.4f and nll0 %.4f'%(mass,nll,nll0)
                if nll*-1 < 100 and nll*-1 >= 0.:
                    self.histo.SetBinContent(self.histo.FindBin(mass-realMass), -1*(nll))
                    self.histo.SetBinError  (self.histo.FindBin(mass-realMass),   0.001 )
                #histo.Fill(mass, -1*(nll))

def runCmd(cmd):
    #os.system('echo combine '+cmd)
    os.system('mkdir -p '+cmd[0])
    os.chdir(cmd[0])
    os.system('combine '+' '.join(cmd[1]))
    return True

def runCombine(njobs):
    pool = Pool(njobs)
    tasks = []
    afsbase = '/afs/cern.ch/user/m/mdunser/wmassSW/inputForCombine/masses/'
    for mi in range(202):
        for doSyst,expSig in itertools.product([0,1],[0.01]):
            syss = 'syst%d'%doSyst
            asis = '' if expSig == 1 else 'expSig%s'%(str(expSig).replace('.','p'))
            cdir = afsbase+syss+asis
            cdirp= afsbase+syss+asis+'Wplus'
            cdirm= afsbase+syss+asis+'Wminus'
            df   = afsbase+'datacardMass{mi}.txt -m {mi}'.format(mi=mi)
            dfp  = afsbase+'datacardMass{mi}plus.txt -m {mi}'.format(mi=mi)
            dfm  = afsbase+'datacardMass{mi}minus.txt -m {mi}'.format(mi=mi)
            opt  = '-M MaxLikelihoodFit --saveNLL -t -1 -S %d'%doSyst
#combine datacardMass95plus.txt  -M GenerateOnly -m 8888 -t -1 --expectSignal=1 --saveToys -S 1
#combine datacardMass95minus.txt -M GenerateOnly -m 9999 -t -1 --expectSignal=1 --saveToys -S 1
            asif = '--toysFile %shiggsCombineTest.GenerateOnly.mH%s.123456.root'  %(afsbase,str(expSig))
            asifp= '--toysFile %shiggsCombineTest.GenerateOnly.mH1%s.123456.root'%(afsbase,str(expSig))
            asifm= '--toysFile %shiggsCombineTest.GenerateOnly.mH-1%s.123456.root'%(afsbase,str(expSig))
            tasks.append([cdir,[df,opt,asif]])
            ## also do the jobs for + and minus
            tasks.append([cdirp,[dfp,opt,asifp]])
            tasks.append([cdirm,[dfm,opt,asifm]])
    print 'submitting %d jobs on %d cores'%(len(tasks), njobs)
    if len(tasks)/njobs > 10: print 'this might take a while'
    pool.map(runCmd, tasks)

def getWidth(pol2, y=0.5):
    p0 = pol2.Value(0)
    p1 = pol2.Value(1)
    p2 = pol2.Value(2)

    p = p1/p2
    q = (p0/p2 - y/p2)

    x1 = -1.*p/2 + math.sqrt(p**2/4 - q)
    x2 = -1.*p/2 - math.sqrt(p**2/4 - q)
    return x2,x1
    

def loadMassWeightIDs(weightfile):
    f = open('weights.info','r')
    ret = {}
    for l in f.readlines():
        ret[int(l.split()[1])] = float(l.split()[3])
    return ret




if __name__ == '__main__':
    from optparse import OptionParser
    usage = """%prog [opts] inputDir

    foobar
    """
    parser = OptionParser(usage=usage)
    parser.add_option("-j", "--njobs", dest="njobs", type="int",
                      default=5, help="number of jobs. default %default");
    parser.add_option("-r", "--runCombine", dest="runCombine", action='store_true',
                      default=False,
                      help="Use N threads");

    global opts, args
    (opts, args) = parser.parse_args()

    if opts.runCombine:
        runCombine(opts.njobs)
    #print asdag

    global realMass, delta, masses
    masses = loadMassWeightIDs('weights.info')
    realMass = 80.386 ## the real mass is 80.385, but massID 95 has 80.386
    delta    =  0.050


    allFlavors = []
    syst1expSig1Both  = nllStruct('inputForCombine/masses/syst1/syst1.root'            , 1, 1         ); allFlavors.append(syst1expSig1Both  )
    syst0expSig1Both  = nllStruct('inputForCombine/masses/syst0/syst0.root'            , 0, 1         ); allFlavors.append(syst0expSig1Both  )
    syst1expSig1Plus  = nllStruct('inputForCombine/masses/syst1Wplus/syst1Wplus.root'  , 1, 1, 'plus' ); allFlavors.append(syst1expSig1Plus  )
    syst0expSig1Plus  = nllStruct('inputForCombine/masses/syst0Wplus/syst0Wplus.root'  , 0, 1, 'plus' ); allFlavors.append(syst0expSig1Plus  )
    syst1expSig1Minus = nllStruct('inputForCombine/masses/syst1Wminus/syst1Wminus.root', 1, 1, 'minus'); allFlavors.append(syst1expSig1Minus )
    syst0expSig1Minus = nllStruct('inputForCombine/masses/syst0Wminus/syst0Wminus.root', 0, 1, 'minus'); allFlavors.append(syst0expSig1Minus )
    ## all with expected signal 0.01
    syst1expSig0p01Both  = nllStruct('inputForCombine/masses/syst1expSig0p01/syst1expSig0p01.root'            , 1, 0.01         ); allFlavors.append(syst1expSig0p01Both  )
    syst0expSig0p01Both  = nllStruct('inputForCombine/masses/syst0expSig0p01/syst0expSig0p01.root'            , 0, 0.01         ); allFlavors.append(syst0expSig0p01Both  )
    syst1expSig0p01Plus  = nllStruct('inputForCombine/masses/syst1expSig0p01Wplus/syst1expSig0p01Wplus.root'  , 1, 0.01, 'plus' ); allFlavors.append(syst1expSig0p01Plus  )
    syst0expSig0p01Plus  = nllStruct('inputForCombine/masses/syst0expSig0p01Wplus/syst0expSig0p01Wplus.root'  , 0, 0.01, 'plus' ); allFlavors.append(syst0expSig0p01Plus  )
    syst1expSig0p01Minus = nllStruct('inputForCombine/masses/syst1expSig0p01Wminus/syst1expSig0p01Wminus.root', 1, 0.01, 'minus'); allFlavors.append(syst1expSig0p01Minus )
    syst0expSig0p01Minus = nllStruct('inputForCombine/masses/syst0expSig0p01Wminus/syst0expSig0p01Wminus.root', 0, 0.01, 'minus'); allFlavors.append(syst0expSig0p01Minus )

    ROOT.gStyle.SetOptStat(0)

    c = ROOT.TCanvas('foo', 'bar', 800, 600)
    leg = ROOT.TLegend(0.7, 0.1, 0.9, 0.4)
    leg.SetNColumns(2)
    for i,f in enumerate(allFlavors):
        if not i:
            f.histo.GetYaxis().SetRangeUser(-0.,10)
            f.histo.SetTitle('a majestic plot')
            f.histo.GetXaxis().SetTitle('m_{test}-m_{W}')
            f.histo.GetYaxis().SetTitle('-nll')
        f.histo.Draw('pe same')
        leg.AddEntry(f.histo, f.fname, 'pe')

    leg.Draw('same')


    #frp_hist11 = hist11.Fit('pol2','S')
    #frp_hist00 = hist00.Fit('pol2','S')
    #fr_hist11 = frp_hist11.Get()
    #fr_hist00 = frp_hist00.Get()

    #width11 = getWidth(fr_hist11)
    #width00 = getWidth(fr_hist00)

    #print width11
    #print width00

    
    

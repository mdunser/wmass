import ROOT, optparse, os, itertools, math
from multiprocessing import Pool


class nllStruct:
    def __init__(self, infile, syst, expSig, sign = 'both', is2D = False):
        self.name  = ('with' if syst else 'without') + ' systematics ' + sign + ' charge'
        self.ifile = ROOT.TFile(infile)
        self.tree  = self.ifile.Get('limit')
        self.syst  = syst
        self.expS  = expSig
        self.is2D  = is2D
        self.strd  = '1D' if not self.is2D else '2D'
        self.sign  = sign
        self.hname = 'syst%sexpSig%sW%s'%(str(self.syst), str(self.expS).replace('.','p'), self.sign)
        self.histo = ROOT.TH1F(self.hname, self.hname, 1000, -1.*delta, delta)
        self.color = ROOT.kBlue if sign == 'both' else ROOT.kOrange if sign == 'minus' else ROOT.kBlack
        if self.is2D: 
            self.color = ROOT.kYellow+1 if sign == 'both' else self.color
        if self.expS != 1: self.color = self.color+3
        self.style = (24 if syst else 20) if sign == 'both' else (26 if syst else 22) if sign == 'plus' else (25 if syst else 21)
        self.histo.SetMarkerColor(self.color); self.histo.SetMarkerStyle(self.style); self.histo.SetMarkerSize(0.9);
        self.wname = 'w+w-' if sign=='both' else 'w+' if sign == 'plus' else 'w-'
        self.sname = 'sys.' if syst else 'no sys'
        self.ename = '%.1f'%(100*expSig)+'%' if expSig != 1 else ''
        self.dname = '1D' if not is2D else '2D'
        self.fname = ' '.join([self.wname,self.sname,self.ename,self.dname])
        self.fillHisto()
        self.getWidth()

    def fillHisto(self):
        nll0norm = 0.
        for m in self.tree:
            if m.quantileExpected == 0.5 and m.mh == 95:
                nll0norm = m.nll0
        for i,m in enumerate(self.tree):
            #if i%4: continue
            if m.quantileExpected == 0.5:
                nll  = m.nll
                nll0 = m.nll0
                mass = masses[m.mh]
                fillval = +1*(nll+nll0)-nll0norm
                if 70 < m.mh < 120:
                    #print 'filling mass %.4f with nll %.4f and nll0 %.4f'%(mass,nll,nll0)
                #if True: #nll*-1 < 100 and nll*-1 >= 0.:
                    self.histo.SetBinContent(self.histo.FindBin(mass-realMass), fillval)
                    self.histo.SetBinError  (self.histo.FindBin(mass-realMass),   0.001 )
                #histo.Fill(mass, -1*(nll))
        
    def getWidth(self, y=0.5):
        self.poly2 = ROOT.TF1('poly2','[0]+[1]*x+[2]*x*x',-1.*delta, delta)
        self.poly2.FixParameter(0,0)
        self.poly2.SetLineColor(self.color); self.poly2.SetLineStyle(1 if self.style < 24 else 2)
        self.frp = self.histo.Fit('poly2','S Q')
        ##frp = self.histo.Fit('pol2','S')
        self.fr = self.frp.Get()
        print 'fit chi2/ndf', self.fr.Chi2(), self.fr.Ndf()
        p0 = self.fr.Value(0)
        p1 = self.fr.Value(1)
        p2 = self.fr.Value(2)
        print p0, p1, p2
    
        p = p1/p2
        q = (p0/p2 - y/p2)
    
        x1 = -1.*p/2 + math.sqrt(p**2/4 - q)
        x2 = -1.*p/2 - math.sqrt(p**2/4 - q)
        print 'mW = %.4f + %.4f - %.4f'%(p0+realMass,x1,x2)
        self.width = (abs(x1)+abs(x2))/2
        return x2,x1
    


def runCmd(cmd):
    #os.system('echo combine '+cmd)
    os.system('mkdir -p '+cmd[0])
    os.chdir(cmd[0])
    os.system('combine '+' '.join(cmd[1]))
    return True

def runCombine(njobs, do2D = False):
    pool = Pool(njobs)
    tasks = []
    str2d = '' if not do2D else '2D'
    fac2d = 2 if do2D else 1 
    afsbase = '/afs/cern.ch/user/m/mdunser/wmassSW/inputForCombine/masses/'
    for mi in range(70,120):
        for doSyst,expSig in itertools.product([0,1],[1]):
            expSigString = '' if expSig == 1 else '_scale%s'%str(expSig).replace('.','p')
            syss = 'syst%d'%doSyst
            asis = '' if expSig == 1 else 'expSig%s'%(str(expSig).replace('.','p'))
            cdir = afsbase+syss+asis+str2d
            cdirp= afsbase+syss+asis+'Wplus'+str2d
            cdirm= afsbase+syss+asis+'Wminus'+str2d
            df   = afsbase+'datacardMass{mi}{s2d}{es}.txt -m {mi}'.format(mi=mi, es=expSigString, s2d=str2d)
            dfp  = afsbase+'datacardMass{mi}plus{s2d}{es}.txt -m {mi}'.format(mi=mi, es=expSigString, s2d=str2d)
            dfm  = afsbase+'datacardMass{mi}minus{s2d}{es}.txt -m {mi}'.format(mi=mi, es=expSigString, s2d=str2d)
            opt  = '-M MaxLikelihoodFit --saveNLL -t -1 -S %d'%doSyst
#combine datacardMass95plus.txt  -M GenerateOnly -m 8888 -t -1 --expectSignal=1 --saveToys -S 1
#combine datacardMass95minus.txt -M GenerateOnly -m 9999 -t -1 --expectSignal=1 --saveToys -S 1
            asif = '--toysFile %shiggsCombineTest.GenerateOnly.mH%s.123456.root' %(afsbase,str(expSig*fac2d)   )
            asifp= '--toysFile %shiggsCombineTest.GenerateOnly.mH%s.123456.root' %(afsbase,str(fac2d*10+expSig))
            asifm= '--toysFile %shiggsCombineTest.GenerateOnly.mH-%s.123456.root'%(afsbase,str(fac2d*10+expSig))
            tasks.append([cdir,[df,opt,asif]])
            ## also do the jobs for + and minus
            tasks.append([cdirp,[dfp,opt,asifp]])
            tasks.append([cdirm,[dfm,opt,asifm]])
    print 'submitting %d jobs on %d cores'%(len(tasks), njobs)
    if len(tasks)/njobs > 10: print 'this might take a while'
    #for i in tasks: print i
    pool.map(runCmd, tasks)


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
    parser.add_option("-2", "--do2D", dest="do2D", action='store_true',
                      default=False,
                      help="perform 2D fits instead of 1D. default: %default");
    parser.add_option("-s", "--summary", dest="doSummary", action='store_true',
                      default=False,
                      help="make the summary plot that mike wants. default: %default");

    global opts, args
    (opts, args) = parser.parse_args()

    if opts.runCombine:
        runCombine(opts.njobs, opts.do2D)
    #print asdag

    global realMass, delta, masses
    masses = loadMassWeightIDs('weights.info')
    realMass = 80.386 ## the real mass is 80.385, but massID 95 has 80.386
    delta    =  0.020


    allFlavors = []
    if not opts.do2D and not opts.doSummary:
        titlestring = '1D'
        syst1expSig1Both  = nllStruct('inputForCombine/masses/syst1/syst1.root'            , 1, 1         ); allFlavors.append(syst1expSig1Both  )
        syst0expSig1Both  = nllStruct('inputForCombine/masses/syst0/syst0.root'            , 0, 1         ); allFlavors.append(syst0expSig1Both  )
        syst1expSig1Plus  = nllStruct('inputForCombine/masses/syst1Wplus/syst1Wplus.root'  , 1, 1, 'plus' ); allFlavors.append(syst1expSig1Plus  )
        syst0expSig1Plus  = nllStruct('inputForCombine/masses/syst0Wplus/syst0Wplus.root'  , 0, 1, 'plus' ); allFlavors.append(syst0expSig1Plus  )
        syst1expSig1Minus = nllStruct('inputForCombine/masses/syst1Wminus/syst1Wminus.root', 1, 1, 'minus'); allFlavors.append(syst1expSig1Minus )
        syst0expSig1Minus = nllStruct('inputForCombine/masses/syst0Wminus/syst0Wminus.root', 0, 1, 'minus'); allFlavors.append(syst0expSig1Minus )
        ## ## all with expected signal 0.01
        ## syst1expSig0p01Both  = nllStruct('inputForCombine/masses/syst1expSig0p01/syst1expSig0p01.root'            , 1, 0.01         ); allFlavors.append(syst1expSig0p01Both  )
        ## syst0expSig0p01Both  = nllStruct('inputForCombine/masses/syst0expSig0p01/syst0expSig0p01.root'            , 0, 0.01         ); allFlavors.append(syst0expSig0p01Both  )
        ## syst1expSig0p01Plus  = nllStruct('inputForCombine/masses/syst1expSig0p01Wplus/syst1expSig0p01Wplus.root'  , 1, 0.01, 'plus' ); allFlavors.append(syst1expSig0p01Plus  )
        ## syst0expSig0p01Plus  = nllStruct('inputForCombine/masses/syst0expSig0p01Wplus/syst0expSig0p01Wplus.root'  , 0, 0.01, 'plus' ); allFlavors.append(syst0expSig0p01Plus  )
        ## syst1expSig0p01Minus = nllStruct('inputForCombine/masses/syst1expSig0p01Wminus/syst1expSig0p01Wminus.root', 1, 0.01, 'minus'); allFlavors.append(syst1expSig0p01Minus )
        ## syst0expSig0p01Minus = nllStruct('inputForCombine/masses/syst0expSig0p01Wminus/syst0expSig0p01Wminus.root', 0, 0.01, 'minus'); allFlavors.append(syst0expSig0p01Minus )
    elif opts.do2D:
        titlestring = '2D'
        syst1expSig1Both2D  = nllStruct('inputForCombine/masses/syst12D/syst12D.root'            , 1, 1, 'both' , opts.do2D); allFlavors.append(syst1expSig1Both2D )
        syst0expSig1Both2D  = nllStruct('inputForCombine/masses/syst02D/syst02D.root'            , 0, 1, 'both' , opts.do2D); allFlavors.append(syst0expSig1Both2D )
        syst1expSig1Plus2D  = nllStruct('inputForCombine/masses/syst1Wplus2D/syst1Wplus2D.root'  , 1, 1, 'plus' , opts.do2D); allFlavors.append(syst1expSig1Plus2D )
        syst0expSig1Plus2D  = nllStruct('inputForCombine/masses/syst0Wplus2D/syst0Wplus2D.root'  , 0, 1, 'plus' , opts.do2D); allFlavors.append(syst0expSig1Plus2D )
        syst1expSig1Minus2D = nllStruct('inputForCombine/masses/syst1Wminus2D/syst1Wminus2D.root', 1, 1, 'minus', opts.do2D); allFlavors.append(syst1expSig1Minus2D)
        syst0expSig1Minus2D = nllStruct('inputForCombine/masses/syst0Wminus2D/syst0Wminus2D.root', 0, 1, 'minus', opts.do2D); allFlavors.append(syst0expSig1Minus2D)
    elif opts.doSummary:
        titlestring = 'summary'
        syst1expSig1Plus    = nllStruct('inputForCombine/masses/syst1Wplus/syst1Wplus.root'  , 1, 1, 'plus' , False); allFlavors.append(syst1expSig1Plus  )
        syst0expSig1Plus    = nllStruct('inputForCombine/masses/syst0Wplus/syst0Wplus.root'  , 0, 1, 'plus' , False); allFlavors.append(syst0expSig1Plus  )
        syst1expSig1Both    = nllStruct('inputForCombine/masses/syst1/syst1.root'            , 1, 1, 'both' , False); allFlavors.append(syst1expSig1Both  )
        syst0expSig1Both    = nllStruct('inputForCombine/masses/syst0/syst0.root'            , 0, 1, 'both' , False); allFlavors.append(syst0expSig1Both  )
        syst1expSig1Both2D  = nllStruct('inputForCombine/masses/syst12D/syst12D.root'        , 1, 1, 'both' , True ); allFlavors.append(syst1expSig1Both2D )
        syst0expSig1Both2D  = nllStruct('inputForCombine/masses/syst02D/syst02D.root'        , 0, 1, 'both' , True ); allFlavors.append(syst0expSig1Both2D )

    ROOT.gStyle.SetOptStat(0)


    c = ROOT.TCanvas('foo', 'bar', 800, 600)
    topr = True; errors = []
    leg = ROOT.TLegend(0.65, 0.7 if topr else 0.1, 0.9, 0.9 if topr else 0.3)
    leg.SetTextSize(0.02)
    leg.SetNColumns(2)
    leg.SetFillColor(ROOT.kWhite)
    for i,f in enumerate(allFlavors):
        print '%-30s %s'%(f.name, f.width)
        if not i:
            f.histo.GetYaxis().SetRangeUser(-0.1,2.5)
            f.histo.SetTitle('a majestic plot - %s'%(titlestring))
            f.histo.GetXaxis().SetTitle('m_{test}-m_{W}')
            f.histo.GetYaxis().SetTitle('nll+(nll0-nll0_{nom})')
        if not i%2:
            try:
                errors.append('#sigma_{PDF} %s %s: %.3f'%(f.sign, f.strd, math.sqrt(f.width**2 - allFlavors[i+1].width**2) ))
            except:
                print 'did not work'
        f.histo.Draw('pe same')
        leg.AddEntry(f.histo, f.fname, 'pe')
    leg.Draw('same')
    line = ROOT.TLine()
    line.SetLineStyle(7);line.SetLineColor(ROOT.kGray+2);
    #line.DrawLine(-1.*delta,0. ,delta,0. );    
    #line.DrawLine(-1.*delta,0.5,delta,0.5);    
    if opts.doSummary:
        line.SetLineStyle(2)
        line.DrawLine(-1.*delta,0.5,-1.*syst1expSig1Both2D.width,0.5)    
        line.DrawLine(-1.*syst1expSig1Both2D.width,-0.1,-1.*syst1expSig1Both2D.width,0.5);    
        line.DrawLine(-1.*delta,0.5,-1.*syst1expSig1Both  .width,0.5);    
        line.DrawLine(-1.*syst1expSig1Both.width,-0.1,-1.*syst1expSig1Both.width,0.5);    
        line.DrawLine(-1.*syst1expSig1Plus.width,-0.1,-1.*syst1expSig1Plus.width,0.5);    
    c.SaveAs('completeStep2_%s.pdf'%(titlestring))
    c.SaveAs('completeStep2_%s.png'%(titlestring))
    c.SaveAs('completeStep2_%s.eps'%(titlestring))

    lat = ROOT.TLatex(); lat.SetTextSize(0.025)
    for i,err in enumerate(errors):
        lat.DrawLatexNDC(0.43, 0.80-i*0.05, err)


#!/usr/bin/env python
import sys, os, copy
import os.path as osp
from copy import deepcopy as dc
from ROOT import TFile, TH1

class histoStruct():
    def __init__(self, rootf, massID, pm, vname, vn, scale=1):
        self.massID = massID
        self.vname  = vname
        self.vn     = vn
        self.rootf  = TFile(rootf,'READ')
        self.pm     = pm
        self.pms    = 'pos'  if pm > 0 else 'neg'
        self.pmsa   = 'plus' if pm > 0 else 'minus'
        self.scale  = scale
        self.scales = '' if scale == 1 else '_scale%s'%str(self.scale).replace('.','p')
        self.setHistos()
        #self.saveInFile()

    def setHistos(self):
        tmp_mthist  = self.rootf.Get('%s_mtEta_m%d_p0'%(self.vname,self.massID)).ProjectionX().Clone('W%sMT' %self.pmsa); tmp_mthist .Scale(self.scale)
        tmp_etahist = self.rootf.Get('%s_mtEta_m%d_p0'%(self.vname,self.massID)).ProjectionY().Clone('W%sETA'%self.pmsa); tmp_etahist.Scale(self.scale)
        setattr(self, 'W%sMT' %self.pmsa, tmp_mthist )
        setattr(self, 'W%sETA'%self.pmsa, tmp_etahist)
        getattr(self, 'W%sMT' %self.pmsa).SetTitle('W%sMT' %self.pmsa); getattr(self, 'W%sMT' %self.pmsa).SetLineColor(1) ;
        getattr(self, 'W%sETA'%self.pmsa).SetTitle('W%sETA'%self.pmsa); getattr(self, 'W%sETA'%self.pmsa).SetLineColor(1);
        for i in range(1,self.vn+1):
            mname = 'W%sMT_v%dUp' %(self.pmsa,i); ename = 'W%sETA_v%dUp' %(self.pmsa,i);
            tmp_mthist  = self.rootf.Get('%s_mtEta_m%d_p%d'%(self.vname,self.massID, i)).ProjectionX().Clone(mname); tmp_mthist .Scale(self.scale)
            tmp_etahist = self.rootf.Get('%s_mtEta_m%d_p%d'%(self.vname,self.massID, i)).ProjectionY().Clone(ename); tmp_etahist.Scale(self.scale)
            setattr(self, mname, tmp_mthist )
            setattr(self, ename, tmp_etahist)
            getattr(self, mname).SetTitle(mname); getattr(self, mname).SetLineColor(2); self.makeDown(getattr(self,mname), 'W%sMT' %(self.pmsa), i);
            getattr(self, ename).SetTitle(ename); getattr(self, ename).SetLineColor(2); self.makeDown(getattr(self,ename), 'W%sETA'%(self.pmsa), i);
        

    def makeUnrolledHisto(self):
        th2 = self.rootf.Get('%s_mtEta_m%d_p0'%(self.vname,self.massID))
        tot_nbins = th2.ProjectionX().GetXaxis().GetNbisX()*th2.ProjectionY().GetXaxis().GetNbinsX()
        tmp_2dhisto = ROOT.TH1F('W%sMTETAunrolled'%self.pmsa, 'W%sMTETAunrolled'%self.pmsa, tot_nbins, 0.5, tot_nbins+0.5)
        xbins = th2.ProjectionX().GetXaxis().GetNbinsX(), ybins = th2.ProjectionY().GetXaxis().GetNbinsX()
        for x in range(xbins+1):
            for y in range(ybins+1):
                cont = th2.GetBinContent(x,y)
        for i in range(1,tot_nbins+1):
            tmp_2dhisto.SetBinContent(i, th2.GetBinContent())
        return tmp_2dhisto

    def makeDown(self, uphisto, typ, ind):
        tmp_nom = getattr(self, typ).Clone('%s_ratio%d'%(typ,ind))
        tmp_dn  = getattr(self, typ).Clone('%s_v%dDown'%(typ,ind))
        tmp_dn.SetLineColor(3)
        tmp_dn.SetTitle('%s_v%dDown'%(typ,ind))
        tmp_up  = uphisto.Clone('%s_uphisto%d'%(typ,ind))
        tmp_nom.Divide(tmp_up)
        tmp_dn.Multiply(tmp_nom)
        setattr(self, tmp_dn.GetName(), tmp_dn)

    def saveInFile(self, other):
        ofname = 'inputForCombine/masses/WplusWminus_massID%s%s.root'%(self.massID, self.scales)
        outfile = TFile(ofname, 'RECREATE')
        outfile.cd()
        otherpmsa = 'minus' if self.pmsa == 'plus' else 'plus'
        getattr(self, 'W%sMT' %self.pmsa).Write(); getattr(other, 'W%sMT' %otherpmsa).Write()
        data_obs = getattr(self, 'W%sMT' %self.pmsa).Clone('data_obs')
        data_obs.Scale(1/data_obs.Integral())
        data_obs.Write()
        getattr(self, 'W%sETA'%self.pmsa).Write(); getattr(other, 'W%sETA'%otherpmsa)
        for i in range(1,self.vn+1):
            mname1 = 'W%sMT_v%dUp' %(self.pmsa, i); ename1 = 'W%sETA_v%dUp' %(self.pmsa, i);
            mname2 = 'W%sMT_v%dUp' %(otherpmsa, i); ename2 = 'W%sETA_v%dUp' %(otherpmsa, i);
            getattr(self , mname1).Write(); getattr(self , mname1.replace('Up','Down')).Write()
            getattr(other, mname2).Write(); getattr(other, mname2.replace('Up','Down')).Write()
            getattr(self , ename1).Write(); getattr(self , ename1.replace('Up','Down')).Write()
            getattr(other, ename2).Write(); getattr(other, ename2.replace('Up','Down')).Write()
        outfile.Close()
        
    def makeDatacard(self, other):
        tmp_file = open('inputForCombine/datacardTemplate.txt','r')
        tmp_data = tmp_file.read()
        tmp_file.close()
        tmp_data = tmp_data.replace('XXX', str(self.massID)+self.scales)
        tmp_data = tmp_data.replace('YYY', 'W%s'%self .pmsa)
        tmp_data = tmp_data.replace('ZZZ', 'W%s'%other.pmsa)
        tmp_data = tmp_data.replace('AAA', str(getattr(self , 'W%sMT' %self .pmsa).Integral(1,getattr(self , 'W%sMT' %self .pmsa).GetNbinsX() )))
        tmp_data = tmp_data.replace('BBB', str(getattr(other, 'W%sMT' %other.pmsa).Integral(1,getattr(other, 'W%sMT' %other.pmsa).GetNbinsX() )))
        for i in range(1, self.vn+1):
            tmp_data += 'v%d    shape    1.6    1.6\n'%i
        outfile = open('inputForCombine/masses/datacardMass%d%s.txt'%(self.massID,self.scales),'w')
        outfile.write(tmp_data)
        outfile.close()

    def makeDatacardSingleSign(self):
        tmp_file = open('inputForCombine/datacardTemplateSingle.txt','r')
        tmp_data = tmp_file.read()
        tmp_file.close()
        tmp_data = tmp_data.replace('XXX', str(self.massID)+self.scales)
        tmp_data = tmp_data.replace('YYY', 'W%s'%self .pmsa)
        tmp_data = tmp_data.replace('AAA', str(getattr(self , 'W%sMT' %self .pmsa).Integral(1,getattr(self , 'W%sMT' %self .pmsa).GetNbinsX() )))
        for i in range(1, self.vn+1):
            tmp_data += 'v%d    shape    1.6\n'%i
        outfile = open('inputForCombine/masses/datacardMass%d%s%s.txt'%(self.massID,self.pmsa,self.scales),'w')
        outfile.write(tmp_data)
        outfile.close()

    
    
    

if __name__ == '__main__':
    from optparse import OptionParser
    usage = """%prog [opts] inputDir
    """
    parser = OptionParser(usage=usage)
    parser.add_option("-s", "--scaleIntegrals", dest="scale", type="float",
                      default=1., help="scale all integrals by this factor default %default");
    (opts, args) = parser.parse_args()
    
    #finalStruct = {}
    for m in range(202):
        print 'at mass', m
        neg = histoStruct('/afs/cern.ch/work/m/mdunser/public/wmass/minusOutputAll20Masses/Wminus_output_all.root', m, -1, 'ct10', 52, opts.scale)
        pos = histoStruct('/afs/cern.ch/work/m/mdunser/public/wmass/plusOutputAll200Masses/Wplus_output_all.root' , m, +1, 'ct10', 52, opts.scale)
        neg.saveInFile  (pos)
        neg.makeDatacard(pos)
        pos.makeDatacardSingleSign()
        neg.makeDatacardSingleSign()
        
    


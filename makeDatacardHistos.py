#!/usr/bin/env python
import sys, os, copy
import os.path as osp
from copy import deepcopy as dc
from ROOT import TFile, TH1

class histoStruct():
    def __init__(self, rootf, massID, pm, vname, vn):
        self.massID = massID
        self.vname  = vname
        self.vn     = vn
        self.rootf  = TFile(rootf,'READ')
        self.pm     = pm
        self.pms    = 'pos'  if pm > 0 else 'neg'
        self.pmsa   = 'plus' if pm > 0 else 'minus'
        self.setHistos()
        #self.saveInFile()

    def setHistos(self):
        setattr(self, 'W%sMT' %self.pmsa, self.rootf.Get('%s_mtEta_m%d_p0'%(self.vname,self.massID)).ProjectionX().Clone('W%sMT' %self.pmsa))
        setattr(self, 'W%sETA'%self.pmsa, self.rootf.Get('%s_mtEta_m%d_p0'%(self.vname,self.massID)).ProjectionY().Clone('W%sETA'%self.pmsa))
        getattr(self, 'W%sMT' %self.pmsa).SetTitle('W%sMT' %self.pmsa); getattr(self, 'W%sMT' %self.pmsa).SetLineColor(1) ;
        getattr(self, 'W%sETA'%self.pmsa).SetTitle('W%sETA'%self.pmsa); getattr(self, 'W%sETA'%self.pmsa).SetLineColor(1);
        for i in range(1,self.vn+1):
            mname = 'W%sMT_v%dUp' %(self.pmsa,i); ename = 'W%sETA_v%dUp' %(self.pmsa,i);
            setattr(self, mname, self.rootf.Get('%s_mtEta_m%d_p%d'%(self.vname,self.massID, i)).ProjectionX().Clone(mname))
            setattr(self, ename, self.rootf.Get('%s_mtEta_m%d_p%d'%(self.vname,self.massID, i)).ProjectionY().Clone(ename))
            getattr(self, mname).SetTitle(mname); getattr(self, mname).SetLineColor(2); self.makeDown(getattr(self,mname), 'W%sMT' %(self.pmsa), i);
            getattr(self, ename).SetTitle(ename); getattr(self, ename).SetLineColor(2); self.makeDown(getattr(self,ename), 'W%sETA'%(self.pmsa), i);

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
        ofname = 'inputForCombine/masses/WplusWminus_massID%s.root'%(self.massID)
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
        tmp_data = tmp_data.replace('XXX', str(self.massID))
        tmp_data = tmp_data.replace('YYY', 'W%s'%self .pmsa)
        tmp_data = tmp_data.replace('ZZZ', 'W%s'%other.pmsa)
        tmp_data = tmp_data.replace('AAA', str(getattr(self , 'W%sMT' %self .pmsa).Integral(1,getattr(self , 'W%sMT' %self .pmsa).GetNbinsX() )))
        tmp_data = tmp_data.replace('BBB', str(getattr(other, 'W%sMT' %other.pmsa).Integral(1,getattr(other, 'W%sMT' %other.pmsa).GetNbinsX() )))
        for i in range(1, self.vn+1):
            tmp_data += 'v%d    shape    1.6    1.6\n'%i
        outfile = open('inputForCombine/masses/datacardMass%d.txt'%(self.massID),'w')
        outfile.write(tmp_data)
        outfile.close()

    def makeDatacardSingleSign(self):
        tmp_file = open('inputForCombine/datacardTemplateSingle.txt','r')
        tmp_data = tmp_file.read()
        tmp_file.close()
        tmp_data = tmp_data.replace('XXX', str(self.massID))
        tmp_data = tmp_data.replace('YYY', 'W%s'%self .pmsa)
        tmp_data = tmp_data.replace('AAA', str(getattr(self , 'W%sMT' %self .pmsa).Integral(1,getattr(self , 'W%sMT' %self .pmsa).GetNbinsX() )))
        for i in range(1, self.vn+1):
            tmp_data += 'v%d    shape    1.6\n'%i
        outfile = open('inputForCombine/masses/datacardMass%d%s.txt'%(self.massID,self.pmsa),'w')
        outfile.write(tmp_data)
        outfile.close()

    
    
    

if __name__ == '__main__':
    from optparse import OptionParser
    usage = """%prog [opts] inputDir
    """
    parser = OptionParser(usage=usage)
    parser.add_option("-m", "--maxEntries", dest="maxEntries", type="int",
                      default=-1, help="Max entries to process");
    (opts, args) = parser.parse_args()
    
    #finalStruct = {}
    for m in range(202):
        print 'at mass', m
        neg = histoStruct('/afs/cern.ch/work/m/mdunser/public/wmass/minusOutputAll20Masses/Wminus_output_all.root', m, -1, 'ct10', 52)
        pos = histoStruct('/afs/cern.ch/work/m/mdunser/public/wmass/plusOutputAll200Masses/Wplus_output_all.root' , m, +1, 'ct10', 52)
        #neg.saveInFile  (pos)
        neg.makeDatacard(pos)
        pos.makeDatacardSingleSign()
        neg.makeDatacardSingleSign()
    


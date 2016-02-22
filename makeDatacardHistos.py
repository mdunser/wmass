#!/usr/bin/env python
import sys, os, copy
import os.path as osp
from copy import deepcopy as dc
from ROOT import TFile, TH1F, TH2F

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
        self.setUnrolledHistos()
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
        

    def setUnrolledHistos(self):
        th2 = self.rootf.Get('%s_mtEta_m%d_p0'%(self.vname,self.massID))
        th2.RebinY(10)
        xbins = th2.GetXaxis().GetNbins() 
        ybins = th2.GetYaxis().GetNbins()
        tot_nbins = xbins*ybins
        tmp_unrolled = TH1F('W%sMTETAunrolled'%self.pmsa, 'W%sMTETAunrolled'%self.pmsa, tot_nbins, 0.5, tot_nbins+0.5)
        for x in range(1,xbins+1):
            for y in range(1,ybins+1):
                cont = th2.GetBinContent(x,y)
                erro = th2.GetBinError  (x,y)
                if cont < 0: print 'negative content at x,y: %d , %d'%(x,y)
                tmp_unrolled.SetBinContent(y+(x-1)*ybins, cont)
                tmp_unrolled.SetBinError  (y+(x-1)*ybins, erro)
        tmp_unrolled.Scale(self.scale);
        setattr(self, 'W%sMTETAunrolled' %self.pmsa, tmp_unrolled)
        print 'total number of bins', tot_nbins
        nfilled = 0
        for i in range(1,tmp_unrolled.GetXaxis().GetNbins()+1):
            if tmp_unrolled.GetBinContent(i) > 0: nfilled+=1
        print 'total number of filled bins', nfilled
        for i in range(1,self.vn+1):
            name = 'W%sMTETAunrolled_v%dUp' %(self.pmsa,i)
            tmp_th2 = self.rootf.Get('%s_mtEta_m%d_p%d'%(self.vname,self.massID, i)).Clone(name);
            tmp_th2.RebinY(10)
            tmp_unrolledi = tmp_unrolled.Clone(name); tmp_unrolledi.SetTitle(name); tmp_unrolledi.Reset();
            for x in range(1,xbins+1):
                for y in range(1,ybins+1):
                    cont = tmp_th2.GetBinContent(x,y)
                    erro = tmp_th2.GetBinError  (x,y)
                    if cont < 0: print 'negative content at variation %d x,y: %d , %d'%(i,x,y)
                    tmp_unrolledi.SetBinContent(y+(x-1)*ybins, cont)
                    tmp_unrolledi.SetBinError  (y+(x-1)*ybins, erro)
            tmp_unrolledi.Scale(self.scale);
            setattr(self, name, tmp_unrolledi )
            getattr(self, name).SetTitle(name); getattr(self, name).SetLineColor(2); self.makeDown(getattr(self,name), 'W%sMTETAunrolled' %(self.pmsa), i);

    def makeDown(self, uphisto, typ, ind):
        tmp_nom = getattr(self, typ).Clone('%s_ratio%d'%(typ,ind))
        tmp_dn  = getattr(self, typ).Clone('%s_v%dDown'%(typ,ind))
        tmp_dn.SetLineColor(3)
        tmp_dn.SetTitle('%s_v%dDown'%(typ,ind))
        tmp_up  = uphisto.Clone('%s_uphisto%d'%(typ,ind))
        tmp_nom.Divide(tmp_up)
        tmp_dn.Multiply(tmp_nom)
        setattr(self, tmp_dn.GetName(), tmp_dn)

    def saveInFile(self, other): ## called on the negative one. hard coded bs here...
        ofname = 'inputForCombine/masses/WplusWminus2D_massID%s%s.root'%(self.massID, self.scales)
        outfile = TFile(ofname, 'RECREATE')
        outfile.cd()
        #otherpmsa = 'minus' if self.pmsa == 'plus' else 'plus'
        for var in ['MT','ETA','MTETAunrolled']:
            pmratio = 1.4; nbinsx = getattr(self, 'W%s%s' %(self.pmsa ,var)).GetNbinsX();
            scale   = getattr(other, 'W%s%s' %(other.pmsa ,var)).Integral(1,nbinsx) / (pmratio*getattr(self, 'W%s%s' %(self.pmsa ,var)).Integral(1,nbinsx))
            getattr(self, 'W%s%s' %(self.pmsa ,var)).Scale(scale); ## scale the histogram to correct +/- ratio
            getattr(self, 'W%s%s' %(self.pmsa ,var)).Write(); getattr(other, 'W%s%s'%(other.pmsa,var)).Write();
            data_obs = getattr(self, 'W%sMT' %self.pmsa).Clone('data_obs')
            data_obs.Scale(1/data_obs.Integral())
            data_obs.Write()
            for v in ['Up','Down']:
                for i in range(1,self.vn+1):
                    sname = 'W%s%s_v%d%s'%(self .pmsa,var,i,v)
                    oname = 'W%s%s_v%d%s'%(other.pmsa,var,i,v)
                    nbinsx = getattr(self, sname).GetNbinsX();
                    scale  = getattr(other, oname).Integral(1,nbinsx) / (pmratio*getattr(self, sname).Integral(1,nbinsx))
                    getattr(self, sname).Scale(scale); ## scale the histogram to correct +/- ratio
                    getattr(self , sname).Write();
                    getattr(other, oname).Write();
        outfile.Close()
        
    def makeDatacard(self, other, do2d = False):
        str2d = '' if not do2d else '2D'
        var = 'MT' if not do2d else 'MTETAunrolled'
        tmp_file = open('inputForCombine/datacardTemplate.txt','r')
        tmp_data = tmp_file.read()
        tmp_file.close()
        tmp_data = tmp_data.replace('XXX', str(self.massID)+self.scales)
        tmp_data = tmp_data.replace('CCC', '2D' if do2d else '2D')
        tmp_data = tmp_data.replace('DDD', 'ETAunrolled' if do2d else '')
        tmp_data = tmp_data.replace('YYY', 'W%s'%self .pmsa)
        tmp_data = tmp_data.replace('ZZZ', 'W%s'%other.pmsa)
        tmp_data = tmp_data.replace('AAA', str(getattr(self , 'W%s%s' %(self .pmsa,var)).Integral(1,getattr(self , 'W%s%s' %(self .pmsa,var)).GetNbinsX() )))
        tmp_data = tmp_data.replace('BBB', str(getattr(other, 'W%s%s' %(other.pmsa,var)).Integral(1,getattr(other, 'W%s%s' %(other.pmsa,var)).GetNbinsX() )))
        for i in range(1, self.vn+1):
            tmp_data += 'v%d    shape    1.6    1.6\n'%i
        outfile = open('inputForCombine/masses/datacardMass%d%s%s.txt'%(self.massID,self.scales,str2d),'w')
        outfile.write(tmp_data)
        outfile.close()

    def makeDatacardSingleSign(self, do2d = False):
        str2d = '' if not do2d else '2D'
        var = 'MT' if not do2d else 'MTETAunrolled'
        tmp_file = open('inputForCombine/datacardTemplateSingle.txt','r')
        tmp_data = tmp_file.read()
        tmp_file.close()
        tmp_data = tmp_data.replace('XXX', str(self.massID)+self.scales)
        tmp_data = tmp_data.replace('CCC', '2D' if do2d else '2D')
        tmp_data = tmp_data.replace('DDD', 'ETAunrolled' if do2d else '')
        tmp_data = tmp_data.replace('YYY', 'W%s'%self .pmsa)
        tmp_data = tmp_data.replace('AAA', str(getattr(self , 'W%s%s' %(self .pmsa,var)).Integral(1,getattr(self , 'W%s%s' %(self .pmsa,var)).GetNbinsX() )))
        for i in range(1, self.vn+1):
            tmp_data += 'v%d    shape    1.6\n'%i
        outfile = open('inputForCombine/masses/datacardMass%d%s%s%s.txt'%(self.massID,self.pmsa,self.scales,str2d),'w')
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
    for m in range(70,121):
        print 'at mass', m
        neg = histoStruct('/afs/cern.ch/work/m/mdunser/public/wmass/minusOutputAll20Masses/Wminus_output_all.root', m, -1, 'ct10', 52, opts.scale)
        pos = histoStruct('/afs/cern.ch/work/m/mdunser/public/wmass/plusOutputAll200Masses/Wplus_output_all.root' , m, +1, 'ct10', 52, opts.scale)
        neg.saveInFile  (pos)
        neg.makeDatacard(pos,False)
        pos.makeDatacardSingleSign(False)
        neg.makeDatacardSingleSign(False)
        
    


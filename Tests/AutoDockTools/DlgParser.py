#############################################################################
#
# Author: Ruth HUEY, Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2002
#
#############################################################################


# $Header: /Users/mp/facil/autodock/git-luna/autodock-cvstar/Tests/AutoDockTools/Attic/DlgParser.py,v 1.2 2005/08/30 16:52:59 rhuey Exp $
#
# $Id: DlgParser.py,v 1.2 2005/08/30 16:52:59 rhuey Exp $
#
#
#
#
#
#
#

"""
This Object parses the result of an AutoDock job and returns a dictionary. 

"""
import os
from string import find, join, replace, split, rfind, strip
import re
import Numeric

from AutoDockTools.ResultParser import ResultParser


class DlgParser(ResultParser):
    """ reads log from a AutoDock docking and return structured data"""

    keywords = ResultParser.keywords + [
        'coords',
        'vdw_energies',
        'estat_energies',
        'total_energies',   #vdw_energies+estat_energies
        'inhib_constant',
        'intermol_energy',  #(1) 
        'internal_energy',  #(2) NB: 1+2->final docked energy
        'torsional_energy', #(3) NB: 1+3->free energy of binding
        'run',
        'parameter_file',   #AD4 parameter library filename
        'include_1_4_interactions',   #AD4 internal energy switch
        ]
        

    def __init__(self, dlgFile=None):
        """selected dlgFile,ok sets which docked conformations to show"""
        ResultParser.__init__(self)
        self.filename = dlgFile
        #set up dict here
        self.getReDict()
        if dlgFile:
            self.filename = os.path.basename(dlgFile)
            self.parse(dlgFile)



    def parse(self, filename):
        """
        uses keys:
            '|           AutoDock' 
            'DPF> outlev'
            'DPF> '
            'INPUT-PDBQ: '
            'DOCKED: '
            'State='
            '^Seeds'
            '^   [1-9]
        after parsing: 
            self.outlev is switch to correctly find coord fields
            self.clist is list of dictionaries for states,
            self.dpfLines are for building DockingParameters object
            self.histogramlines show results histogram
            self.clusterlines have clustering info
            self.modelList is list of dicts of docked: coords plus energy info 
            self.wroteAll is a flag set by 'write_all_cluster_members' or
                    by 'write_all'
        """
        self.filename = filename
        self.version = 3.05
        #reset
        for item in ['clist','clusterlines','dpfLines','histogramlines',\
                     'ligLines','modelList']:
            setattr(self, item, [])

        self.wroteAll = 0
        self.getReDict()
        try:
            dlgptr = open(filename, 'r')
        except:
            raise IOError
        allLines = dlgptr.readlines()
        self.allLines = allLines
        self.match(allLines)
        if self.version==4.0 and hasattr(self, 'energyDict'):
            #print 'updating with energyDict', self.energyDict
            #do an update of energies not written in std section in 4.0
            for c in self.clist:
                run = c['run']
                #print 'update run=', run
                if not self.energyDict.has_key(run):
                    continue
                for k,v in self.energyDict[run].items():
                    c[k] = v
                    #print 'k,v=', k, v
                
        #at the end, check if clist items have estat_energies and not
        #binding_energy and docking_energy
        #for c in self.clist:
            #energyVal = c['estat_energies'][0]
            #if len(c['estat_energies']):
                #if not c.has_key('binding_energy'):
                    #c['binding_energy'] = energyVal
                #if not c.has_key('docking_energy'):
                    #c['binding_energy'] = energyVal

    def getReDict(self):
        if hasattr(self, 'reDict'):
            for k, d in self.reDict.items():
                d['lines'] = []
            return
        self.reDict = {}
        #had to add MODEL|^USER|^ATOM key for outlev -1
        #seems to slow stuff down alot
        self.reKeys = [
            '                 |            AutoDock', #used to correct autodock4 output 
            'DPF> outlev',       #format switch
            'DPF> write_all',    #wroteAll switch
            'INPUT-PDBQ: USER    NEWDPF',         #extra info in case of clustering dlg
            'DPF> ',             #lines for dpo
            'INPUT-PDBQ: ',      #lines for ligand
            'INPUT-PDBQT: ',     #lines for ligand, v4.0
            'DOCKED: ',          #lines for models
            '^MODEL|^USER|^ATOM|^ENDMDL',       #outlev -1 modellines
            'State=',            #lines for conformations
            '^ [1-9]|^  [1-9]|^   [1-9]',   #cluster and histogram lines
            '^Seeds',            #lines for seed
            'Total number of atoms found in PDBQT file', #to get ligand at #
            '^Atom: ID',         #for non-bond table
            ]
        self.reFuncs = [ 
            self.setADVersion,
            self.setOutlev,
            self.setWroteAll,
            self.getNewDpfInfo,
            self.processDpfLines,
            self.processLigLines,
            self.processLigLinesV4,
            self.getModelLines,
            self.getShortModelLines,
            self.getDlgStates,
            self.getClusterInfo,
            self.getSeedInfo,
            self.getLigandAtomCount,
            self.getNonBondTable,
            ]

        for i in range(len(self.reKeys)):
            k = self.reKeys[i]
            dict =  self.reDict[k] = {}
            dict['re'] = re.compile(k)
            dict['lines'] = []
            dict['func'] = self.reFuncs[i]


    def match(self, allLines):
        self.getReDict()
        self.tested = 1
        for i in range(len(allLines)):
            item = allLines[i]
            #if find item, mark it found + don't test
            for k in self.reKeys:
                d = self.reDict[k]
                m = d['re'].match(item)
                if m:
                    d['lines'].append(item) 
                    break
        for k in self.reKeys:
            d = self.reDict[k]
            lines = d['lines']
            if lines==[]:
                input_key = 'INPUT-PDBQ: '
                if k==input_key:
                    if self.version!=4.0:
                        print "!!no lines found for key=", k, "!!"
                else:
                    print "!!no lines found for key=", k, "!!"
            else:
                apply(d['func'], (lines,), {})


    def setADVersion(self,lines):
        if len(lines):
            for l in lines:
                if find(l, 'AutoDock')>-1:
                    ll = split(l)
                    self.version = float(ll[2])
                    break
        else:
            self.version = 3.0


    def setOutlev(self,lines):
        if len(lines):
            ll = split(lines[0])
            self.outlev = int(ll[2])
        else:
            self.outlev = None


    def setWroteAll(self, lines):
        if len(lines):
            self.wroteAll = 1


    def getLigandAtomCount(self, lines):
        self.ligand_atom_count = int(split(lines[0])[-2])
        #print "parsed ligand_atom_count =", self.ligand_atom_count

    def getNonBondTable(self, lines, echo=False):
        ind = self.allLines.index(lines[0])
        ct = self.ligand_atom_count
        nb_lines = self.allLines[ind+2:ind+2+ct]
        self.nb_array = Numeric.zeros((ct, ct))
        if echo:
            for l in nb_lines:
                print l,
        for i in range(ct):
            l = nb_lines[i][10:-1]
            for j in range(ct):
                jind = j*2 + 1
                if l[jind]=='X':
                    self.nb_array[i][j]=1
            if echo:
                print
                print l
                print self.nb_array[i]


    def getSeedInfo(self, lines):
        seeds = []
        for l in lines:
            ilist = split(l)
            slist = []
            for item in ilist[1:]:
                slist.append(int(item))
            seeds.append(tuple(slist))
            #seeds.append(tuple((int(ilist[1]),int(ilist[2]))))
        if len(seeds)==len(self.clist):
            for i in range(len(seeds)):
                conf = self.clist[i]
                if not conf.has_key('run'):
                    continue
                    #outlev -1 prints States in order
                    #s = seeds[i]
                    #conf['rseed1'] = s[0]
                    #if len(s)==2:
                        #conf['rseed2'] = s[1]
                #else:
                    #run = conf['run']
                    ##run is 1-based; seeds indexed starting w/0
                    #s = seeds[run-1]
                    #conf['rseed1'] = s[0]
                    #if len(s)==2:
                        #conf['rseed2'] = s[1]
                run = conf['run']
                #run is 1-based; seeds indexed starting w/0
                s = seeds[run-1]
                conf['rseed1'] = s[0]
                if len(s)==2:
                    conf['rseed2'] = s[1]


    def getClusterInfo(self, lines):
        cl = self.clusterlines
        hl = self.histogramlines
        for l in lines:
            if find(l, 'RANKING')>-1:
                cl.append(l[:-1])
            elif find(l, '#')>-1:
                hl.append(l[:-1])
        if len(cl):
            self.getClusterRecord(cl)
        else:
            self.clusterRecord = None


    def getClusterRecord(self, cl):
        #print 'in getClusterRecord'
        clRecList = []
        curList = []
        #curList gets list of conf info
        curInd = int(split(cl[0])[0])
        ctr = 1
        for l in cl:
            ll = split(l)
            #when built, newList is
            #[Rank,SubRank,Run,DockedEnergy,ClusterRMSD,RefREMSD]
            newList = map(lambda x:int(x),ll[:3])
            #3/29/05
            if self.wroteAll and self.version!=4.0:
                #print "setting run number to ", ctr
                newList[2] = ctr
                ctr = ctr + 1
            newList2 = map(lambda x:float(x),ll[3:-1])
            newList.extend(newList2)
            if newList[0]==curInd:
                curList.append(newList)
            else:
                clRecList.append(curList)
                curList = [newList]
                curInd = newList[0]
        clRecList.append(curList)
        self.clusterRecord = clRecList
        
            

    def processDpfLines(self, lines):
        for l in lines:
            if len(l[:-1])>5:
                self.dpfLines.append(l[5:-1])

    
    def getNewDpfInfo(self, lines):
        if not len(lines):
            return
        dpfLines = []
        keys = []
        for l in lines:
            ll = split(l)
            ind = ll.index('NEWDPF')
            k = ll[ind+1]
            #only add each key once
            if k not in keys:
                dpfLines.append(join(ll[ind+1:]))
                keys.append(k)
        if len(dpfLines):
            self.dpfLines.extend(dpfLines)



    def processLigLines(self, lines):
        #do not use this for version 4.0
        if self.version==4.0:
            return
        ligLINES = []
        foundRun = 0
        for l in lines:
            #in clustering dlg, multiple copies of input-pdbq are present
            if find(l, 'Run')>-1 and foundRun:
                break
            elif find(l, 'Run')>-1:
                foundRun = 1
            else:
                ligLINES.append(l[12:-1])
        #check here to remove lines of just spaces
        nl = []
        for l in ligLINES:
            if len(strip(l)):
                nl.append(l)
        self.ligLines = nl
        #print "in processLigLines:len(ligLines)=", len(nl)
        #self.ligLines = ligLINES

    def processLigLinesV4(self, lines):
        #use this only for version 4.0
        if self.version!=4.0:
            #print "not version 4.0"
            return

        ligLINES = []
        foundRun = 0
        for l in lines:
            #in clustering dlg, multiple copies of input-pdbq are present
            if find(l, 'Run')>-1 and foundRun:
                break
            elif find(l, 'Run')>-1:
                foundRun = 1
            else:
                ligLINES.append(l[13:-1])
        #check here to remove lines of just spaces
        nl = []
        for l in ligLINES:
            if len(strip(l)):
                nl.append(l)
        self.ligLines = nl


    def getModelLines(self, lines):
        if not len(lines):
            return
        modelList = []
        if find(lines[0], 'DOCKED')==0:
            ind = 8
        elif find(lines[0], 'INPUT-PDBQ')==0:
            ind = 12
        #preprocess lines to remove DOCKED or INPUT-PDBQ
        nlines = []
        for l in lines:
            if self.version!=4.0:
                nlines.append(l[ind:-1])
            else:
                if find(l, 'ATOM')>-1:
                    newLine = l[ind:64] + l[66:-1]
                    nlines.append(newLine)
                else:
                    nlines.append(l[ind:-1])
        curMod = [nlines[0]]
        for l in nlines[1:]:
            if find(l, 'MODEL')>-1:
                modelList.append(curMod)
                curMod = [l]
            else:
                curMod.append(l)
        modelList.append(curMod)
        self.modelList = modelList
        self.makeModels(modelList)


    def updateEnergies(self, lines):
        #this is a terrible hack to try and retrieve energies 
        #not written in DOCKED:  USER section for AutoDock4.0
        energyDict = {}
        ulines = []
        for l in lines:
            if find(l, 'USER')==0:
                ulines.append(l)
        for l in ulines:
            if find(l, 'Run')>-1:
                k = int(split(l)[3])
                energyDict[k] = {}
            elif find(l, 'Estimated Free Energy of Binding')>-1:
                energyDict[k]['binding_energy'] = float(split(l)[7])
            elif find(l, 'Final Intermolecular Energy')>-1:
                energyDict[k]['intermol_energy'] = float(split(l)[6])
            elif find(l, 'Final Internal Energy of Ligand')>-1:
                energyDict[k]['internal_energy'] = float(split(l)[8])
            elif find(l, 'Final Docked Energy ')>-1:
                energyDict[k]['docking_energy'] = float(split(l)[5])
            elif find(l, 'Torsional Free Energy ')>-1:
                #this line has parenthesis: eg '9.03)\012'
                energyDict[k]['torsional_energy'] = float(split(l[:-2])[6])
        self.energyDict = energyDict



    def getShortModelLines(self, lines):
        #print 'in getShortModelLines', len(lines)
        if not len(lines):
            return
        if len(self.modelList): 
            if self.version==4.0:
                self.updateEnergies(lines)
            #print 'returning early'
            return
        modelList = []
        curMod = [lines[0][:-1]]
        endmdl = 0
        for l in lines[1:]:
            #use endmdl because cluster dlg format
            # has MODEL then lines w/o MODEL then MODEL
            #again in desc of 1 model...
            if find(l, 'ENDMDL')>-1:
                curMod.append(l)
                modelList.append(curMod)
                endmdl = 1
            elif endmdl:
                curMod = [l[:-1]]
                endmdl = 0
            else:
                curMod.append(l[:-1])
        #don't append after for loop 
        #because it happens in the else
        #print 'len(modelList)=', len(modelList)
        self.modelList = modelList
        self.makeModels(modelList)


    def getDlgStates(self, lines):
        #print "in getDlgStates: len(self.clist)=", len(self.clist)
        if len(self.reDict['^MODEL|^USER|^ATOM|^ENDMDL']['lines']):
            #print 'not building states because models present'
            return
        if len(self.clist):
            return
        for l in lines:
            # in test-1 State= + 17 items: 3 trans, 4quat + ndihe(10) torsions
            xx = split(l)
            # remove possible punctuation
            for ind in range(len(xx)):
                if xx[ind][-1]==',': xx[ind] = xx[ind][:-1]
                if xx[ind][-1]=='.': xx[ind] = xx[ind][:-1]
            #transList = xx[1:4]
            trans = []
            for p in [0,1,2]:
                trans.append(float(xx[p+1]))
            trans = tuple(trans)

            quat = []
            for n in xx[4:8]:
                quat.append(float(n))
            quat = tuple(quat)

            angList = []
            #NB: here torsions are in the same line
            for n in xx[8:]:
                angList.append(float(n))

            #BUILD A DICTIONARY and put it in clist
            d = {}
            d['trn_x'] = trans[0]
            d['trn_y'] = trans[1]
            d['trn_z'] = trans[2]
            d['qtn_nx'] = quat[0]
            d['qtn_ny'] = quat[1]
            d['qtn_nz'] = quat[2]
            d['qtn_ang_deg'] = quat[3]
            d['num_torsions'] = len(angList)
            d['torsion_values'] = angList
            self.clist.append(d)
            #print "added ", len(self.clist), "th conformation"


    def makeModels(self, modelList):
        #print "in makeModels with ", len(modelList), ' models to build'
        clist = []
        for curMod in modelList:
            clist.append(self.makeModel(curMod))
        self.clist = clist
        #print "end of  makeModels: len(self.clist)=", len(self.clist)



    def makeModel(self, lines):
        coords = []
        vdW = []
        Elec = []
        d = {}
        corr = 0
        if self.outlev==-1:
            corr = -1
        binding_energy2 = None
        for l in lines:
            ll = split(l)
            #if find(l, 'MODEL')>-1:
            #    d['num'] = int((ll)[1])
            if find(l, 'Run')>-1 and find(l, 'Rank')==-1:
                d['run'] = int((ll)[3])
            elif find(l, 'Estimated Free Energy of Binding')>-1:
                d['binding_energy'] = float((ll)[7])
            elif find(l, 'Estimated Inhibition Constant')>-1:
                if find(l, 'N/A') < 0:
                    d['inhib_constant'] = float((ll)[6])
            elif find(l, 'Final Docked Energy')>-1:
                d['docking_energy'] = float((ll)[5])
            elif find(l, 'Final Intermolecular Energy')>-1:
                d['intermol_energy'] = float((ll)[6])
            elif find(l, 'Final Internal Energy of Ligand')>-1:
                d['internal_energy'] = float((ll)[8])
            elif find(l, 'Torsional Free Energy')>-1:
                d['torsional_energy'] = float((ll)[6])
            elif find(l, 'NEWDPF tran0')>-1:
                d['trn_x'] = float(ll[3])
                d['trn_y'] = float(ll[4])
                d['trn_z'] = float(ll[5])
            elif find(l, 'NEWDPF quat0')>-1:
                d['qtn_nx'] = float(ll[3])
                d['qtn_ny'] = float(ll[4])
                d['qtn_nz'] = float(ll[5])
                d['qtn_ang_deg'] = float(ll[6])
            elif find(l, 'NEWDPF dihe0')>-1:
                angList = []
                for n in ll[3:]:
                    angList.append(float(n))
                d['torsion_values'] = angList
                d['num_torsions'] = len(angList)
            elif find(l, 'ATOM')>-1:
                coords.append([float(l[30:38]),float(l[38:46]),float(l[46:54])])
                try:
                    vdW.append(float(l[54:60]))
                except:
                    vdW.append(0.0)
                try:
                    Elec.append(float(l[60:66]))
                except:
                    Elec.append(0.0)
                if len(l)>77:
                    try:
                        binding_energy2 = float(l[70:76])
                    except:
                        pass
            elif find(l, 'HETA')>-1:
                coords.append([float(l[30:38]),float(l[38:46]),float(l[46:54])])
                try:
                    vdW.append(float(l[54:60]))
                except:
                    vdW.append(0.0)
                try:
                    Elec.append(float(l[60:66]))
                except:
                    Elec.append(0.0)
            d['coords'] = coords
            d['vdw_energies'] = vdW
            d['estat_energies'] = Elec
            d['total_energies'] = Numeric.array(Numeric.array(vdW)+Numeric.array(Elec)).tolist()
            if binding_energy2 and not d.has_key('binding_energy'):
                d['binding_energy'] = binding_energy2
                d['docking_energy'] = binding_energy2
        return d



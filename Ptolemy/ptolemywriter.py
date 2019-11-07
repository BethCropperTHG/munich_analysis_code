import os
import sys
import numpy as np

asym = "ASYMPTOPIA = 100"


def ptolemywrite(target, reaction, elab, energy, incoming_potential, outgoing_potential, savedir, ptolemydir):


    #need reaction constituents, masses needed for calculating daughter nucleus
    #I'm keeping the 'H's there for now, because they might be necessary for optical model stuff I may or may not do later
    if reaction[1] == 'd':
        incomingparticle = '2H'
    elif reaction[1] == 'p':
        incomingparticle = '1H'

    if reaction[3] == 'd':
        outgoingparticle = '2H'
    elif reaction[3] == 'p':
        outgoingparticle = '1H'
    elif reaction[3] == 't':
        outgoingparticle = '3H'

    if reaction[1] == '3':
        incomingparticle = '3HE'
        outgoingparticle = '4HE'


    #get the daughter nucleus
    #pretty standard, chopping up strings to get a calculation for the mass, and since this is neutron transfer the element stays the same
    daughter = '%s%s'%(int(target[0:3])+int(incomingparticle[0])-int(outgoingparticle[0]), target[3:5])
    reaction_no_brackets = reaction[1:-1]
    e_kev = round(energy*1000,2)

    #make directories. One for all the files for this peak, one for input, output files
    directory_name = 'Ptolemy_%s_%s_%s_elab%s_excitation%s'%(target,reaction_no_brackets,daughter,elab,e_kev)
    os.mkdir('%s%s'%(savedir,directory_name))
    os.chdir('%s%s'%(savedir,directory_name))
    os.mkdir('%s%s/input_files'%(savedir,directory_name))
    os.mkdir('%s%s/output_files'%(savedir,directory_name))
    os.chdir(savedir)

    """
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~WRITE PTOLEMY FILES AND RUN PTOLEMY~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    """

    #loop over l, l = 6 highest in shell model up to 126
    for l in range(0,7):

        #get parity, p = (-1)^l. This is important because If I don't have the parity right ptolemy will get annoyed. LAZY CODE can't even be bothered to work it out!
        p = None
        if l % 2 == 0: #even
            p = '+' 
        else: #odd
            p = '-'

        #need two mss as well! This is because ptolemy needs the j value
        for s in range(0,2):    

            #define j
            j = l + s - 0.5

            #don't want state of negative angular momentum!
            if j<0: continue

            #can't use this, needs to be in format x/2, also need to save a float j for later
            floatj = j
            j = '%s/2'%(int(j*2))

            #also needs to be for different principal quantum numbers. 
            for nodes in range(0,3):

                #pick the right n, l, j

                #first do this for p,d, neutron addition can populate the 50-82 or 82-126 subshell
                #check the nuclear shell model to see these states
                #we can't differentiate spins so I picked a random one for ls that could have either
                if (incomingparticle == '1H' and outgoingparticle == '2H') or (incomingparticle == '3HE' and outgoingparticle == '4HE'):

                    if l == 0 and (j != '1/2' or nodes != 2): continue                   # 2s1/2
                    if l == 1 and (j != '1/2' or nodes != 2): continue #j arbitrary      # 2p1/2, 2p3/2
                    if l == 2 and ((j != '3/2' and j!= '5/2') or nodes != 1): continue #j arbitrary      # 1d3/2, 1d5/2
                    if l == 3 and (j != '5/2' or nodes != 1): continue #j arbitrary      # 1f3/2, 1f5/2
                    if l == 4 and (j != '7/2' or nodes != 0): continue                   # 0g7/2
                    if l == 5 and (j != '11/2' or nodes != 0): continue #j arbitrary     # 0h9/2,0h11/2
                    if l == 6 and (j != '13/2' or nodes != 0): continue                  # 0i13/2 

                if (incomingparticle == '2H' and outgoingparticle == '1H') or (incomingparticle == '2H' and outgoingparticle == '3H'):
                    
                    if l == 0 and (j != '1/2' or nodes != 2): continue              # 2s1/2
                    if l == 1 and (j != '1/2' or nodes != 1): continue #j arbitrary # 1p1/2, 1p3/2
                    if l == 2 and ((j != '3/2' and j!= '5/2') or nodes != 1): continue #j arbitrary      # 1d3/2, 1d5/2
                    if l == 3 and (j != '5/2' or nodes != 0): continue #j arbitrary # 0f5/2, 0f7/2
                    if l == 4 and (j != '7/2' or nodes != 0): continue #j arbitrary # 0g7/2, 0g9/2
                    if l == 5 and (j != '11/2' or nodes != 0): continue             # 0h11/2
                    if l == 6 : continue                                                # no l = 6 for (d,p)

                if (incomingparticle == '2H' or incomingparticle == '1H'):
                    projectile_wavefunction = 'NODES = 0\nR = 1   A = 0.5   WAVEFUNCTION = av18   L = 0'
                    paramset = 'dpsb'
                elif (incomingparticle == '3HE'):
                    projectile_wavefunction = 'wavefunction phiffer nodes=0 l=0 jp=1/2 spfacp=1.6 v=200.93 r=0.88 a=0.69 param1=0.79 param2=0.87 rc=2.0'
                    paramset = 'alpha3'
                else:
                    raise ValueError('Warning: you are trying to use projectile wavefunctions for a beam you don\'t have the wavefuctions for') 

                #write the ptolemy file using parameters defined so far
                ptolemyfile = """reset
r0target
print 0
REACTION: %s%s%s(%s%s %s) ELAB=%s
PARAMETERSET %s labangles r0target maxlextrap=0

PROJECTILE
%s

%s
;
TARGET
nodes=%s l=%s jp=%s r0=1.28 a=0.65 vso=6 rso0=1.10 aso=0.65 rc0=1.3 
;
INCOMING
%s;
OUTGOING
%s;
LABANGLES
ANGLEMIN=0 ANGLEMAX=60 ANGLESTEP=1
;
writens crosssec

end
"""%(target,reaction,daughter,j,p,energy,elab,paramset,projectile_wavefunction,asym,nodes,l,j,incoming_potential,outgoing_potential)

                #now want to save this to a file in the new directory
                os.chdir('%s%s/input_files'%(savedir,directory_name))


                #need to define the name so we can write to different files for each parameter set
                infilename = 'Ptolemy_%s_%s_%s_elab%s_excitation%s_nequals%s_jpequals%s%s.in'%(target,reaction_no_brackets,daughter,elab,e_kev,nodes+1,floatj,p)

                #open the file and write to it
                f = open(infilename,'w')
                f.write(ptolemyfile)
                f.close()

                #ptolemy doesn't like brackets basically, so need to make backslashes in front of certain characters in some of the variables.

                directory_name2 = 'Ptolemy_%s_%s_%s_elab%s_excitation%s'%(target,reaction_no_brackets,daughter,elab,e_kev)
                infileptolemy = '<%s>'%infilename
                outfileptolemy = '%s%s/output_files/Ptolemy_%s_%s_%s_elab%s_excitation%s_nequals%s_jpequals%s%s.out'%(savedir,directory_name2,target,reaction_no_brackets,daughter,elab,e_kev,nodes+1,floatj,p)

                #run ptolemy for this file, specifying outfile path in the outfile directory
                os.system('%sptolemy %s %s'%(ptolemydir, infileptolemy, outfileptolemy))

                #now need to check if the asymptopia is OK
                f = open(outfileptolemy)
                for line in f:
                    if "IN FUTURE RUNS INCREASE ONE OR BOTH ASYMPTOPIA" in line:
                        raise ValueError("Asymptopia too low in ", outfileptolemy)
                    if "IN FUTURE RUNS INCREASE ASYMPTOPIA" in line:
                        raise ValueError("Asymptopia too low in ", outfileptolemy)

                f.close()

                os.chdir(savedir)

    """
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PTCLEAN~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    """


    #at this point all of the ptolemy output files are pretty messy, thankfully somebody's written a scipt to turn it into 2-column data
    #first of all go into the right directory, then run ptclean
    os.chdir('%s%s/output_files'%(savedir,directory_name))
    os.system('%sptclean'%ptolemydir)

    #want to remove the non-cleaned ones because they get in the way and are unnecessary. can remove this if we need them!
    #define the directory name, then loop over files in the directory
    indir = '%s%s/output_files'%(savedir,directory_name)
    ''''''
    for root, dirs, filenames in os.walk(indir):
        for f in filenames:
            #pick out the ones that don't say 'clean' at the end, and remove them
            if f[-5:-1] != 'clea': #yes I know it doesn't have the 'n'
                os.remove(f)
                pass
    ''''''
    return;

def pttestwrite(target, reaction, elab, energy, incoming_potential, outgoing_potential, savedir, ptolemydir, l, inname, outname):



    #need reaction constituents, masses needed for calculating daughter nucleus
    #I'm keeping the 'H's there for now, because they might be necessary for optical model stuff I may or may not do later
    if reaction[1] == 'd':
        incomingparticle = '2H'
    elif reaction[1] == 'p':
        incomingparticle = '1H'

    if reaction[3] == 'd':
        outgoingparticle = '2H'
    elif reaction[3] == 'p':
        outgoingparticle = '1H'
    elif reaction[3] == 't':
        outgoingparticle = '3H'

    if reaction[1] == '3':
        incomingparticle = '3HE'
        outgoingparticle = '4HE'


    #get the daughter nucleus
    #pretty standard, chopping up strings to get a calculation for the mass, and since this is neutron transfer the element stays the same
    daughter = '%s%s'%(int(target[0:3])+int(incomingparticle[0])-int(outgoingparticle[0]), target[3:5])
    reaction_no_brackets = reaction[1:-1]
    e_kev = round(energy*1000,2)

    #make directories. One for all the files for this peak, one for input, output files
    directory_name = 'Ptolemy_%s_%s_%s_elab%s_excitation%s'%(target,reaction_no_brackets,daughter,elab,e_kev)
    os.chdir(savedir)

    """
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~WRITE PTOLEMY FILES AND RUN PTOLEMY~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    """
    p = None
    if l % 2 == 0: #even
        p = '+' 
    else: #odd
        p = '-'

        #need two mss as well! This is because ptolemy needs the j value
    for s in range(0,2):    

        #define j
        j = l + s - 0.5

        #don't want state of negative angular momentum!
        if j<0: continue

        #can't use this, needs to be in format x/2, also need to save a float j for later
        floatj = j
        j = '%s/2'%(int(j*2))

        #also needs to be for different principal quantum numbers. 
        for nodes in range(0,3):

            #pick the right n, l, j
            #first do this for p,d, neutron addition can populate the 50-82 or 82-126 subshell
            #check the nuclear shell model to see these states
            #we can't differentiate spins so I picked a random one for ls that could have either
            if (incomingparticle == '1H' and outgoingparticle == '2H') or (incomingparticle == '3HE' and outgoingparticle == '4HE'):

                if l == 0 and (j != '1/2' or nodes != 2): continue                   # 2s1/2
                if l == 1 and (j != '1/2' or nodes != 2): continue #j arbitrary      # 2p1/2, 2p3/2
                if l == 2 and (j != '3/2' or nodes != 1): continue #j arbitrary      # 1d3/2, 1d5/2
                if l == 3 and (j != '5/2' or nodes != 1): continue #j arbitrary      # 1f3/2, 1f5/2
                if l == 4 and (j != '7/2' or nodes != 0): continue                   # 0g7/2
                if l == 5 and (j != '11/2' or nodes != 0): continue #j arbitrary     # 0h9/2,0h11/2
                if l == 6 and (j != '13/2' or nodes != 0): continue                  # 0i13/2 

            if (incomingparticle == '2H' and outgoingparticle == '1H') or (incomingparticle == '2H' and outgoingparticle == '3H'):
                    
                if l == 0 and (j != '1/2' or nodes != 2): continue              # 2s1/2
                if l == 1 and (j != '1/2' or nodes != 1): continue #j arbitrary # 1p1/2, 1p3/2
                if l == 2 and ((j != '3/2' and j!= '5/2') or nodes != 1): continue #j arbitrary      # 1d3/2, 1d5/2
                if l == 3 and (j != '5/2' or nodes != 0): continue #j arbitrary # 0f5/2, 0f7/2
                if l == 4 and (j != '7/2' or nodes != 0): continue #j arbitrary # 0g7/2, 0g9/2
                if l == 5 and (j != '11/2' or nodes != 0): continue             # 0h11/2
                if l == 6 : continue                                                # no l = 6 for (d,p)

            if (incomingparticle == '2H' or incomingparticle == '1H'):
                    projectile_wavefunction = 'NODES = 0\nR = 1   A = 0.5   WAVEFUNCTION = av18   L = 0'
                    paramset = 'dpsb'
            elif (incomingparticle == '3HE'):
                    projectile_wavefunction = 'wavefunction phiffer nodes=0 l=0 jp=1/2 spfacp=1.6 v=200.93 r=0.88 a=0.69 param1=0.79 param2=0.87 rc=2.0'
                    paramset = 'alpha3'
            else:
                raise ValueError('Warning: you are trying to use projectile wavefunctions for a beam you don\'t have the wavefuctions for') 

                #write the ptolemy file using parameters defined so far
            ptolemyfile = """reset
r0target
print 0
REACTION: %s%s%s(%s%s %s) ELAB=%s
PARAMETERSET %s labangles r0target maxlextrap=0

PROJECTILE
%s

%s
;
TARGET
nodes=%s l=%s jp=%s r0=1.28 a=0.65 vso=6 rso0=1.10 aso=0.65 rc0=1.3 
;
INCOMING
%s;
OUTGOING
%s;
LABANGLES
ANGLEMIN=0 ANGLEMAX=60 ANGLESTEP=1
;
writens crosssec

end
"""%(target,reaction,daughter,j,p,energy,elab,paramset,projectile_wavefunction,asym,nodes,l,j,incoming_potential,outgoing_potential)

      #now want to save this to a file in the new directory
            os.chdir('%s%s/input_files'%(savedir,directory_name))


            #need to define the name so we can write to different files for each parameter set
            infilename = 'Ptolemy_%s_%s_%s_elab%s_excitation%s_nequals%s_jpequals%s%s_%s_%s.in'%(target,reaction_no_brackets,daughter,elab,e_kev,nodes+1,floatj,p,outname, inname)

            #open the file and write to it
            f = open(infilename,'w')
            f.write(ptolemyfile)
            f.close()

                #ptolemy doesn't like brackets basically, so need to make backslashes in front of certain characters in some of the variables.

            directory_name2 = 'Ptolemy_%s_%s_%s_elab%s_excitation%s'%(target,reaction_no_brackets,daughter,elab,e_kev)
            infileptolemy = '<%s>'%infilename
            outfileptolemy = '%s%s/output_files/Ptolemy_%s_%s_%s_elab%s_excitation%s_nequals%s_jpequals%s%s_%s_%s_.out'%(savedir,directory_name2,target,reaction_no_brackets,daughter,elab,e_kev,nodes+1,floatj,p,outname,inname)

                #run ptolemy for this file, specifying outfile path in the outfile directory
            os.system('%sptolemy %s %s'%(ptolemydir, infileptolemy, outfileptolemy))

                #now need to check if the asymptopia is OK
            f = open(outfileptolemy)
            for line in f:
                if "IN FUTURE RUNS INCREASE ONE OR BOTH ASYMPTOPIA" in line:
                    raise ValueError("Asymptopia too low in ", outfileptolemy)
                if "IN FUTURE RUNS INCREASE ASYMPTOPIA" in line:
                    raise ValueError("Asymptopia too low in ", outfileptolemy)

            f.close()

            os.chdir(savedir)

    """
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PTCLEAN~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    """


    #at this point all of the ptolemy output files are pretty messy, thankfully somebody's written a scipt to turn it into 2-column data
    #first of all go into the right directory, then run ptclean
    os.chdir('%s%s/output_files'%(savedir,directory_name))
    os.system('%sptclean'%ptolemydir)

    #want to remove the non-cleaned ones because they get in the way and are unnecessary. can remove this if we need them!
    #define the directory name, then loop over files in the directory
    indir = '%s%s/output_files'%(savedir,directory_name)
    ''''''
    for root, dirs, filenames in os.walk(indir):
        for f in filenames:
            #pick out the ones that don't say 'clean' at the end, and remove them
            if f[-5:-1] != 'clea': #yes I know it doesn't have the 'n'
                os.remove(f)
                pass
    ''''''
    return;





        

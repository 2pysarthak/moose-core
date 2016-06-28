import sys
import moose
print( '[INFO] Using moose from %s' % moose.__file__ )
import pylab
# import numpy as np
import numpy
import matplotlib.pyplot as plt
from matplotlib import cm
import rdesigneur as rd
import xml.etree.ElementTree as ET
import itertools
from scipy import stats

fname = 'bis4'
displayMoogli = False
moogliDistance = 10
moogliDt = 1.0
mootliSequence = '01234'
displayScatter = False
scatterParamToUseInStats = 0
numSeg = 4

# Stim amplitude is unitless
# Stim Width is unitless, defined as multiple of diffusion length.
# Stim Vel is unitless, defined in terms of diffusion length by  time units of diff constt.
# diffConst is defined in terms of square of diffusion length by unit time.
# diffLength here is in SI units: m^2/s
#
# Equations here are:
# Equations here are:
# Adot = 1 -6A + 5A^2 - A^3, or spread out as:
# Adot = k0a + k1a.A + k2a.A.A + k3a.A.A.A + k4a.Ca.A(1+A+10*B) - k5a.A.B
# Bdot = k1b.A - k2b.B
#

params = {
    'k0a':0.1,  # Constant
    'k1a':-5.0,  # Coeff for A
    'k2a':5.0,  # Coeff for A^2
    'k3a':-1.0,  # Coeff for A^3
    'k4a':10.0,  # turnon of A by A and Ca
    'k5a':-5.0,  # Turnoff of A by B
    'k1b':0.01,  # turnon of B by A
    'k2b':-0.01,   # Decay rate of B
    'diffusionLength':1.0e-6,  # Diffusion characteristic length, used as voxel length too.
    'dendDiameter': 10e-6,  # Diameter of section of dendrite in model
    'dendLength': 100e-6,   # Length of section of dendrite in model
    'diffConstA':5.0,       # Diffusion constant of A
    'diffConstB':2.0,       # Diffusion constant of B
    'stimWidth' :1.0,        # Stimulus width in seconds
    'stimAmplitude':1.0,    # Stimulus amplitude, arb units. From FHN review
    'blankVoxelsAtEnd':25,  # of voxels to leave blank at end of cylinder
    # 'blankVoxelsAtEnd':1,
    'preStimTime':10.0,     # Time to run before turning on stimulus.
    'postStimTime':40.0,    # Time to run after stimulus. ~3x decay time
    'settleTime':20.0,    # Settling time of response, after stimulus. 
                          # To include in analysis of total response over 
                          # entire dendrite.
    'fnumber':1,          # Number to append to fname
}
numSpine = 8

def sp( arg, term ):
    return str( params[arg] ) + term

def makeChemProto( name, stimAmpl = 1, diffLength = 1e-6, preStim = 10.0 ):
    # Parameters

    sw = params['stimWidth']
    dca = params['diffConstA'] * diffLength * diffLength
    dcb = params['diffConstB'] * diffLength * diffLength
    num = params['dendLength']/diffLength

    # Objects
    chem = moose.Neutral( '/library/' + name )
    compt = moose.CubeMesh( chem.path + '/dend')
    # compt = moose.CylMesh( chem.path + '/dend')
    # # for i in range(int(num)):
    # compt.x0 = 0
    # compt.x1 = params['dendLength']
    # compt.r0 = 2*params['dendDiameter']
    # compt.r1 = params['dendDiameter']
    # compt.diffLength = diffLength
    A = moose.Pool( compt.path + '/A' )
    B = moose.Pool( compt.path + '/B' )
    Z = moose.BufPool( compt.path + '/Z' )
    Ca = moose.BufPool( compt.path + '/Ca' )
    phase = moose.BufPool( compt.path + '/phase' )
    vel = moose.BufPool( compt.path + '/vel' )
    ampl = moose.BufPool( compt.path + '/ampl' )
    Adot = moose.Function( A.path + '/Adot' )
    Bdot = moose.Function( B.path + '/Bdot' )
    CaStim = moose.Function( Ca.path + '/CaStim' )
    A.diffConst = dca
    B.diffConst = dcb

    # Equations

    Adot.expr = 'x3*(' + sp('k0a', '+ ' ) + sp('k1a','*x1 + ' ) + sp( 'k2a', '*x1*x1 + ') + sp( 'k3a', '*x1*x1*x1 + ') +sp('k4a','*x0*x1/(1+x1+10*x2) + ' ) + sp( 'k5a', '*x1*x2') + ')'

    Bdot.expr = 'x2*(' + sp('k1b', '*x0*x0 + ') + sp('k2b', '*x1' ) +  ')'
    CaStim.expr = 'x2 * exp( -((x0 - t)^2)/(2* ' + str(sw*sw) + ') )'

    print Adot.expr
    print Bdot.expr
    print CaStim.expr

    # Connections
    Adot.x.num = 4
    moose.connect( Ca, 'nOut', Adot.x[0], 'input' )
    moose.connect( A, 'nOut', Adot.x[1], 'input' )
    moose.connect( B, 'nOut', Adot.x[2], 'input' )
    moose.connect( Z, 'nOut', Adot.x[3], 'input' )
    moose.connect( Adot, 'valueOut', A, 'increment' )

    Bdot.x.num = 3
    moose.connect( A, 'nOut', Bdot.x[0], 'input' )
    moose.connect( B, 'nOut', Bdot.x[1], 'input' )
    moose.connect( Z, 'nOut', Bdot.x[2], 'input' )
    moose.connect( Bdot, 'valueOut', B, 'increment' )

    CaStim.x.num = 3
    moose.connect( phase, 'nOut', CaStim.x[0], 'input' )
    moose.connect( vel, 'nOut', CaStim.x[1], 'input' )
    moose.connect( ampl, 'nOut', CaStim.x[2], 'input' )
    moose.connect( CaStim, 'valueOut', Ca, 'setN' )

    return compt

# def makeCompt( name, parent, dx, dy, dia ):
#     RM = 1.0
#     RA = 1.0
#     CM = 0.01
#     EM = -0.065
#     pax = 0
#     pay = 0
#     if ( parent.className == "Compartment" ):
#         pax = parent.x
#         pay = parent.y
#     compt = moose.Compartment( name )
#     compt.x0 = pax
#     compt.y0 = pay
#     compt.z0 = 0
#     compt.x = pax + dx
#     compt.y = pay + dy
#     compt.z = 0
#     compt.diameter = dia
#     clen = np.sqrt( dx * dx + dy * dy )
#     compt.length = clen
#     compt.Rm = RM / (np.pi * dia * clen)
#     compt.Ra = RA * 4.0 * np.pi * clen / ( dia * dia )
#     compt.Cm = CM * np.pi * dia * clen
#     if ( parent.className == "Compartment" ):
#         moose.connect( parent, 'raxial', compt, 'axial' )
#     return compt


def makePassiveSoma( name, length, diameter ):
    # segmentLength = 1e-6
    # segmentDia = 1e-6
    # shaftLength = 1e-6
    # shaftDia = 0.2e-6
    # headLength = 0.5e-6
    # headDia = 0.5e-6
    
    # cell = moose.Neutral( '/library/' + name )
    # model = moose.element( '/library' )
    # prev = makeCompt( '/library/cell/soma', 
    #         model, 0.0, segmentLength, segmentDia )
    # dend = prev
    # for i in range( 0, numSeg ):
    #     nameDend = '/library/cell/dend' + str( i )
    #     dend = makeCompt( nameDend, dend, 0.0, segmentLength, segmentDia )
    #     # name = '/model/cell/shaft' + str( i )
    #     # shaft = makeCompt( name, dend, shaftLength, 0.0, shaftDia )
    #     # name = '/model/cell/head' + str( i )
    #     # head = makeCompt( name, shaft, headLength, 0.0, headDia )
    # return cell
    elecid = moose.Neuron( '/library/' + name )
    # soma = moose.Compartment(elecid.path + '/soma')
    '''This will not work: since we cannot connect messages wherein source and destination are the same.
    Error: SharedFinfo::addMsg: MessageId 5[13]
    Source Element == DestElement == dend
    Recommend that you individually set up messages for the components of the SharedFinfo, 
    to ensure that the direction of messaging is consistent.
    0: Error: Shell::handleAddMsg: Unable to make/connect Msg: Single from dend to dend
    '''
    prev = moose.Compartment( elecid.path + '/soma')
    parent = elecid
    # segmentLength = 1e-6
    # segmentDia = 1e-6
    segmentDia = diameter
    segmentLength = length
    shaftLength = 1e-6
    shaftDia = 0.2e-6
    headLength = 0.5e-6
    headDia = 0.5e-6

    dx = 0.0
    dy = segmentLength
    dia = segmentDia

    RM = 1.0
    RA = 1.0
    CM = 0.01
    EM = -0.065
    pax = 0
    pay = 0
    if ( parent.className == "Compartment" ):
        pax = parent.x
        pay = parent.y
    prev = moose.Compartment( name )
    prev.x0 = pax
    prev.y0 = pay
    prev.z0 = 0
    prev.x = pax + dx
    prev.y = pay + dy
    prev.z = 0
    prev.diameter = dia
    clen = numpy.sqrt( dx * dx + dy * dy )
    prev.length = clen
    prev.Rm = RM / (numpy.pi * dia * clen)
    prev.Ra = RA * 4.0 * numpy.pi * clen / ( dia * dia )
    prev.Cm = CM * numpy.pi * dia * clen
    
    # dia = segmentDia
    parent = prev
    for i in range(4):
        dx = 0.0
        dy = segmentLength
        dia = dia - 1.0
        name = '/library/cell/dend' + str( i )
        # dend = makeCompt( name, dend, 0.0, segmentLength, segmentDia )
        RM = 1.0
        RA = 1.0
        CM = 0.01
        EM = -0.065
        pax = 0
        pay = 0
        if ( parent.className == "Compartment" ):
            pax = parent.x
            pay = parent.y
        compt = moose.Compartment( name )
        compt.x0 = pax
        compt.y0 = pay
        compt.z0 = 0
        compt.x = pax + dx
        compt.y = pay + dy
        compt.z = 0
        compt.diameter = dia
        clen = numpy.sqrt( dx * dx + dy * dy )
        compt.length = clen
        # compt.Rm = RM / (numpy.pi * dia * clen)
        # compt.Ra = RA * 4.0 * numpy.pi * clen / ( dia * dia )
        # compt.Cm = CM * numpy.pi * dia * clen
        if ( parent.className == "Compartment" ):
            moose.connect( parent, 'raxial', compt, 'axial' )
        parent = compt

    # dend_vec = moose.vec( elecid.path + '/soma', n=4, dtype='Compartment' )   
    # dend.diameter = diameter
    # dend.length = length
    # dend.x0 = 0
    # dend.x = length
    # dend.Em = -65e-3 # Leak potential
    # dend.initVm = -65e-3 # Initial membrane potential
    # dend.Rm = 5e9 # Total membrane resistance of the compartment
    # dend.Cm = 1e-12 # Total membrane capacitance of the compartment
    # dend.Ra = 1e6
    # ls = 0.0
    # n = 1.0
    # for ind,dend in enumerate(dend_vec):
    #     dend.diameter = diameter/n
    #     # dl = ls+length/2
    #     dend.x0 = ls
    #     dend.x = ls + length/2
    #     # dend.y0 = ls
    #     # dend.y = dl
    #     dend.length = length/2
    #     # dend.length = np.sqrt( 2*(dl**2))
    #     ls = ls + length/2
    #     # n = n + 1.0
    #     n = n * 2
    #     print("dend path: ", dend.path)
    #     print("dend diameter: ", dend.diameter)
    #     print("dend start: ", dend.x0)
    #     print("dend end: ", dend.x)
    #     if ind!=0:
    #         moose.connect(dend_vec[ind-1], 'raxial', dend_vec[ind], 'axial')
    # moose.connect(dend, 'raxial', dend1, 'axial')
    # moose.connect(dend_vec[0], 'raxial', dend_vec[1], 'axial')
    # moose.connect(dend_vec[1], 'raxial', dend_vec[2], 'axial')
    # moose.connect(dend_vec[2], 'raxial', dend_vec[3], 'axial')
    # moose.connect(dend_vec, 'raxial', '')
    
    

    # dendComp = moose.Compartment( elecid.path + '/soma')
    # dend1 = moose.Compartment( elecid.path + '/dend[1]' )
    # dend2 = moose.Compartment( elecid.path + '/dend[2]' )
    # dend3 = moose.Compartment( elecid.path + '/dend[3]' )
    # dend4 = moose.Compartment( elecid.path + '/dend[4]' )
    # ls = length
    # n = 1.0
    # dendList = [dend1, dend2, dend3, dend4]
    # for ind,dend in enumerate(dendList):
    #     dend.diameter = diameter/n
    #     dl = ls+length/4
    #     dend.x0 = 0
    #     dend.x = dl
    #     dend.y0 = ls
    #     dend.y = dl
    #     dend.length = np.sqrt( 2*(dl**2))
    #     ls = ls + length/4
    #     n = n + 1.0
    #     print("dend path: ", dend.path)
    #     print("dend diameter: ", dend.diameter)
    #     print("dend start: ", dend.x0)
    #     print("dend end: ", dend.x)
    #     # print("dend raxial: ", dend.raxial)
    # moose.connect(dendComp, 'raxial', dend1, 'axial')
    # moose.connect(dend1, 'raxial', dend2, 'axial')
    # moose.connect(dend2, 'raxial', dend3, 'axial')
    # moose.connect(dend3, 'raxial', dend4, 'axial')
    return elecid
    
def extractParms( tab, preStim, endStim, dt ):
    #aoc = sum(tab) - min(tab) * len( tab )
    #aoc = sum(tab) - tab[0] * len( tab )
    start = int( preStim / dt )
    end = int( endStim / dt )
    temp = tab[start:end]
    if displayMoogli:
        aoc = 0.0
        peak = 0.0
    else:
        aoc = sum(temp) - tab[start-1] * len( temp )
        peak = max( temp )
    return [aoc, peak] 

def runTrial( diffusionLength, v, dist, blanks, preStim, postStim ):
    settleTime = params['settleTime']
    vel = moose.vec( '/model/chem/dend/vel' )
    vel.nInit = v * diffusionLength
    #B = moose.vec( '/model/chem/dend/B' )
    #B.nInit = 0.004
    #print 'runTrial(', v, dist, blanks, ')'
    runtime = preStim + postStim + (dist+ blanks*2)/float(v)

    moose.reinit()
    if not displayMoogli:
        moose.start( runtime )
    tabs = moose.vec( '/model/graphs/plot0' )
    #print len(tabs), len( tabs[0].vector )
    stride = int(dist) / numSpine
    #print 'tab number = ', blanks, blanks + dist/2, blanks + dist - 1
    ret = []
    ps = preStim - 4 * params[ 'stimWidth' ]
    for i in range( numSpine ):
        data = tabs[(blanks + i*stride)].vector
        ret.extend( extractParms( data, ps, runtime + settleTime - postStim, tabs[0].dt ) )
    #print ret
    aoc = 0.0
    peak = 0.0
    spatialSum = np.array( [0.0, 0.0] ) # aoc, peak
    settleTime = params['settleTime']
    
    for i in tabs:
        data = i.vector
        spatialSum += np.array( extractParms( data, ps, runtime + settleTime - postStim, tabs[0].dt ) )

    ret.extend( [spatialSum[0], spatialSum[1] ] )
    return ret

def writeXML( seqDtRange, drange, labels, allSimOutput ):
    seqList, seqScore = makeSequence( numSpine )
    sh = allSimOutput.shape
    print "SHAPES = ", sh
    print "seqList = ", seqList
    print "seqScore = ", seqScore
    print 'LEN Labels = ', len(labels)

    assert( len(sh) == 4 and 
            sh[0] == len(seqDtRange) and 
            sh[1] == len(drange) and 
            sh[2] == len(seqList) and 
            sh[3] == len(labels) )
    root = ET.Element( 'twoDplot' )
    yaxis = ET.SubElement( root, 'yaxis' )
    yaxis.set( "title", "seqDt" )
    yaxis.text = ''.join( str(seqDt) + ' ' for seqDt in seqDtRange ) + '\n'
    xaxis = ET.SubElement( root, 'xaxis' )
    xaxis.set( "title", "distance" )
    xaxis.text = ''.join( str(d) + ' ' for d in drange ) + '\n'
    parameters = ET.SubElement( root, 'parameters' )
    p = []
    for j in params:
        p.append( ET.SubElement( parameters, j ) )
        p[-1].text = str( params[j] )

    seqOutput = ET.SubElement( root, 'seqOutput' )
    for iseq in range( len( seqList ) ):
        seqData = ET.SubElement( seqOutput, 'stimSequence' )
        seqData.set( 'index', str( iseq ) )
        seqData.set( 'score', str( seqScore[iseq] ) )
        seqData.text = ''.join( str(j) + ' ' for j in seqList[iseq] )

    simOutput = ET.SubElement( root, 'simOutput' )
    for idt in range( len( seqDtRange ) ):
        dtData = ET.SubElement( simOutput, 'dtData' )
        dtData.set( 'dt', str( seqDtRange[idt]) )
        for idistance in range( len( drange ) ):
            distanceData = ET.SubElement( dtData, 'distanceData' )
            distanceData.set( 'distance', str( drange[idistance] ) )
            for ilab in range( len( labels ) ):
                labelData = ET.SubElement( distanceData, 'labelData' )
                labelData.set( 'title', labels[ilab] )
                y = allSimOutput[idt,idistance,:,ilab]
                labelData.text = ''.join( str(j) + ' ' for j in y )

    tree = ET.ElementTree( root )
    tree.write( fname + '.' + str( int(params['fnumber']) ) +  '.xml' )

def convertSeq( arg ):
    x = int( arg )
    ret = [0] * numSpine
    for i in range( numSpine ):
        ret[-i - 1 ] = x % 10
        x /= 10
    return [ret,]

def makeSequence( numSpines ):
    x = range( numSpines )
    a = list( itertools.permutations( x ) )
    skipWeight = 0.64

    linearity = []
    allSlopes = []
    b = []
    #patterns = [ [1,2,3,4,0],[1,4,2,0,3],[2,0,4,3,1],[4,3,2,1,0],[0,1,2,3,4]]
    #patterns = a[0::2]
    if displayMoogli:
        patterns = convertSeq( moogliPattern )
    else:
        patterns = a[0::1]

    for i in patterns:   # Do every other sample.
        slope, intercept, r_value, p_value, std_err = stats.linregress(x,i)
        '''
        if slope > 0.1:
            b.append( i )
            linearity.append( r_value**2 )
            #print len( b ), i, linearity[-1]
        '''
        b.append( i )
        #linearity.append( np.sign( slope ) * r_value**2 )
        linearity.append( slope * r_value**2 )
    return b, linearity

def scanDistDtGrid():
    seqList, seqScore = makeSequence( numSpine )
    #moose.setClock(18, 0.02 )
    #print "DT = ", moose.element( '/model/graphs/plot0').dt
    seqDtRange = ( 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0 ) # Interval between stimuli, in sec
    drange = ( 5, 10, 15, 20, 25, 30 ) # Distance covered, in diffLen.
    if displayMoogli:
        seqDtRange = ( moogliDt, )
        drange = ( moogliDistance, )
        # seqDtRange = ( 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0 ) # Interval between stimuli, in sec
        # drange = ( 5, 10, 15, 20, 25, 30 )
    allSimOutput = []

    ampl = moose.vec( '/model/chem/dend/ampl' )
    print("ampl: ", len(ampl))
    ampl.nInit = params['stimAmplitude']
    blanks = params['blankVoxelsAtEnd']
    preStim = params['preStimTime']
    postStim = params['postStimTime']
    diffusionLength = params['diffusionLength']
    
    for seqDt in seqDtRange:
        temp2 = []
        temp3 = []
        for d in drange:
            # compt = moose.Compartment( '/model/chem/soma')
            # compt.diameter = params['dendDiameter']/8 
            dend = moose.element( '/model/chem/dend')
            # print("dend diameter: ", dend.diameter)
            moose.le( '/model/chem/dend' )

            phase = moose.vec( '/model/chem/dend/phase' )
            ampl = moose.vec( '/model/chem/dend/ampl' )

            # Old version. Assumes we activate density * d compartments.
            #ampl.nInit = float(v)/float(d)  

            ampl.nInit = params['stimAmplitude']
            Z = moose.vec( '/model/chem/dend/Z' )
            stride = int(d) / numSpine
            Z.nInit = 0
            phase.nInit = 10000.0

            temp = []
            slopes = []
            for seq in seqList:
                print '.',
                sys.stdout.flush()
                for j in range( numSpine ):
                    k = (blanks + j * stride) 
                    Z[ k ].nInit = 1
                    phase[ k ].nInit = preStim + seq[j] * seqDt
                temp.append( runTrial( diffusionLength, seqDt, d, blanks, preStim, postStim ))
                #print temp
            simOut = np.array( temp )
            temp3.append( simOut )
            print seqDt, d #temp[-1], temp[-2], 
        allSimOutput.append( temp3 )
    return allSimOutput, seqDtRange, drange


def main():
    global displayMoogli
    global moogliPattern
    global moogliDt
    global moogliDistance
    # global rdes
    if len( sys.argv ) == 2:
        print "Usage: ", sys.argv[0], " --param value --moogli distance dt sequence"
        quit()

    displayMoogli = False
    for ii in range( len( sys.argv ) ):
        if sys.argv[ii][:2] == '--':
            argName = sys.argv[ii][2:]
            if argName in params:
                params[argName] *= float(sys.argv[ii+1])
            if argName == 'moogli':
                displayMoogli = True
                moogliDistance = float( sys.argv[ii+1] )
                moogliDt = float( sys.argv[ii+2] )
                moogliPattern = sys.argv[ii+3]

    # moose.le('/library')
    # moose.le('/library/soma')
                
    # moose.le( '/model/chem' )

    diffusionLength = params['diffusionLength']
    dendLength = params['dendLength']
    diffusionLength = params['diffusionLength']
    library = moose.Neutral( '/library' )
#def makeChemProto( name, stimAmpl = 1, diffLength = 1e-6, ):
    makeChemProto( 'spread', stimAmpl = 1,
            diffLength = diffusionLength )
    makePassiveSoma( 'cell', params['dendLength'], params['dendDiameter'] )

    rdes = rd.rdesigneur(
        useGssa = False,
        turnOffElec = True,
        chemPlotDt = 0.02,
        #chemDt = 0.01,
        #diffDt = 0.01,
        diffusionLength = diffusionLength,
        cellProto = [['cell', 'soma']],
        chemProto = [['spread', 'spread']],
        chemDistrib = [['spread', 'soma', 'install', '1' ]],
        plotList = [['soma', '1', 'dend/A', 'n', '# of A'],
            ['soma', '1', 'dend/B', 'n', '# of B'],
            ['soma', '1', 'dend/Ca', 'n', '# of Ca']
            ],
        moogList = [
            ['soma', '1', 'dend/A', 'n', 'A n', 0, 4],
            ['soma', '1', 'dend/B', 'n', 'B n', 0, 2.5],
            ['soma', '1', 'dend/Ca', 'n', 'Ca n', 0.0, 1]
        ]
    )
    rdes.buildModel()
    # return

    allSimOutput, seqDtRange, drange = scanDistDtGrid()

    if displayMoogli:
        preStim = params['preStimTime']
        postStim = params['postStimTime']
        blanks = params['blankVoxelsAtEnd']
        runtime =  postStim + preStim + (drange[0]+ blanks*2)/float(seqDtRange[0])
        rdes.displayMoogli( 0.1, runtime, 0.0 )
    else:
        labels = ['aoc0', 'peak0', 'aoc1', 'peak1', 'aoc2', 'peak2', 'aoc3', 'peak3', 'aoc4', 'peak4', 'aocSum', 'peakSum' ]
        writeXML( seqDtRange, drange, labels, np.array( allSimOutput) )

if __name__ == '__main__':
    main()

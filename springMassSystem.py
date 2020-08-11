# Generate a large spring-mass system in Abaqus
# Developed by Rui Zhang (rui.zhang@utdallas.edu), Nov. 25, 2019

## parameters
# User provided parameters
modelname = 'model-1'   # must be consistent with input file

# DO NOT MODIFY ANYTHING BELOW
## Read particles and connectivity
fileHandle = open(modelname + '.geometry', 'r')
numParticles = int(fileHandle.readline())
particleCoord = []
particleLabel = []
particleMass = []
particleInertia = []
for i in range(numParticles):
    line = fileHandle.readline()
    line = line.strip()
    data = line.split()
    particleLabel.append(int(data[0]))
    particleCoord.append([float(data[1]), float(data[2]), float(data[3])])
    particleMass.append(float(data[4]))
    particleInertia.append(float(data[5]))
print 'Total number of particles = ', numParticles
numSprings = int(fileHandle.readline())
springConn = []
springLabel = []
springStiff = []
for i in range(numSprings):
    line = fileHandle.readline()
    line = line.strip()
    data = line.split()
    springLabel.append(int(data[0]))
    springConn.append([int(data[1]), int(data[2])])
    springStiff.append(float(data[3]))
print 'Total number of springs = ', numSprings
fileHandle.close()
## Create Abaqus model
from abaqus import *
from abaqusConstants import *
from caeModules import *
myModel = mdb.Model('Model-1')
## Create parts
partName = []
for i in range(numParticles):
    M = particleMass[i]
    I = particleInertia[i]
    pn = 'Part-' + str(particleLabel[i])
    partName.append(pn)
    p = myModel.Part(name=pn, dimensionality=THREE_D, 
        type=DEFORMABLE_BODY)
    p.ReferencePoint(point=particleCoord[i])
    r = p.referencePoints
    refPoints=(r[1], )
    region=p.Set(referencePoints=refPoints, name='Set-1')
    p.engineeringFeatures.PointMassInertia(name='Inertia-1', 
        region=region, mass=M, i11=I, i22=I, i33=I, 
        alpha=0.0, composite=0.0)
print 'Particles generated'
## Create assembly
a = myModel.rootAssembly
instanceName = []
for i in range(numParticles):
    insn = partName[i] + '-1'
    instanceName.append(insn)
    p = myModel.parts[partName[i]]
    a.Instance(name=insn, part=p, dependent=ON)
## Create springs
region = []
for i in range(numSprings):
    p1 = springConn[i][0] - 1
    p2 = springConn[i][1] - 1
    r1 = a.instances[instanceName[p1]].referencePoints
    refPoints1=(r1[1], )
    rgn1pair0=regionToolset.Region(referencePoints=refPoints1)
    r1 = a.instances[instanceName[p2]].referencePoints
    refPoints1=(r1[1], )
    rgn2pair0=regionToolset.Region(referencePoints=refPoints1)
    region=((rgn1pair0, rgn2pair0),)
    sn = 'Spring-'+str(springLabel[i])
    a.engineeringFeatures.TwoPointSpringDashpot(name=sn, regionPairs=region, 
        axis=NODAL_LINE, springBehavior=ON, springStiffness=springStiff[i], 
        dashpotBehavior=OFF, dashpotCoefficient=0.0)
print 'Springs generated'
session.viewports['Viewport: 1'].setValues(displayedObject=a)

"""""""""""""""""""""""""""
Torsion Calculator:

Takes in input of PDB structure file and allows explicit morphing of structure in xyz plane according to 
proteins possible degrees of freedom given by Ramachandran plot.
"""""""""""""""""""""""""""
from numpy import sin, cos, tan, arccos, arcsin, cross, dot, array
from chimera import openModels, Molecule
def getMols():
	mlist = openModels.list(modelTypes = [Molecule])
	return mlist

"""""""""""""""""""""""""""
Main(molecule,startingRes,startingBond,endingRes,endingBond,vector_plane_volume2Match)
molecule	--- Chimera molecule object
startingRes --- number of res 
startingBond--- string denoting bond pair and 
				direction (A to B)
				example "P,O"
"""""""""""""""""""""""""""
strDict={}
strDict['o5']=u"O5'"


def getEndBonds(mol):
	endList=[]
	resList=[]
	torBond=[]
	for res in mol.residues:
		resList.append(res)
	torBond.append(resList[57].atomsMap[u"O3'"][0])
	torBond.append(resList[57].atomsMap[u"P"][0])
	endList.append(resList[61].atomsMap[u"O5'"][0])
	endList.append(resList[64].atomsMap[u"O5'"][0])
	endList.append(resList[67].atomsMap[u"O5'"][0])
	endList.append(resList[70].atomsMap[u"O5'"][0])
	endList.append(resList[81].atomsMap[u"O5'"][0])
	endList.append(resList[55].atomsMap[u"O5'"][0])
	endList.append(resList[47].atomsMap[u"O5'"][0])
	endList.append(resList[86].atomsMap[u"O5'"][0])
	return torBond, endList

"""""""""""""""""""""""""""""""""
Tools for spacial coordinate,
we have spherical coordinate system where the conventional Phi (z=pcos(phi)) and Theta

(x=psin(phi)cos(theta), y = psin(phi)cos(theta)) is transformed in a function to dependent 

on the torsion angles between the bonds on the backbone chain
"""""""""""""""""""""""""""""""""
#TORSION BOND is just a pair of specific atom objects that we 
#want to calculate torsion for
from sympy import *

def MorphTo(mol,resChainIndecies,xyzPoint):
	residues=mol.residues[resChainIndecies[0]:resChainIndecies[1]+1]
	if resChainIndecies[0]<resChainIndecies[1]:
		torsionBondsList=directionLow2High(residues)
	else:
		torsionBondsList=directionHigh2Low(residues)
	"""Either get torsion phi and psi or just calculate via dihedral angles"""
	"""Will start by seeing how much we need to change the existing Tau to get to point, then add that to existing torsion angle"""
	n=0
	xformCrds={}
	for elem in torsionBondsList:
		xformCrds[str(elem[0])]=elem[0].xformCoord()
		xformCrds[str(elem[1])]=elem[1].xformCoord()
	TauPlusDict={}
	endAtom=torsionBondsList[len(torsionBondsList)-1][1]
	while n<len(torsionBondsList)-2:
		print(str(n)+'  out of  '+str(len(torsionBondsList)))
		print('distance 1-----'+str(distance2(xformCrds[str(torsionBondsList[8][0])],xformCrds[str(torsionBondsList[8][1])])))
		tauValue=GetClosestPoint(torsionBondsList[n],endAtom,xyzPoint,xformCrds)
		TauPlusDict[str(n)]=tauValue
		xformCrds=UpdateTorList(torsionBondsList,n,tauValue,xformCrds)
		print('distance 2-----'+str(distance2(xformCrds[str(torsionBondsList[8][0])],xformCrds[str(torsionBondsList[8][1])])))
		dtown=distance2(endAtom.xformCoord(),xyzPoint)
		print('DISTANCE_______'+str(dtown))
		if dtown<1:
			break
		n+=1
	return TauPlusDict




"""Helper functions for MorphTo"""
def directionLow2High(residueList):
	torsionBondsList=[]
	count=0
	for elem in residueList:
		if count!=0:
			torsionBondsList.append([torsionBondsList[(count*4)-1][1],elem.atomsMap[u"O3'"][0]])
		torsionBondsList.append([elem.atomsMap[u"O3'"][0],elem.atomsMap[u"C3'"][0]])
		torsionBondsList.append([elem.atomsMap[u"C4'"][0],elem.atomsMap[u"C5'"][0]])
		torsionBondsList.append([elem.atomsMap[u"C5'"][0],elem.atomsMap[u"O5'"][0]])
		torsionBondsList.append([elem.atomsMap[u"O5'"][0],elem.atomsMap[u'P'][0]])
		count+=1
	return torsionBondsList

def directionHigh2Low(residues):
	torsionBondsList=[]
	for elem in residueList:
		if count!=0:
			torsionBondsList.append([torsionBondsList[(count*4)-1][0],elem.atomsMap[u"O3'"][0]])
		torsionBondsList.append([elem.atomsMap[u"C3'"][0],elem.atomsMap[u"O3'"][0]])
		torsionBondsList.append([elem.atomsMap[u"C5'"][0],elem.atomsMap[u"C4'"][0]])
		torsionBondsList.append([elem.atomsMap[u"O5'"][0],elem.atomsMap[u"C5'"][0]])
		torsionBondsList.append([elem.atomsMap[u'P'][0],elem.atomsMap[u"O5'"][0]])
		count+=1
	return torsionBondsList

def GetClosestPoint(torsionBond,endAtom,xyzPoint,xformCrds):
	"""Omega=getOmega(torsionBond[0].xformCoord(),torsionBond[1].xformCoord(), endAtom.xformCoord())
	torsionVect=vector(torsionBond[0].xformCoord(), torsionBond[1].xformCoord())
	Phi=getPhi(torsionVect)
	Tau=getTau(torsionBond, endAtom, Omega)
	Theta=getTheta(torsionVect)
	R=float(distance2(torsionBond[1].xformCoord(),endAtom.xformCoord()))
	"""
	from sympy import cos, sin
	print('torsionBond '+str(torsionBond))
	print('endAtom '+str(endAtom))
	print('xyzPoint '+str(xyzPoint))
	print('__________________________')
	print('xformCrds'+str(xformCrds))
	Omega=getOmega(xformCrds[str(torsionBond[0])],xformCrds[str(torsionBond[1])], xformCrds[str(endAtom)])
	torsionVect=vector(xformCrds[str(torsionBond[0])],xformCrds[str(torsionBond[1])])
	Phi=getPhi(torsionVect)
	Tau=getTau(torsionBond, endAtom, Omega,xformCrds)
	Theta=getTheta(torsionVect)
	R=float(distance2(xformCrds[str(torsionBond[0])],xformCrds[str(endAtom)]))
	Y=symbols('Y')
	
	zPlus=R*(cos(Omega))*(cos(Phi))+R*(sin(Omega))*(cos(Y))*(sin(Phi))
	xPlus=R*cos(Omega)*sin(Phi)*cos(Theta)+R*sin(Omega)*cos(Phi)*(-cos(Y))*cos(Theta)+R*sin(Omega)*sin(Theta)*(-sin(Y))
	yPlus=R*cos(Omega)*sin(Phi)*sin(Theta)+R*sin(Omega)*cos(Phi)*(-cos(Y))*sin(Theta)+R*sin(Omega)*cos(Theta)*(sin(Y))
	x1=xformCrds[str(torsionBond[1])][0]+xPlus
	y=xformCrds[str(torsionBond[1])][1]+yPlus
	z=xformCrds[str(torsionBond[1])][2]+zPlus
	distance1=(xyzPoint[0]-x1)**2+(xyzPoint[1]-y)**2+(xyzPoint[2]-z)**2
	print('this is dist______'+str(distance1))
	deriv=diff(distance1,Y)
	print('this is deriv'+str(deriv))
	
	solutionList=solve(deriv)
	if len(solutionList)==0:
		print('zero solutions!!!')
		return Tau
	for sol in solutionList:
		solList=[]
		zPlus=R*(cos(Omega))*(cos(Phi))+R*(sin(Omega))*(cos(sol))*(sin(Phi))
		xPlus=R*cos(Omega)*sin(Phi)*cos(Theta)+R*sin(Omega)*cos(Phi)*(-cos(sol))*cos(Theta)+R*sin(Omega)*sin(Theta)*(-sin(sol))
		yPlus=R*cos(Omega)*sin(Phi)*sin(Theta)+R*sin(Omega)*cos(Phi)*(-cos(sol))*sin(Theta)+R*sin(Omega)*cos(Theta)*(sin(sol))
		x=xformCrds[str(torsionBond[1])][0]+xPlus
		y=xformCrds[str(torsionBond[1])][1]+yPlus
		z=xformCrds[str(torsionBond[1])][2]+zPlus
		solList.append([distance2([x,y,z], [xformCrds[str(endAtom)][0],xformCrds[str(endAtom)][1],xformCrds[str(endAtom)][2]]),sol])
	minnie=solList[0][0]
	if len(solList)==1:
		return solList[0][1]
	for elem in solList:
		if elem[0]<minnie:
			solution=elem[1]
			minnie=elem[0]
	return solution

def UpdateTorList(torsionBondList,n,tauValue,xformCrds):
	jank=torsionBondList[n+1:]
	for elem in jank:
		xformCrds[str(elem[0])]=torsionToXYZ(torsionBondList[n],elem[0],tauValue,xformCrds)
		xformCrds[str(elem[1])]=torsionToXYZ(torsionBondList[n],elem[1],tauValue,xformCrds)
	return xformCrds

def torsionToXYZ(torsionBond, endAtom,tauValue,xformCrds):
	"""Omega=getOmega(torsionBond[0].xformCoord(),torsionBond[1].xformCoord(), endAtom.xformCoord())
	torsionVect=vector(torsionBond[0].xformCoord(), torsionBond[1].xformCoord())
	Phi=getPhi(torsionVect)
	Tau=tauValue
	Theta=getTheta(torsionVect)
	R=float(distance2(torsionBond[1].xformCoord(),endAtom.xformCoord()))
	"""
	Omega=getOmega(xformCrds[str(torsionBond[0])],xformCrds[str(torsionBond[1])], xformCrds[str(endAtom)])
	torsionVect=vector(xformCrds[str(torsionBond[0])],xformCrds[str(torsionBond[1])])
	Phi=getPhi(torsionVect)
	#Tau=getTau(torsionBond, endAtom, Omega,xformCrds)
	Tau=tauValue
	Theta=getTheta(torsionVect)
	R=float(distance2(xformCrds[str(torsionBond[0])],xformCrds[str(endAtom)]))
	zPlus=R*(cos(Omega))*(cos(Phi))+R*(sin(Omega))*(cos(Tau))*(sin(Phi))
	
	xPlus=R*cos(Omega)*sin(Phi)*cos(Theta)+R*sin(Omega)*cos(Phi)*(-cos(Tau))*cos(Theta)+R*sin(Omega)*sin(Theta)*(-sin(Tau))
	                                                             #-1
	yPlus=R*cos(Omega)*sin(Phi)*sin(Theta)+R*sin(Omega)*cos(Phi)*(-cos(Tau))*sin(Theta)+R*sin(Omega)*cos(Theta)*(sin(Tau))
	
	"""
	xPlus=-R*sin(Omega)*sin(Tau)*cos(Phi)*cos(Theta)+R*cos(Omega)*sin(Phi)*cos(Theta)+R*sin(Omega)*sin(Phi)*sin(Theta)*sin(Tau)
	yPlus=-R*sin(Omega)*(-sin(Tau))*sin(Phi)*cos(Theta)+R*cos(Omega)*sin(Phi)*sin(Theta)+R*sin(Omega)*cos(Phi)*sin(Theta)*(-cos(Tau))
	"""
	xo=xformCrds[str(torsionBond[1])][0]
	yo=xformCrds[str(torsionBond[1])][1]
	zo=xformCrds[str(torsionBond[1])][2]

	x=xformCrds[str(torsionBond[1])][0]+xPlus
	y1=xformCrds[str(torsionBond[1])][1]+yPlus
	z=xformCrds[str(torsionBond[1])][2]+zPlus

	x2=xformCrds[str(endAtom)][0]
	y2=xformCrds[str(endAtom)][1]
	z2=xformCrds[str(endAtom)][2]

	y=chimera.Point()
	y[0]=x
	y[1]=y1
	y[2]=z
	"""
	print('DISTANCE '+str(distance2([xo,yo,zo],[x2,y2,z2])))
	print('_____________________________')
	print('_____________________________')
	print('x difference '+str((((xo-x2)**2)**(.5)-((x-xo)**2)**(.5))))
	print('_____________________________')
	print('y difference '+str((((yo-y2)**2)**(.5)-((y-yo)**2)**(.5))))
	print('_____________________________')
	print('z dist difference '+str((((zo-z2)**2)**(.5)-((z-zo)**2)**(.5))))
	print('_____________________________')
	"""
	return y
	
def getTau(torsionBond, endAtom, omega,xformCrds):
	torsionVect=vector(xformCrds[str(torsionBond[0])],xformCrds[str(torsionBond[1])])
	coneVector=vector(xformCrds[str(torsionBond[1])],xformCrds[str(endAtom)])
	needle=orthoganolPart(torsionVect,coneVector)
	"""debugging part compare needle and reference vector to see what angle we are getting"""
	xo=xformCrds[str(torsionBond[1])][0]
	yo=xformCrds[str(torsionBond[1])][1]
	zo=xformCrds[str(torsionBond[1])][2]
	#bar2([xo,yo,zo],[xo+needle[0],yo+needle[1],zo+needle[2]])
	reference=orthoganolPart(torsionVect,[0.0,0.0,1.0])
	#bar2([xo,yo,zo],[xo+reference[0],yo+reference[1],zo+reference[2]])
	dotted=dot(reference,needle)
	Tau=dotted/float(distance(reference)*distance(needle))
	Tau=math.acos(Tau)
	return Tau


def getPhi(bondVector):
	z=[0.0,0.0,1.0]
	dotted=dot(z,bondVector)
	Phi=dotted/float(distance(z)*distance(bondVector))
	Phi=math.acos(Phi)
	return Phi

def getTheta(bondVector):
	x=[1.0,0.0,0.0]
	bondVector=orthoganolPart([0.0,0.0,1.0],bondVector)
	dotted=dot(x,bondVector)
	theta=dotted/float(distance(x)*distance(bondVector))
	theta=math.acos(theta)
	return theta

def vector(startCoord, endCoord):
	vect=[]
	vect.append(endCoord[0]-startCoord[0])
	vect.append(endCoord[1]-startCoord[1])
	vect.append(endCoord[2]-startCoord[2])
	return np.array(vect)

def getOmega(startBondCoord,endBondCoord,atomCoord):
	from numpy import arccos
	vect1_2=vector(startBondCoord, endBondCoord)
	vect2_i=vector(endBondCoord,atomCoord)
	omega=dot(vect1_2,vect2_i)
	omega=omega/float(distance(vect1_2)*distance(vect2_i))
	omega=numpy.array(omega)
	return math.acos(omega)

def distance(xyzList):
	return (xyzList[0]**2+xyzList[1]**2+xyzList[2]**2)**(0.5)


def orthoganolPart(base, vector2Ortho):
	w=array(base)
	v=array(vector2Ortho)
	w2=v-(dot(v,w))/float(dot(w,w))*w
	return w2


def distance2(xyzList1, xyzList2):
	distance=((xyzList1[0]-xyzList2[0])**2+(xyzList1[1]-xyzList2[1])**2+(xyzList1[2]-xyzList2[2])**2)**(0.5)
	return distance



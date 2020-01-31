import sys
import codecs
import numpy as np
import math
import random
from numpy import matrix
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt

########################################################################################
########################### Definition of auxilary functions ###########################
########################################################################################


# Calculate random outcome for Monte Carlo simulation

def R(state, N , const , axis):

	# First we calculate the probabilities of the state being in the +1-measurement state ,p11, and the -1-measurement state, p22
	# We do this by tracing over the density matrix, corresponding to the +1 and -1 eigenstates of the measurement, times the state matrix
	rho=np.matrix( [[state.item(0,0),state.item(0,1)],[state.item(1,0),state.item(1,1)]] )
	eta=4.*const/float(N)
	axismatrixplus= np.matrix( [[(1+axis.item(2))/2,(axis.item(0)-axis.item(1)*1j)/2],[(axis.item(0)+axis.item(1)*1j)/2,(1-axis.item(2))/2]] )
	axismatrixminus= np.matrix( [[(1-axis.item(2))/2,(-axis.item(0)+axis.item(1)*1j)/2],[(-axis.item(0)-axis.item(1)*1j)/2,(1+axis.item(2))/2]] ) 
	p11 = np.absolute(np.trace(np.dot(rho,(axismatrixplus+((1.-eta))*axismatrixminus))))
	p22 = np.absolute(np.trace(np.dot(rho,((eta))*axismatrixminus)))

	# We let the coin flip decides which state is measurement and then return the corresponding random number. 

	coin_flip = np.random.choice(2, p = [np.absolute(p11), np.absolute(p22)]) #aalthough p11 and p22 are real the're formatted as complex 
	if coin_flip == 0:
		r = 1
	else:
		r = -1
	
	return r;


# Generate random number and update the state accordingly

def BU(state, N, const, axis): 

	# This function updates the density matrix according to a measurement of state in direction axis with measurement strength tau and time interval dt.

	# First, we generate a random outcome corresponding to the state and the axis. 
	randomr=R(state,N,const,axis) 


	# Now, we calculate the corresponding measurement matrix M. 

	rho=np.matrix( [[state.item(0,0),state.item(0,1)],[state.item(1,0),state.item(1,1)]] )
	axismatrixplus= np.matrix( [[(1+axis.item(2))/2,(axis.item(0)-axis.item(1)*1j)/2],[(axis.item(0)+axis.item(1)*1j)/2,(1-axis.item(2))/2]] )
	axismatrixminus= np.matrix( [[(1-axis.item(2))/2,(-axis.item(0)+axis.item(1)*1j)/2],[(-axis.item(0)-axis.item(1)*1j)/2,(1+axis.item(2))/2]] ) 

	# Series in cosines
	#g=math.sqrt(const/float(N))
	eta=4.*const/float(N)
	if randomr == 1:
		M = axismatrixplus + (math.sqrt(1.-eta))*axismatrixminus
	else:
		M = math.sqrt(eta)*axismatrixminus
	# Series in angles
	#if randomr == 1:
	#	M = math.cos(math.pi/4.-const/math.sqrt(float(N)))*axismatrixplus + math.cos(math.pi/4.+const/math.sqrt(float(N)))*axismatrixminus
	#else:
	#	M = math.sin(math.pi/4.-const/math.sqrt(float(N)))*axismatrixplus + math.sin(math.pi/4.+const/math.sqrt(float(N)))*axismatrixminus
	# The updated state is the normalized prodect M * rho * M

	newstateprime = np.dot(M,np.dot(rho,M))

	newstate=newstateprime/np.trace(newstateprime)

	return np.matrix([newstate,randomr])


# Arctanget

def newatan(y,x):
	return math.atan2(y,x)+math.pi*(1-np.sign(math.atan2(y,x)))

# geometric phase corresponding a jump, in our chosen gauge

def phasejump(thin,phin,thf,phf):
	return math.atan2(math.sin(phin-phf)*math.tan(thin/2.)*math.tan(thf/2.),1.+math.cos(phin-phf)*math.tan(thin/2.)*math.tan(thf/2.))

# Wheighted average and corresponding standard deviation

def weighted_avg_and_std(values, weights):
    average = np.average(values, weights=weights)
    variance = np.average((values-average)**2, weights=weights)  
    return (average, math.sqrt(variance))


########################################################################################
########################### Prepare and define variables ###########################
########################################################################################

#Define Pauli matrices

sigmax= np.matrix( [[0,1.],[1.,0]] )
sigmay= np.matrix( [[0,-1.j],[1.j,0]] )
sigmaz= np.matrix( [[1.,0],[0,-1.]] )

# Define latitude of measurement

thetam = math.pi/4.

# Give initial state: Eigenstate of initial measurmenet axis inaxis to eigenvalue +1

inaxis = np.array( [math.sin(thetam),0,math.cos(thetam)] )
instate= np.matrix( [[(1+inaxis.item(2))/2,(inaxis.item(0)-inaxis.item(1)*1j)/2],[(inaxis.item(0)+inaxis.item(1)*1j)/2,(1-inaxis.item(2))/2]] )

# Define the number of loops 

n= 1

# Define resolution of one loop

m=500
dt=1./m
N=n*m

# Define the measurement strengths we want to simulate 

valuesconst=np.logspace(-1.6,math.log((N/4.-.1))/math.log(10.),30)


# Prepare the histogram

# Resolution of the phase
res=180
minphase=-math.pi
maxphase=math.pi

histogram1 = []
for k in range(len(valuesconst)):
	histogram1.append([0] * res)

# Add point to histogram

def addpoint1(posi,x,weight):
	result=histogram1
	gridphase=np.linspace(minphase,maxphase,res)
	gridphase=gridphase.tolist()
	gridphase.append(x)
	gridphase.sort()
	pos_phase=gridphase.index(x)
	result[posi][pos_phase] +=weight
	return result

########################################################################################
################################# Perform the Monte Carlo ##############################
########################################################################################

# Do the iteration
allphases=[]
alltheta=[]
allphi=[]

# number of realizations

NR=100

results=[]

# loop over measurement strengths

for l in range(valuesconst.size):
	print("running loop "+ str(l+1) + " out of " + str(valuesconst.size))

	expphase=[]
	const=valuesconst[l]
	weights=[]
	singlecos=[]
	singlesin=[]
	total=[]
	lastthetas=[]
	lastphis=[]

	# loop over the realizations

	for k in range(NR):
	
		# In the lists x,y and z we save the data of the realization

		x=[]
		y=[]
		z=[]
		t=[]
		r=[]

		x.append(np.trace(np.dot(instate,sigmax)))
		y.append(np.trace(np.dot(instate,sigmay)))
		z.append(np.trace(np.dot(instate,sigmaz)))
		t.append(0)
		r.append(0)
		phi=[]
		theta=[]
		phi.append(0.)
		theta.append(thetam)
		phase=0.


		# loop over the quasicontinous measurement sequence

		state= instate
		for i in range(1,N+1):
			
			axis = np.matrix( [math.sin(thetam)*math.cos(2*math.pi*i/m),math.sin(thetam)*math.sin(2*math.pi*i/m),math.cos(thetam)] )
			state=BU(state, m, const, axis).item(0)
			x.append(np.trace(np.dot(state,sigmax)))
			y.append(np.trace(np.dot(state,sigmay)))
			z.append(np.trace(np.dot(state,sigmaz)))
			theta.append(math.acos(np.real(np.trace(np.dot(state,sigmaz)))))
			phi.append(math.atan2(np.real(np.trace(np.dot(state,sigmay))),np.real(np.trace(np.dot(state,sigmax)))))
			t.append(float(i)/float(m))
			r.append(BU(state, m, const, axis).item(1))

			phase +=phasejump(theta[i-1],phi[i-1],theta[i],phi[i])

		lasttheta=theta[N]
		lastphi=phi[N]
		lastphasejump = phasejump(lasttheta,lastphi,thetam,0.)
		phase += lastphasejump
		# correct if phase outside of interval from -pi to pi
		if phase>math.pi:
			phase -= math.pi*2
		if phase<-math.pi:
			phase += math.pi*2
		total.append(phase)
		lastthetas.append(lasttheta)
		lastphis.append(lastphi)
		addpoint1(l,phase,np.absolute(math.cos(thetam/2.)*math.cos(lasttheta/2.)+math.sin(thetam/2.)*math.sin(lasttheta/2.)*math.cos(lastphi) + 1.j*math.sin(thetam/2.)*math.sin(lasttheta/2.)*math.sin(lastphi))**2)
		singlecos.append(math.cos(2.*phase))
		singlesin.append(math.sin(2.*phase))


		#all weighted:
		weights.append(np.absolute(math.cos(thetam/2.)*math.cos(lasttheta/2.)+math.sin(thetam/2.)*math.sin(lasttheta/2.)*math.cos(lastphi) + 1.j*math.sin(thetam/2.)*math.sin(lasttheta/2.)*math.sin(lastphi))**2)

		#exp all weighted:
		expphase.append(np.exp(2.*phase*1.j))

	alltheta.append(lastthetas)
	allphi.append(lastphis)
	allphases.append(total)

	averagingcos=weighted_avg_and_std(singlecos,weights)
	averagingsin=weighted_avg_and_std(singlesin,weights)
	averaging2=np.average(expphase,weights=weights)
	propav=sum(weights)/float(NR)


	results.append([const,averagingcos[0],averagingcos[1],averagingsin[0],averagingsin[1],np.angle(averaging2)/2.,propav,propav*(averagingcos[0]**2+averagingsin[0]**2)])
# The table has rows like this for all different measurement strengths: 
# [measurement strength, wheighted average of cos(phase), wheighted std of cos(phase), 
# wheighted average of sin(phase), wheighted std of sin(phase), phase/2 of average of exp(2 i phase), 
# average probability of succesfull final projection,  average prob of succ of final projection times abs( average of exp(2 i phase)) ]
np.savetxt('plotting_phasespaper.txt',results)

# histogram of phases
np.savetxt('whistrogram_paper.txt',histogram1)
# all final phases
np.savetxt('allphases.txt',allphases)
# all final polar angles
np.savetxt('alltheta.txt',alltheta)
# all final azimutal angles
np.savetxt('allphi.txt',allphi)

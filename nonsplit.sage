# -*- coding: latin-1 -*-

########################################################################
#
# This code computes integral points on the modular curve X_ns^+(p),
# for rational primes p = 7, 11, 13, ....
# The corresponding main method is
#
# - integralPoints_on_XnsPlus_P(p,...).
#
# Authors: Aurelien Bajolet, Yuri Bilu, Benjamin Matschke
# Licence: Creative commons BY-NC.
#
# The code is mostly explained in the authors' paper
# "Computing Integral Points on $X_{ns}^+(p)$" (second arxiv version).
# The first arXiv version was based on pari/GP code.
#
########################################################################

pathData = pathData = "data/";

########################################################################
### Parallel computing settings: #######################################
########################################################################

#maxCPUs = 1;
#maxCPUs = 128;
maxCPUs = 12; #Adjusted for plafrim

numCPUs = min(maxCPUs,sage.parallel.ncpus.ncpus()); #the latter is number of available cpus.
parallelIterator = "fork"; #Parallel computing.
##parallelIterator = "multiprocessing"; #Parallel computing, but is not generic enough for us.

#numCPUs = 1;
#parallelIterator = "reference"; #Restricts computation to single cpu.

########################################################################
### Numerics: ##########################################################
########################################################################

precision = 200; #Used as the basic precision, gets automatically increased locally, if needed.
RIF = RealIntervalField(precision);
CIF = ComplexIntervalField(precision);
CC = ComplexField(precision);
RR = RealField(precision);

def myRRtoRIF(x,precisionLoss = 5):
	'''
	Transforms an elmeent from RR to an element in RIF. 
	We assume that at most the last precisionLoss number of bits of x are wrong.
	More precisely we assume that x lies in the following interval:
	  [x-x*2^(precisionLoss-x.precision()),x+x*2^(precisionLoss-x.precision())]
	'''
	return RealIntervalField(x.precision())(x-x*2^(precisionLoss-x.precision()),x+x*2^(precisionLoss-x.precision()));

def myCCtoCIF(z,precisionLoss = 5):
	'''
	Transformas an element from CC to an element in CIF 
	using the myRRtoRIF method for the real and imaginary parts separately.
	'''
	return ComplexIntervalField(z.prec())(myRRtoRIF(z.real(),precisionLoss=precisionLoss),myRRtoRIF(z.imag(),precisionLoss=precisionLoss));

def myRIF_radius(x):
	'''
	Returns an upper bound (in RIF) for the ratius of x in RIF.
	'''
	RIFprec = RealIntervalField(x.prec());
	c = RIFprec(x.center());
	return RIFprec(abs(x-c).upper());

def myRIF_center(x):
	'''
	Returns a center (in RIF) of x in RIF.
	It does NOT necessarily contain the exact center of x;
	it is only an approximation.
	'''
	return RealIntervalField(x.prec())(x.center());

def integralPoints_in_realInterval(realInterval):
	'''
	Returns the list of integers that lie in a real interval.
	'''
	return [ZZ(x) for x in range(realInterval.lower().ceil(),realInterval.upper().floor()+1)];

def LLL_gram_rif(Q,precision=20):
	'''
	Suppose Q is a positive definite matrix with RIF coefficients.
	Returns an integral matrix U, such that
	 U.transpose() * Q * U is a roughly LLL-reduced Gram matrix.
	'''
	n = Q.nrows();
	Qlist = Q.list();
	m = max([abs(x) for x in Qlist]);
	#print "debug15",Q.dimensions(), [(2^precision/m*x).center().ceil() for x in Q.list()]
	QZ = matrix(ZZ,n,n,[(2^precision/m*x).center().round() for x in Qlist]);
	while not QZ.is_positive_definite():
		QZ += identity_matrix(n); #Faster alternative would be to add the smallest eigenvalue (plus something small and positive) times the identity matrix to QZ in order to make it positive definite.
	U = QZ.LLL_gram();
	#print "debug13",U;
	return U;

def myRIFtoString(x):
	'''
	Converts an element of RIF to a string.
	'''
	return "["+str(x.lower().numerical_approx(digits=5))+","+str(x.upper().numerical_approx(digits=5))+"]";

def contains_only_finite_rifs(someIterable):
	'''
	Checks whether the given iterable contains only finite real intervals.	
	'''
	for x in someIterable:
		if x.is_NaN():
			return False;
		if x.lower()==-Infinity:
			return False;
		if x.upper()==+Infinity:
			return False;
	return True;

def lowerAndUpperBoundsForSmallestAndLargestEigenvalueOfSelfadjointRealMatrix(M,verbose=False):
	'''
	We assume that M.base_ring() is RealField(prec) or RIF(prec) for some natural number prec.
	We assume that M is self-adjoint, or at least that all roots are real.
	The result will be an element in RR(prec), where prec is the precision of the base ring of M.
	'''
	
	#Convert matrix elements into real intervals:
	prec = M.base_ring().precision();
	if M.base_ring().name().startswith('IntervalRealIntervalField'):
		Mi = M;
	elif M.base_ring().name().startswith('RealField'):
		Mi = myRRtoRIF_matrix(M);
	else:
		print "Error: M.base_ring() is neither RR nor RIF.";
		return;

	#Do Newton algorithm for the minimal polynomial of M, starting left of all roots.
	#This works so easily because all roots are real!
	Q = Mi.charpoly();
	Qd = Q.derivative();
	cauchyBoundForRoots = 1 + max([abs(c) for c in Q.coefficients()]);
	lamdaMinQ = -cauchyBoundForRoots;
	lamdaMaxQ = +cauchyBoundForRoots;
	if verbose:
		print "lamdaMinQ obtained from Cauchy bound:",lamdaMinQ;
		print "lamdaMaxQ obtained from Cauchy bound:",lamdaMaxQ;

	while True:
		lamdaMinQNew = lamdaMinQ - Q(lamdaMinQ)/Qd(lamdaMinQ);
		lamdaMinQNew = RIF(lamdaMinQNew.lower()); #If we don't do this, then sometimes the error (diameter of the real interval) would explode exponentially, and adjusting the precision slightly will not help.
		if verbose:
			print "New lamdaMinQNew:",lamdaMinQNew;
		if lamdaMinQNew.lower()<=lamdaMinQ.lower():
			break;
		else:
			lamdaMinQ = lamdaMinQNew;

	while True:
		lamdaMaxQNew = lamdaMaxQ - Q(lamdaMaxQ)/Qd(lamdaMaxQ);
		lamdaMaxQNew = RIF(lamdaMaxQNew.upper()); #If we don't do this, then sometimes the error (diameter of the real interval) would explode exponentially, and adjusting the precision slightly will not help.
		if verbose:
			print "New lamdaMaxQNew:",lamdaMaxQNew;
		if lamdaMaxQNew.upper()>=lamdaMaxQ.upper():
			break;
		else:
			lamdaMaxQ = lamdaMaxQNew;

	#If M is the exact height pairing matrix, and M' the above machine approximation,
	#and assuming we can compute the smallest eigenvalue of M' precisely within the machine precision,
	#this is the same as the smallest eigenvalue of M up to an error which is bounded by
	#the spectral radius of the matrix M-M', by a theorem of Ronn ["Bounds on Eigenvalues of Interval Matrices"].
	#This spectral radius is trivially bounded by the maximal absolute value of an entry of M-M' times the dimension of M.

	lamdaError = Mi.nrows()*max([x.absolute_diameter() for x in Mi.list()]); 
	#Now we replace the computed lamda by a lower bound of it:
	lamdaMinQ -= lamdaError;
	lamdaMaxQ += lamdaError;
	return (lamdaMinQ.lower(),lamdaMaxQ.upper());

def rifSignToString(x):
	'''
	Returns the sign (as a string) of the center of the real interval x.
	'''
	sign = x.center().sign();
	if sign>0:
		return "+";
	elif sign<0:
		return "-";
	else:
		return "0";

class PrecisionError(Exception):
    '''
    Exception raised for errors coming from too low precision.

    Attributes:
        msg  -- explanation of the error
    '''

    def __init__(self, msg):
        self.msg = msg

########################################################################
### Group actions:######################################################
########################################################################

class GroupAction:
	'''
	Class for group actions of a group G on a set M.
	'''
	
	def __init__(self,G,M,action,M_to_hashable=lambda x:x):
		'''
		Here, action(g,m) is a map G x M -> M.
		It may be a left or a right action.
		G and M are any Iterables (of finite cardinality!).
		If elements in M are not hashable, an injective map
		 M_to_hashable from M to some set of hashable objects must be given.
		'''
		self.G = G;
		self.M = M;
		self.M_to_hashable = M_to_hashable;
		self.action = action;
		self.orbitOf_m = None;
		self.orbits = None;
		self.representatives = None;
	
	def Orbits(self):
		'''
		Here, action(g,m) is a map G x M -> M.
		It may be a left or a right action.
		G and M are any Iterables (of finite length!), G should be a group.
		The output is a list of orbits, the orbits being lists of elements of M.
		'''
		if self.orbits != None:
			return self.orbits;

		G = self.G;
		M = self.M;
		action = self.action;
		M_to_hashable = self.M_to_hashable;

		Mcopy = {};
		for m in M:
			Mcopy[M_to_hashable(m)] = m;			
		
		orbits = [];
		while len(Mcopy)>0:
			m_hash,m = Mcopy.popitem();
			orbit = [m];
			for g in G:
				gm_hash = M_to_hashable(action(g,m));
				if gm_hash in Mcopy:
					orbit.append(Mcopy.pop(gm_hash));
			orbit.sort();
			orbits.append(orbit);
		orbits.sort();			
		#orbitOf_m_asLists = {};
		#for m in M:
		#	orbitOf_m_asLists[m] = orbitOf_m[m].list();
		#	orbitOf_m_asLists[m].sort();		
		#return orbitOf_m_asLists;
		self.orbits = orbits;
		orbitOf_m = {};
		for orbit in orbits:
			for m in orbit:
				orbitOf_m[M_to_hashable(m)] = orbit;
		self.orbitOf_m = orbitOf_m;
		return orbits;

	def OrbitRepresentative(self,orbit):
		'''
		Returns a representative of an orbit (namely the first stored one).
		'''
		return orbit[0];

	def OrbitRepresentatives(self):
		'''
		Returns a list of representatives, one for each orbit,
		in the order of the stored orbits.		
		'''
		if self.representatives != None:
			return self.representatives;
			
		orbits = self.Orbits();
		representatives = [self.OrbitRepresentative(orbit) for orbit in orbits];
		self.representatives = representatives;
		return representatives;

	def OrbitOfM(self,m):
		'''
		Returns the orbit that contains m.
		'''
		if self.orbitOf_m == None:
			self.Orbits();
		return self.orbitOf_m[self.M_to_hashable(m)];

	def hashedRepresentativeOfOrbitOfM(self,m):
		'''
		Returns the hashable image of 'the' representative of the
		orbit of m. The representative is the first element in the
		orbit, according to the order in the memory.
		'''
		orbit = self.OrbitOfM(m);
		representative = self.OrbitRepresentative(orbit);
		hashedRepresentative = self.M_to_hashable(representative);
		return hashedRepresentative;

	def GroupElementSendingM1toM2(self,m1,m2):
		'''
		Note: This method may be quite slow!
		Returns a group element that sends m1 to m2.
		'''
		G = self.G;
		M = self.M;
		action = self.action;
		for g in G:
			if m2 == action(g,m1):
				return g;	

########################################################################
### Python helper functions: ###########################################
########################################################################

def removeDuplicates(aList):
	'''
	Takes a list and removes all duplicate entries.
	list(set(aList)) would be a shorter way,
	however this works only when the entries are hashable!
	'''
	result = [];
	for x in aList:
		if x not in result:
			result.append(x);
	return result;

def mod_small(a,p):
	'''
	Computes a representative of a mod p which is closest to 0 as possible.
	'''
	a = a % p;
	if a > p/2:
		a -= p;
	return a;

def aNonsquareModP(ZmodP):
	'''
	Returns a representative of a non-square residue class mod p.
	'''
	for i in range(-1,ZmodP.order()-1):
		if not ZmodP(i).is_square():
			return ZmodP(i);
	#for i in ZmodP:
	#	if not i.is_square():
	#		return i;
	print "Error: There is no non-square mod p =",p,"!";
	return None;

########################################################################
### Galois theory basic functions: #####################################
########################################################################

def lift_automorphism_from_subfield_to_superfield(k,K,L):
	'''
	k is in Gal(K:Q), L/K a Galois extension.
	Returns an element l in Gal(L:Q) whose restriction to K is k.
	'''

	alpha = K.primitive_element();
	k_alpha = k(alpha);
	for l in L.galois_group():
		if l(alpha) == k_alpha:
			return l;
	print "Error: No lift was found!";
	return None;

def restrict_automorphism_from_superfield_to_subfield(l,L,K):
	'''
	l is in Gal(L:Q), L/K a Galois extension.
	Return an element k in Gal(k:Q) which is the restriction of l to K.
	'''
	
	alpha = K.primitive_element();
	l_alpha = l(alpha);
	for k in K.galois_group():
		if k(alpha) == l_alpha:
			return k;
	print "Error: No restriction was found!";
	return None;

########################################################################
### SL_2(Z) and related: ###############################################
########################################################################

def lift_of_Mp_to_ZZ2overP((a1p,a2p),p):
	'''
	Lifting as in section "General Modular Units" of the paper.
	'''
	a1tilda = (a1p.lift() % p)/p;
	a2tilda = (a2p.lift() % p)/p;
	if a2tilda >= 1/2:
		a2tilda -= 1;
	return (a1tilda,a2tilda);

def lift_SL2ZmodP_to_SL2Z(M,p):
	'''
	Lifts a matrix of SL_2(F_p) to SL_2(Z).
	'''
	M_ZZ = M.lift();
	for i in range(M_ZZ.nrows()):
		for j in range(M_ZZ.ncols()):
			M_ZZ[i,j] = mod_small(M_ZZ[i,j],p);
	D, U, V = M_ZZ.smith_form();
	#We have D = U * M_ZZ * V and D is diagonal.
	DetSignChanger = matrix(ZZ,2,2,[-1,0,0,1]);
	if U.det() == -1:
		U = DetSignChanger * U;
		D = DetSignChanger * D;
	if V.det() == -1:
		V = V * DetSignChanger;
		D = D * DetSignChanger;
	a = mod_small(D[0,0],p);
	b = mod_small(D[1,1],p);
	#We replace D=diag(a,b) by the following matrix,
	#which equals D mod p, but which has determinant 1:
	DD = matrix(2,2,[1,-1,1-b,b]) * matrix(2,2,[1,0,1-a,1]) * matrix(2,2,[1,b,0,1]);
	M_ZmodP = U^(-1) * DD * V^(-1);
	return M_ZmodP.change_ring(ZZ);

def init_SL2(field):
	'''
	Returns a list of all matrices in SL_2(field):
	'''
	SL2_field = [];
	for a in field:
		for b in field:
			for c in field:
				if a.is_zero():
					if b*c == field(-1):
						for d_ in field:
							SL2_field.append(matrix(field,2,2,[a,b,c,d_]));
				else:
					d_ = (1+b*c)/a;
					SL2_field.append(matrix(field,2,2,[a,b,c,d_]));
	return SL2_field;

def convert_realQparameter_to_tau(q_rif):
	'''
	q = exp(2*pi*i*tau).
	Given a real interval q_rif for q not containing 0, this method returns
	an "interval" in CC of the corresponding tau values
	in the standard fundamental domain of the upper half plane modulo
	SL_2(Z), unless q_rif intersects [0,1728], in which case the interval
	might lie slightly outside.
	'''
	debug = False;
	#if debug:
	#	q_rif = abs(q_rif); #just for debugging...
	RIFprec = RealIntervalField(q_rif.prec());
	CIFprec = ComplexIntervalField(q_rif.prec());
	if q_rif.contains_zero():
		raise NotImplementedError("q_rif contains zero: This is not properly implemented yet.");
	if q_rif > 0:
		tau1 = -I*RIFprec(log(q_rif)/(2*pi)).upper();
		tau2 = -I*RIFprec(log(q_rif)/(2*pi)).lower();
	else:
		tau1 = 1/2 - I*RIFprec(log(-q_rif)/(2*pi)).upper();
		tau2 = 1/2 - I*RIFprec(log(-q_rif)/(2*pi)).lower();
	if debug:
		print "q_rif =",q_rif;
		print "exp(2*pi*tau1) =",CIFprec(exp(2*pi*I*tau1));
		print "exp(2*pi*tau2) =",CIFprec(exp(2*pi*I*tau2));
	return CIFprec(tau1),CIFprec(tau2);

########################################################################
### Cyclotomic field: ##################################################
########################################################################

def fundamentalUnits_of_realCyclotomicField(p):
	'''
	Consider K = Q(zeta_p + conj(zeta_p)), the real cyclotomic field.
	Here we compute some canonical choice fundamental units, assuming p<=151.
	These elements are considered as elements in Q(zeta_p), NOT as elements in K!
	They are given by |theta_k|, 2 <= k <= (p-1)/2,
	 where |.| is the standard absolute value of CC, and
	 theta_k = (1-zeta_p^k)/(1-zeta_p).
	Their absolute values are computed by multiplying a suitable power of (-zeta_p) to them.
	'''
	
	if p>151:
		print "Only the case p<=151 is implemented, as only there we know that the class number equals one.";
		return None;

	Qzeta = CyclotomicField(p);
	zeta = Qzeta.zeta(); #generator of roots of unity in Qzeta, note: it may have order 2*p!.
	zeta_p = Qzeta.zeta(p); #p'th root of unity
	zetaOrder = Qzeta.zeta_order();
	#theta = [(1-zeta^k)/(1-zeta) for k in range(2,(p-1)/2+1)];

	thetaPm1 = zeta; #-zeta_p;
	arg_thetaPm1 = myCCtoCIF(CC(thetaPm1)).arg();
	#thetaPm1 and all theta's in the following loop generate the units in Qzeta, in case the class number of K is 1.

	#print 1/CC(CC(thetaPm1).arg()/2/pi);

	fundamental_units = [];

	for k in range(2,(p-1)/2+1):
		thetaK = (1-zeta_p^k)/(1-zeta_p);
		arg_thetaK = myCCtoCIF(CC(thetaK)).arg();
		#print 1/CC(CC(theta).arg()/2/pi);
		for i in range(Qzeta.zeta_order()):
			arg_possible_eta = RIF((arg_thetaK + i*arg_thetaPm1)/2/pi);
			try:
				arg_possible_eta = arg_possible_eta - arg_possible_eta.unique_round();
			except ValueError:
				continue;
			if arg_possible_eta.abs() < 1/2/zetaOrder:
				possible_eta = thetaK * thetaPm1^i;
				if possible_eta.is_real_positive():
					fundamental_units.append(possible_eta);
					break;
		else:
			#print "p =",p;
			print "Error: didn't find the absolute value of theta_k due to insufficient precision.";
			raise PrecisionError("Error: didn't find the absolute value of theta_k due to insufficient precision.");
	#print fundamental_units;
	return fundamental_units;

#fundamentalUnits_of_realCyclotomicField(71);

########################################################################
### Modular curve X_ns^+(p): ###########################################
########################################################################

class X_ns_plus:
	'''
	The moduli curve X_{ns}^+(p).
	'''
	
	def __init__(self,p,d=None):
		'''
		p is a prime >= 7, which determines the moduli curve X_{ns}^+(p).
		d >= 2 is a divisor of (p-1)/2, which determines the subfield
		 of degree d of the real cyclotomic field of order p that we
		 will work with.		
		'''
		
		self.p = p;

		global precision;
		self.stdPrecision = precision;
		
		if d==None:
			#Choose a divisor d of (p-1)/2, such that
			# - it is not too large (we need to be able to compute fundamental units in K)
			# - it is not too small (otherwise we get too many wrong candidates).
			print "(p-1)/2 =",factor((p-1)/2);
			d = ZZ((p-1)/2);
			#d = ((p-1)//2).prime_divisors()[0];
		else:
			#Check whether given d is valid:
			if d<=1 or ((p-1)/2) % d != 0:
				raise ValueError("Error: The input d is not valid!");
		print "d = deg(K) =",d;
		self.d = d;

		#Ring Z/pZ:
		ZmodP = IntegerModRing(p);
		self.ZmodP = ZmodP;
		
		#Space of 2x2-matrices over Z/pZ:
		Matrices_ZmodP = MatrixSpace(ZmodP,2,2);
		self.Matrices_ZmodP = Matrices_ZmodP;
		
		#Set of invertable 2x2-matrices over Z/pZ:
		#GL2_ZmodP = [M for M in Matrices_ZmodP if det(M) != 0];
		#Set of 2x2-matrices over Z/pZ with determinant 1:
		SL2_ZmodP = init_SL2(ZmodP);
		self.SL2_ZmodP = SL2_ZmodP;
		print "|SL_2(Z/pZ)| =",len(SL2_ZmodP);
		
		#Cyclotomic field of order p:
		Qzeta = CyclotomicField(p);
		self.Qzeta = Qzeta;
		print "Q(zeta_p) =",Qzeta;
		zetaP = Qzeta.zeta(p); #We take a p'th root of unity, but not necessarily a generator of the roots of unity in Qzeta!
		self.zetaP = zetaP;
		print "zeta_p =",zetaP,"(its order is",zetaP.order(),")";
		
		#The Galois group of Qzeta (as permutation group):
		GalQzeta = Qzeta.galois_group();
		self.GalQzeta = GalQzeta;

		#The same, now as elements in (Z/pZ)^*
		GalQzeta_subsetOfZmodP = [self.elementOfGalQzeta_AsElementIn_ZmodPstar(g) for g in GalQzeta];
		self.GalQzeta_subsetOfZmodP = GalQzeta_subsetOfZmodP;
		print "GalQzeta as a subset of Z/pZ =",GalQzeta_subsetOfZmodP;
		GalQZeta_gen = GalQzeta.gen(0);
		self.GalQZeta_gen = GalQZeta_gen;
		
		#A subgroup of the Galois group GalQzeta of index d, in different representations:
		H, H_subsetOfGalQzeta, H_subsetOfZmodP = self.init_subgroupH_of_GalQzeta();
		self.H = H;
		self.H_subsetOfGalQzeta = H_subsetOfGalQzeta;
		self.H_subsetOfZmodP = H_subsetOfZmodP;

		#for h in H:
		#	print "Element h in H:",h;
		print "Chosen subgroup H of the Galois group of Qzeta has order:",H.order();
		print "H as a subset of Z/pZ =",H_subsetOfZmodP;
		#The fixed field of H, and its inclusion into Qzeta:
		print "From version 04 on we don't compute K anymore.";
		#K, incl_K_to_Qzeta = H.fixed_field();
		#self.K = K;
		#self.incl_K_to_Qzeta = incl_K_to_Qzeta;		
		#print "K =",K;

		#The Galois group of K:
		print "And GalK will be identified with some choice of lift in GalQzeta!"; 
		#We don't compute GalK directly via K anymore,
		# as this is implemented in Sage only up to d<=11:
		#GalK = K.galois_group(); 
		#self.GalK = GalK;
		#Take instead a set of representatives of GalQZeta/H:
		GalK_lift = [GalQZeta_gen^i for i in range(d)];
		self.GalK_lift = GalK_lift;
		
		#Thus we have an extension 0 -> H -> GalQzeta -> GalK -> 0
		#The number m from section "Units" in the paper:
		if (p+1)*H.order() % 3 == 0:
			m = 2;
		else:
			m = 6;
		self.m = m;
		print "m =",m;
		
		#The fundamental units of K:
		if len(H) != 2:
			print "TODO: Improve the algorithm according to Yuri's idea that circumvents the need of computing the fundamental units.";
			K, incl_K_to_Qzeta = H.fixed_field();
			#The number \eta_0 from section "The Principal Relation":
			eta0 = (1-zetaP).norm(K); #The relative norm of zeta_p in K.
			print "eta0 =",eta0;
			fundamental_units_of_K = list(K.units());
			fundamental_units_of_K = [incl_K_to_Qzeta(u) for u in fundamental_units_of_K];
			print "fundamental units =",fundamental_units_of_K;
			#print "TODO: We should improve the basis of fundamental units!"; #Update: As we use ellipsoids, this is not really necessary anymore, except for the original first reduction, which does not use ellipsoids.
			#List of eta0 ... eta(d-1) as in section "The Principal Relation" in the paper:
		else:
			#|H|=2, i.e. K = Q(zeta_p+conj(zeta_p)) is the real cyclotomic field.
			#Thus we can take ready fundamental units if p<=151:
			fundamental_units_of_K = fundamentalUnits_of_realCyclotomicField(p);

			#The number \eta_0 from section "The Principal Relation":
			eta0 = prod([h(1-zetaP) for h in H_subsetOfGalQzeta]); #The relative norm of zeta_p in K.
			print "eta0 =",eta0;
			print "eta0 is real positive:",eta0.is_real_positive();

		etas = [eta0] + fundamental_units_of_K;
		self.etas = etas;
		print "#etas =",len(etas);

		#F_p^2 \ (0,0):
		Mp = [(ZmodP(i),ZmodP(j)) for i in range(p) for j in range(p) if (i,j)!=(0,0)];
		self.Mp = Mp;
		
		#Lift of F_p^2 \ (0,0) to p^(-1)Z^2:
		MpLift = [lift_of_Mp_to_ZZ2overP(a,p) for a in Mp];
		self.MpLift = MpLift;
		#print MpLift;

		#A representative of a nonsquare residue class mod p:
		Xi = aNonsquareModP(ZmodP);
		self.Xi = Xi;
		print "Xi =",Xi,"(nonsquare mod p)";

		#A normalizer of a non-split Cartan subgroup of GL_2(F_p):
		G = [Matrices_ZmodP([a,Xi*b,b,a]) for (a,b) in Mp];
		G += [Matrices_ZmodP([a,Xi*b,-b,-a]) for (a,b) in Mp];
		self.G = G;

		#The elements of G with determinant 1:
		G1 = [g for g in G if det(g)==ZmodP(1)];
		self.G1 = G1;
		
		#The elements of G with determinant +-1:
		Gpm1 = [g for g in G if det(g) in [ZmodP(1),ZmodP(-1)]];
		self.Gpm1 = Gpm1;
		
		#The elements of G with determinant in H:
		GH = [g for g in G if det(g) in H_subsetOfZmodP];
		self.GH = GH;
		
		print "|G_1| =",len(G1);
		print "|G_pm1| =",len(Gpm1);
		print "|G_H| =",len(GH);

		leftAction_Gpm1_on_Mp = GroupAction(Gpm1,Mp,lambda g,m: tuple(g*vector(m)));
		self.leftAction_Gpm1_on_Mp = leftAction_Gpm1_on_Mp;

		rightAction_GH_on_Mp = GroupAction(GH,Mp,lambda g,m: tuple(vector(m)*g));
		self.rightAction_GH_on_Mp = rightAction_GH_on_Mp;
		
		orbits_Mp_mod_GH = rightAction_GH_on_Mp.Orbits();
		self.orbits_Mp_mod_GH = orbits_Mp_mod_GH;
		Orbit0 = orbits_Mp_mod_GH[0];
		self.Orbit0 = Orbit0;
		print "Length of first orbits O of Mp/G_H:",len(Orbit0);
		#reps_Mp_mod_GH = rightAction_GH_on_Mp.OrbitRepresentatives();
		#print "Representatives of Mp/G_H:",reps_Mp_mod_GH;

		rightAction_H_on_GalQzeta = GroupAction(H_subsetOfZmodP,GalQzeta_subsetOfZmodP,lambda h,g: g * h);
		self.rightAction_H_on_GalQzeta = rightAction_H_on_GalQzeta;
		
		#Let G1 act on SL_2(Z/pZ) from the left:
		leftAction_G1_on_SL2ZmodP = GroupAction(G1,SL2_ZmodP,lambda g,M: g*M,M_to_hashable=str);
		self.leftAction_G1_on_SL2ZmodP = leftAction_G1_on_SL2ZmodP;

		#Let G1 act on Mp = F_p^2 \ (0,0): from the left:
		leftAction_G1_on_Mp = GroupAction(G1,Mp,lambda g,m: tuple(g*vector(m)));
		self.leftAction_G1_on_Mp = leftAction_G1_on_Mp;

		#Representatives of cusps as chosen in section "Cusps" in the paper:
		cuspRepresentatives = self.init_cuspRepresentatives();
		self.cuspRepresentatives = cuspRepresentatives;
		print "Cusp representatives:",cuspRepresentatives;

		#Compute optimal system of representatives, see section "Cusps" (and section 2.1 of version 1 of the arXiv paper):
		Sigma = self.init_optimalSystemOfRepresentatives();
		self.Sigma = Sigma;
		print "|Sigma| =",len(Sigma);

		#Eta: d times d matrix whose entries are the logarithms of the
		# absolute values of the Galois conjugates of the etas,
		# computed up to the given precision, and
		#Alpha is Eta's inverse matrix:
		self.Eta = {};
		self.Alpha = {};
		self.init_Eta_and_Alpha(prec=self.stdPrecision);
		print "Eta[0,0] =",self.Eta[self.stdPrecision][0,0];
		print "det(Eta) =",det(self.Eta[self.stdPrecision]);
		print "Alpha[0,0] =",self.Alpha[self.stdPrecision][0,0];

		ZmodPstar_to_GalKIndex = self.init_ZmodPstar_to_GalKIndex();
		self.ZmodPstar_to_GalKIndex = ZmodPstar_to_GalKIndex;
		print "ZmodPstar_to_GalKIndex:",ZmodPstar_to_GalKIndex;

		lift_GalKIndex_to_G = self.init_lift_GalKIndex_to_G();
		self.lift_GalKIndex_to_G = lift_GalKIndex_to_G;
		print "lift_GalKIndex_to_G[0]:",lift_GalKIndex_to_G[0].list();


		self.exp_2pi_I_a2 = {};
		self.gamma_a = {};
		self.log_abs_gamma_a = {};
		self.ell_a = {};
		self.gamma_c_k_by_OrbitRepresentative = {};	
		self.log_abs_gamma_c_k_by_OrbitRepresentative = {};	
		self.ord_c_Uk_by_OrbitRepresentative = {};	
		self.kappa = {};
		self.Theta = {};
		self.deltaMax = {};
		self.thetaMax = {};
		self.B0_old1 = {}; #from first arxiv version
		self.B0 = {};
		self.init_variables_around_bakerBound(prec=self.stdPrecision);

		self.init_expected_integral_points();

	def init_expected_integral_points(self):
		'''
		Computes all (non-cusp) points of the moduli curve that
		are expected to be integral.
		'''

		CIF = ComplexIntervalField(self.stdPrecision);
		
		CM = cm_j_invariants_and_orders(QQ);
		for (discriminant,conductor,j) in CM:
			D = Mod(discriminant,self.p);
			if D.is_zero():
				continue;
			if D.is_square():
				continue;
			print " ================================ ";
			print "j =",j,"(discr ="+str(discriminant)+") should give a point on X_ns^+(p)!";

			E = EllipticCurve_from_j(j);
			tau = E.period_lattice().tau(100);
			q = CIF(exp(2*pi*I*tau));
			print "E =",E;
			print "tau =",tau;
			print "q =",q;
			
		#raise Exception("Debug stop");

		return;
	
	def init_cuspRepresentatives(self):
		'''
		Returns representatives of cusps as chosen in section "Cusps" of the paper.
		'''
		
		ZmodP = self.ZmodP;
		Xi = self.Xi;

		#First construct a dictionary, such that squareRootsModP[x] is a list of all square roots of x mod p:
		p = ZmodP.order();
		squareRootsModP = {};
		for x in ZmodP:
			squareRootsModP[x] = [];
		for x in ZmodP:
			squareRootsModP[x^2].append(x);
		print "Square roots mod p:",squareRootsModP;

		#Representatives of cusps as chosen in section "Cusps" of the paper:
		cuspRepresentatives = [];
		for c in range(1,(p-1)/2+1):
			for b in ZmodP:
				aSq = 1/ZmodP(c) + Xi*b^2;
				if squareRootsModP[aSq] != []:
					a = squareRootsModP[aSq][0];
					g = matrix(ZmodP,2,2,[c*a,b*Xi,c*b,a]);
					gLift = lift_SL2ZmodP_to_SL2Z(g,p);
					cusp = tuple(gLift * vector([1,0]));
					cuspRepresentatives.append(cusp);
					break;
		return cuspRepresentatives;

	def init_optimalSystemOfRepresentatives(self,verbose=False):
		'''
		Computes an optimal system of representatives
		for the orbits of the left group action of G1 on SL_2(Z).
		This optimality is explained in Sect. 2.1 of the first arXiv version of the paper.
		It is in fact not needed for our purposes, any representatives will do.
		'''
		
		p = self.p;
		ZmodP = self.ZmodP;

		#Compute the orbits of left action of G1 on SL_2(Z/pZ):
		SigmaModP = self.leftAction_G1_on_SL2ZmodP.OrbitRepresentatives();

		cuspCorrespondingToOrbitByItsFirstElement = {};
		for cusp in self.cuspRepresentatives:
			cuspModP = (ZmodP(cusp[0]),ZmodP(cusp[1]));
			orbit = self.leftAction_G1_on_Mp.OrbitOfM(cuspModP);
			if cuspCorrespondingToOrbitByItsFirstElement.has_key(orbit[0]):
				raise Exception("Error: cusp mod p does not uniquely determine the cusp! That is, the algorithm must work harder as it does right now!");
			cuspCorrespondingToOrbitByItsFirstElement[orbit[0]] = cusp;
		print "cuspCorrespondingToOrbitByItsFirstElement =", cuspCorrespondingToOrbitByItsFirstElement;
		
		#Make Sigma an optimal system of representatives:
		Sigma = [];
		for sigmaModP in SigmaModP:
			a = sigmaModP[0,0];
			c = sigmaModP[1,0];
			m = (a,c);
			if verbose:
				print "----"
				print "sigmaModP =",sigmaModP.list();
				print "m =",m;
			orbit = self.leftAction_G1_on_Mp.OrbitOfM(m);
			cusp = cuspCorrespondingToOrbitByItsFirstElement[orbit[0]];
			if verbose:
				print "cusp =",cusp;
			cuspModP = (ZmodP(cusp[0]),ZmodP(cusp[1]));
			g = self.leftAction_G1_on_Mp.GroupElementSendingM1toM2(m,cuspModP);
			#gamma = lift_SL2ZmodP_to_SL2Z(g,p);
			#sigmaNonOptimal = lift_SL2ZmodP_to_SL2Z(sigmaModP,p);
			#sigmaAlmostOptimal = gamma * sigmaNonOptimal;
			sigmaAlmostOptimal = lift_SL2ZmodP_to_SL2Z(g*sigmaModP,p);
			if verbose:
				print "sigmaAlmostOptimal =",sigmaAlmostOptimal.list();
			#We want to multiply some gamma in the lift of G1 to sigmaAlmostOptimal, such
			#that it sends Infinity = (1,0) to cusp.
			#So far its image of Infinity differs by cusp only by a multiple of p.
			#Thus gamma reduces to the identity matrix mod p:
			#gamma = [[p*kx+1,p*ky],[p*kz,p*kw+1]]
			#With below notation, we thus need a gamma such that
			# gamma * (a,c)^t = (aStar,cStar)^t
			#This implies that kx,ky,kz,kw can be written in terms of k1 and k2 as below.
			#Moreover, det(gamma)=1, which is equivalent to
			# f1*k1 + f2*k2 = k, with f1,f2,k as below.
			#With xgcd we find h1,h2,g such that f1*h1 + f2*h2 = g,
			#so we can define k1 and k2 as below in terms of h1,h2,g,k.
			#(Note that all possible choices of k1 and k2 are parametrized affinely
			# over some parameter in Z.)

			a = sigmaAlmostOptimal[0,0];
			c = sigmaAlmostOptimal[1,0];
			tmp, s, t = xgcd(a,c);
			if tmp == -1:
				s = -s;
				t = -t;
			aStar = cusp[0];
			cStar = cusp[1];
			aTilda = ZZ((aStar-a)/p); 
			cTilda = ZZ((cStar-c)/p);
			
			f1 = -(a*cTilda*p*s + c*cTilda*p*t + c);
			f2 = a*aTilda*p*s + aTilda*c*p*t + a;
			k = aTilda*s + cTilda*t;
			g, h1, h2 = xgcd(f1,f2);
			k1 = h1 * ZZ(k/g);
			k2 = h2 * ZZ(k/g);
			kx = s*aTilda + c*k1;
			ky = t*aTilda - a*k1;
			kz = s*cTilda + c*k2;
			kw = t*cTilda - a*k2;
			gamma = matrix(2,2,[p*kx+1,p*ky,p*kz,p*kw+1]);
			sigmaOptimal = gamma*sigmaAlmostOptimal;
			if verbose:
				print "sigmaOptimal =",sigmaOptimal.list();
			Sigma.append(sigmaOptimal);		
			continue;

		if verbose:
			print "|Sigma| =",len(Sigma);
		return Sigma;

	def init_subgroupH_of_GalQzeta(self):
		'''
		Returns the subgroup of the galois group of Q(zeta_p) of index d:
		Here, Gal is the cyclic Galois group of Q(zeta_p).
		d is the index of H in this Galois group.
		'''
		d = self.d;
		g = self.GalQzeta.gen(0);
		n = g.order();
		gd = g^d;
		g_ZmodP = self.elementOfGalQzeta_AsElementIn_ZmodPstar(g);
		H_subsetOfGalQzeta = [gd^i for i in range(g.order()/d)];
		H_subsetOfZmodP = [g_ZmodP^(d*i) for i in range(g.order()/d)];
		H = self.GalQzeta.subgroup(H_subsetOfGalQzeta);
		return H, H_subsetOfGalQzeta, H_subsetOfZmodP;

	def init_ZmodPstar_to_GalKIndex(self):
		'''
		There is a canonical map from F_p^* --> F_p^*/H --> Gal(K/Q), the latter being an iso.
		The output is a dictionary which represents this map.
		'''
		ZmodPstar_to_GalKIndex = {};
		for kIndex in range(self.d):
			#k = self.GalK[kIndex];
			k = self.GalK_lift[kIndex];
			#g = lift_automorphism_from_subfield_to_superfield(k,self.K,self.Qzeta);
			g = k;
			i0 = self.elementOfGalQzeta_AsElementIn_ZmodPstar(g);
			for h in self.H_subsetOfZmodP:
				i = i0 * h;
				ZmodPstar_to_GalKIndex[i] = kIndex;
		return ZmodPstar_to_GalKIndex;

	def init_lift_GalKIndex_to_G(self):
		lift_GalKIndex_to_G = {};
		for g in self.G:
			kIndex = self.GL2ZmodP_to_GalKIndex(g);
			lift_GalKIndex_to_G[kIndex] = g;
		return lift_GalKIndex_to_G;		

	def GL2ZmodP_to_GalKIndex(self,M):
		#print "type of M:",type(M), "base ring:",M.base_ring(),"det(M) =",det(M);
		return self.ZmodPstar_to_GalKIndex[det(M)];

	def elementOfGalQzeta_AsElementIn_ZmodPstar(self,g):
		'''
		We assume here that p is an odd prime!
		Gal is the Galois group of Q(zeta_p).
		We take its generator and check how it acts on zeta_p:
		The image of zeta_p under this auto must be zeta_p^i for some i.
		We return ZmodP(i).
		'''
		zeta = self.Qzeta.zeta();
		#print "zeta =",zeta;
		#print "g =",g,type(g);
		gZeta = g(zeta);
		#print "gZeta =",gZeta;
		for i in range(self.Qzeta.zeta_order()):
			if gZeta == zeta^i:
				return self.ZmodP(i);
		print "Error in generatorOfGal_AsElementIn_ZmodPstar().";

	def j_invariant_in_qInterval(self,q_rif):
		'''
		q is assumed to be a real interval!
		The output is a real interval containing all associated j-invariants (possibly more).
		'''
		debug = False;
		RIFprec = RealIntervalField(q_rif.precision());
		CIFprec = ComplexIntervalField(q_rif.precision());
		CCprec = ComplexField(q_rif.precision());
		if q_rif.contains_zero():
			print "q contains zero: This is not properly implemented yet, so we output the trivial interval RR...";
			return RIFprec(-Infinity,Infinity);
		if q_rif > 0:
			tau1 = -I*RIFprec(log(q_rif)/(2*pi)).upper();
			tau2 = -I*RIFprec(log(q_rif)/(2*pi)).lower();
		else:
			tau1 = 1/2 - I*RIFprec(log(-q_rif)/(2*pi)).upper();
			tau2 = 1/2 - I*RIFprec(log(-q_rif)/(2*pi)).lower();
		if debug:
			print "tau1 =",tau1;
			print "tau2 =",tau2;
		j1 = myRRtoRIF(elliptic_j(CCprec(tau1)).real()); #taking CCprec here is important, because otherwise tau gets treated like a symbolic expression and gets rounded to CC.
		j2 = myRRtoRIF(elliptic_j(CCprec(tau2)).real());
		if debug:
			print "j1 =",j1;
			print "j2 =",j2;
		jMin = min(j1,j2);
		jMax = max(j1,j2);
		j = RIFprec(jMin.lower(),jMax.upper());
		return j;

	def init_Eta_and_Alpha(self,prec=None):
		'''
		d times d matrix whose entries are the logarithms of the
		 absolute values of the Galois conjugates of the etas,
		 computed up to the given precision:
		'''
		d = self.d;
		GalK_lift = self.GalK_lift;
		etas = self.etas;
		Qzeta = self.Qzeta;
		if prec == None:
			prec = self.stdPrecision;
		CCprec = ComplexField(prec);
		CIFprec = ComplexIntervalField(prec);
		RIFprec = RealIntervalField(prec);

		precEta = prec+self.stdPrecision;

		firstLoop = True;

		while True:
			CCprecEta = ComplexField(precEta);
			CIFprecEta = ComplexIntervalField(precEta);
			RIFprecEta = RealIntervalField(precEta);

			Eta = matrix(RIFprecEta,d,d);
			for kIndex in range(d):
				#sigma = GalK[kIndex];
				#sigmaLift = lift_automorphism_from_subfield_to_superfield(sigma,K,Qzeta);
				sigmaLift = GalK_lift[kIndex];
				for j in range(d):
					eta = etas[j];
					#etaSigma = sigma(eta);
					etaSigma = sigmaLift(Qzeta(eta));
					
					#etaSigma = CIFprec(etaSigma);
					etaSigma = myCCtoCIF(ComplexField(precEta + 4)(etaSigma));
					Eta[kIndex,j] = log(abs(etaSigma));
				Eta = Eta.change_ring(RIFprecEta);
			if firstLoop:
				self.Eta[prec] = Eta.change_ring(RIFprec);			

			#The inverse of the matrix Eta:
			Alpha = Eta.inverse(); #Will always yield a result in CIF, but maybe involving +-infinity entries.

			maxDiameter = max([x.diameter() for x in Alpha.list()]);

			if maxDiameter > 2^(-prec+1):
				precEta += 200;
				print "maxDiameter for Alpha:",maxDiameter,"is too big, so restart computing Eta with precision",precEta;
				firstLoop = False;
				continue;

			self.Alpha[prec] = Alpha.change_ring(RIFprec);

			#TODO: Should increase the precision for Eta until Alpha has precise enough entries!
			#Otherwise all later sieves will be much slower!

			if True:
				print "precise Eta[0,0] =",Eta[0,0];
				print "precise Alpha[0,0] =",Alpha[0,0];
				print "self.Eta[0,0] =",self.Eta[prec][0,0];
				print "self.Alpha[0,0] =",self.Alpha[prec][0,0];
				#raise Exception("Debug");
			break;

	def	ord_c_Uk(self,cIndex,kIndex):
		'''
		Computes ord_c of U^k according to above Proposition 6.6:
		'''
		ell_a = self.ell_a;
		sigma_c = self.Sigma[cIndex];
		sigma_k = self.lift_GalKIndex_to_G[kIndex];
		ord_c_Uk = 0;
		for a0_in_Mp in self.Orbit0:
			a_in_Mp = self.rightAction_GH_on_Mp.action(sigma_k*sigma_c,a0_in_Mp);
			ord_c_Uk += self.p*self.m*ell_a[a_in_Mp];
		return ord_c_Uk;

	def	gamma_c_k(self,cIndex,kIndex,prec=None):
		'''
		Computes gamma_{c,k} according to section "Approximate Formulas" of the paper.
		'''
		if prec == None:
			prec = self.stdPrecision;
		CIFprec = ComplexIntervalField(prec);
		gamma_a = self.gamma_a[prec];
		sigma_c = self.Sigma[cIndex];
		sigma_k = self.lift_GalKIndex_to_G[kIndex];
		gamma_c_k = CIFprec(1);
		for a0_in_Mp in self.Orbit0:
			a_in_Mp = self.rightAction_GH_on_Mp.action(sigma_k*sigma_c,a0_in_Mp);
			gamma_c_k *= CIFprec(gamma_a[a_in_Mp]);
		gamma_c_k = gamma_c_k^self.m;
		return gamma_c_k;

	def log_abs_gamma_c_k(self,cIndex,kIndex,prec=None):
		'''
		Computes log |gamma_{c,k}| according to section "Approximate Formulas" of the paper.
		'''
		if prec == None:
			prec = self.stdPrecision;
		RIFprec = RealIntervalField(prec);
		log_abs_gamma_a = self.log_abs_gamma_a[prec];
		sigma_c = self.Sigma[cIndex];
		sigma_k = self.lift_GalKIndex_to_G[kIndex];
		log_abs_gamma_c_k = RIFprec(0);
		for a0_in_Mp in self.Orbit0:
			a_in_Mp = self.rightAction_GH_on_Mp.action(sigma_k*sigma_c,a0_in_Mp);
			log_abs_gamma_c_k += RIFprec(log_abs_gamma_a[a_in_Mp]);
		log_abs_gamma_c_k *= self.m;
		return log_abs_gamma_c_k;		

	def	DEPRECIATED_ord_c_Uk(self,cIndex,kIndex):
		'''
		SLOW!!! (if used a lot)
		Computes the order of U^(sigma_k) at c according to the
		 formula before Proposition 6.6.
		NOTE: This can be optimized!
		'''
		result = 0;
		p = self.p;
		ell_a = self.ell_a;
		sigma_c = self.Sigma[cIndex];
		sigma_k = self.lift_GalKIndex_to_G[kIndex];
		for a0_in_Mp in self.Orbit0:
			#a_in_Mp = self.rightAction_GH_on_Mp.action(sigma_c*sigma_k,a0_in_Mp); #Wrong order!!!
			a_in_Mp = self.rightAction_GH_on_Mp.action(sigma_k*sigma_c,a0_in_Mp);
			result += ell_a[a_in_Mp];
		result *= p*self.m;
		return result;
		
	def DEPRECIATED_gamma_c_k(self,cIndex,kIndex):
		'''
		SLOW!!! (if used a lot)
		Computes gamma_{c,k} according to section "Approximate Formulas" of the paper.
		'''
		result = CIF(1);
		p = self.p;
		gamma_a = self.gamma_a;
		sigma_c = self.Sigma[cIndex];
		sigma_k = self.lift_GalKIndex_to_G[kIndex];
		for a0_in_Mp in self.Orbit0:
			#a_in_Mp = self.rightAction_GH_on_Mp.action(sigma_c*sigma_k,a0_in_Mp); #Wrong order!!!
			a_in_Mp = self.rightAction_GH_on_Mp.action(sigma_k*sigma_c,a0_in_Mp);
			result *= CIF(gamma_a[a_in_Mp]); 
		result = result^self.m;
		return result;

	def DEPRECIATED_log_abs_gamma_c_k(self,cIndex,kIndex):
		'''
		SLOW!!! (if used a lot)
		Computes log |gamma_{c,k}| to section "Approximate Formulas" of the paper.
		'''
		result = RIF(0);
		p = self.p;
		log_abs_gamma_a = self.log_abs_gamma_a;
		sigma_c = self.Sigma[cIndex];
		sigma_k = self.lift_GalKIndex_to_G[kIndex];
		for a0_in_Mp in self.Orbit0:
			#a_in_Mp = self.rightAction_GH_on_Mp.action(sigma_c*sigma_k,a0_in_Mp); #Wrong order!!!
			a_in_Mp = self.rightAction_GH_on_Mp.action(sigma_k*sigma_c,a0_in_Mp);
			result += RIF(log_abs_gamma_a[a_in_Mp]);
		result *= self.m;
		return result;

	def delta_c_j(self,cIndex,jIndex,prec=None):
		'''
		Computes delta_{c,j} according to section "Approximate Formulas" of the paper.
		Does NOT depend only on the cooresponding orbit!
		Therefore, storing all delta_c_j in memory would take a lot of space!
		'''
		if prec == None:
			prec = self.stdPrecision;
		RIFprec = RealIntervalField(prec);
		result = RIFprec(0);
		Alpha = self.Alpha[prec];		
		for kIndex in range(self.d):
			result += Alpha[kIndex,jIndex] * self.ord_c_Uk(cIndex,kIndex);			
		result *= RIFprec(-1/self.p);
		return result;

	def theta_c_j(self,cIndex,jIndex,prec=None):
		'''
		Computes theta_{c,j} according to section "Approximate Formulas" of the paper.
		Does NOT depend only on the cooresponding orbit!
		Therefore, storing all theta_c_j in memory would take a lot of space!
		'''
		if prec == None:
			prec = self.stdPrecision;
		RIFprec = RealIntervalField(prec);
		result = RIFprec(0);
		Alpha = self.Alpha[prec];		
		for kIndex in range(self.d):
			result += Alpha[kIndex,jIndex] * self.log_abs_gamma_c_k(cIndex,kIndex,prec);			
		return result;

	def init_variables_around_bakerBound(self,prec,verbose=False):
		#verbose = True;
		
		if prec == None:
			prec = self.stdPrecision;
		RIFprec = RealIntervalField(prec);
		CIFprec = ComplexIntervalField(prec);
		Mp = self.Mp;
		p = self.p;
		d = self.d;
		Sigma = self.Sigma;
		Alpha = self.Alpha[prec];

		exp_2pi_I_a2 = {};
		for j in range(-p,p+1):
			exp_2pi_I_a2[j/p] = CIFprec((2*CIFprec(pi)*I*j/p).exp());
		self.exp_2pi_I_a2[prec] = exp_2pi_I_a2;

		#Initialize gamma[a], log_abs_gamma[a], and ell[a]:
		gamma_a = {};
		log_abs_gamma_a = {};
		ell_a = {};
		for a_in_Mp in Mp:
			a1,a2 = lift_of_Mp_to_ZZ2overP(a_in_Mp,p);
			#According to section "Approximate Formulas" of the paper:
			if a1 != 0:
				gamma_a[a_in_Mp] = CIFprec(exp(CIFprec(pi)*I*a2*(a1-1)));
			else:
				gamma_a[a_in_Mp] = CIFprec(exp(CIFprec(pi)*I*a2*(a1-1))*(1-exp_2pi_I_a2[a2]));
			log_abs_gamma_a[a_in_Mp] = log(abs(gamma_a[a_in_Mp]));
			#print "debug234",type(log_abs_gamma_a[a_in_Mp]);
			#raise Exception();
			ell_a[a_in_Mp] = bernoulli_polynomial(a1,2)/2;
		self.gamma_a[prec] = gamma_a;
		self.log_abs_gamma_a[prec] = log_abs_gamma_a;
		self.ell_a = ell_a; #ell_a does not depend on the prec.

		#self.init_gammaCk_and_logAbsGammaCk_and_ordcUk_by_OrbitRepresentative(prec=prec);
		#print "debug4",self.gamma_c_k_by_OrbitRepresentative.keys();
		#print "debug5",self.log_abs_gamma_c_k_by_OrbitRepresentative.keys();

		#gamma_c_k_by_OrbitRepresentative = self.gamma_c_k_by_OrbitRepresentative[prec];
		#log_abs_gamma_c_k_by_OrbitRepresentative = self.log_abs_gamma_c_k_by_OrbitRepresentative[prec];
		#ord_c_Uk_by_OrbitRepresentative = self.ord_c_Uk_by_OrbitRepresentative;

		kappa = max([sum([abs(Alpha[kIndex,jIndex]) for jIndex in range(d)]) for kIndex in range(d)]);
		self.kappa[prec] = kappa;

		Theta = kappa*self.m*(self.p+1)*len(self.H);
		self.Theta[prec] = Theta;

		#More precise, but it takes forever:
		#deltaMax = max([self.delta_c_j(cIndex,jIndex).abs() for cIndex in range(len(Sigma)) for jIndex in range(d)]);
		#thetaMax = max([self.theta_c_j(cIndex,jIndex).abs() for cIndex in range(len(Sigma)) for jIndex in range(d)]);
		#Quick and dirty:
		representatives = self.rightAction_GH_on_Mp.OrbitRepresentatives();

		#max_gamma_c_k = max([gamma_c_k_by_OrbitRepresentative[rep].abs() for rep in representatives]);
		#max_log_abs_gamma_c_k = max([log_abs_gamma_c_k_by_OrbitRepresentative[rep].abs() for rep in representatives]);
		#max_ord_c_Uk = max([ord_c_Uk_by_OrbitRepresentative[rep].abs() for rep in representatives]);

		max_gamma_c_k = max(gamma_a.values()) ^ (self.m * len(self.Orbit0));
		max_log_abs_gamma_c_k = len(self.Orbit0) * self.m * max(log_abs_gamma_a.values());
		max_ord_c_Uk = len(self.Orbit0) * p * self.m * max(ell_a.values());
		


		AlphaMax = max([Alpha[i,j].abs() for i in range(d) for j in range(d)]);
		deltaMax = RIFprec(1/p * kappa * max_ord_c_Uk).abs();
		thetaMax = RIFprec(kappa * max_log_abs_gamma_c_k).abs();
		self.deltaMax[prec] = deltaMax;
		self.thetaMax[prec] = thetaMax;

		U1 = RIFprec(10^8*deltaMax*9^d*d^6*p^(4*d+2));
		U1 *= prod([myRRtoRIF(self.etas[k].global_height().abs()) for k in range(1,d)]);
		U2 = RIFprec(U1 + thetaMax + 4*kappa*p^3);
		B0_old1 = RIFprec(2*U1*log(U1) + 2*U2); #From first arxiv version.
		self.B0_old1[prec] = B0_old1;
		if verbose:
			print "Baker bound B0_old1 =",B0_old1;
		
		#	
		#New version (second arxiv version) of Baker's bound:
		#
		delta = RIFprec(min([n for n in divisors(ZZ((p-1)/2)) if n>=3]));
		#print "delta:",delta;
		#Ohm is an upper bound for j(P), P an integral point:
		Ohm = RIFprec(30^(delta+5)*delta^(-2*delta+RIFprec(4.5))*p^(6*delta+5)*(log(RIFprec(p)))^2);
		#print "Ohm:",Ohm;
		#Ohm0 is an upper bound for |log(1/q_c(P))|:
		#Ohm0 = RIFprec(log(exp(Ohm)+2079)); #This is the definition for Ohm0 in the paper, however it needs too high precision this way!
		Ohm0 = max(Ohm,1) + 2079; #Take instead this simple estimate from above. 
		#print "Ohm0:",Ohm0;
		#B0 is an upper bound for the exponents b_k in the fundamental relation:
		B0 = deltaMax*Ohm0 + thetaMax + Theta;
		self.B0[prec] = B0;
		if verbose:
			print "New Baker bound B0 =",B0;
		#raise Exception(); #Debug

	def reduction_of_BakerBound_B0_and_qc_at_sigmaC(self,sigma_c,B0=None,j1=0,j2=1,verbose=False):
		'''
		Notation is in large parts as in the paper.
		One difference is that "k" in the paper is here a "j",
		 and the "l" in the paper is here a "k" (because it denotes an element in K).
		'''
		debug = False;
		if debug:
			verbose = True;
		
		if verbose:
			print "Reduction of Baker bound:";
			print "j1 =",j1,"and j2 =",j2;

		RIFstd = RealIntervalField(self.stdPrecision);
		p = self.p;

		if B0 == None:
			B0 = self.B0[self.stdPrecision];
		#A priori, T can be chosen arbitrarily such that T>=2.
		#However it may need to be adjusted.
		T = 10;
		while True:
		
			if verbose:
				print "B0 =",B0;
				print "T =",T;
				print "stdPrecision =",self.stdPrecision;

			neededPrec = (RIFstd(T*B0).log2().upper()/self.stdPrecision).ceil()*self.stdPrecision;
			safetyBits = self.stdPrecision; 

			if verbose:
				print "neededPrec =",neededPrec,", safetyBits =",safetyBits;

			while True:
				precisionsSoFarThatAreUseful = [pr for pr in self.Alpha.keys() if pr>=neededPrec+2*safetyBits];
				if precisionsSoFarThatAreUseful != []:
					prec = min(precisionsSoFarThatAreUseful);
					if verbose:
						print "Need precision",neededPrec+4*safetyBits,"and we can take already used precision",prec;
				else:
					prec = neededPrec + 3*safetyBits; #Do a bit more in case we later need to increase T.
					if verbose:
						print "Need to redo some numerical computations with higher precision",prec;
					self.init_for_higher_precision(prec);

				RIFprec = RealIntervalField(prec);
				B0 = RIFprec(B0.upper());
				if verbose:
					print "More precise B0 =",B0;

				cIndex = self.Sigma.index(sigma_c);
				delta1 = self.delta_c_j(cIndex,j1,prec=prec);
				delta2 = self.delta_c_j(cIndex,j2,prec=prec);
				delta = delta2/delta1;

				if verbose:
					print "delta =",delta;
					print "Absolute diameter of delta:",delta.absolute_diameter();

				#Check precision for delta:
				if not delta.absolute_diameter() < RIFprec(1/(T*B0)):
					#Restart with larger precision:
					if verbose:
						print "Precision of delta was not enough! Redo with higher precision...";
					safetyBits += self.stdPrecision;
					neededPrec += ceil(log(delta.absolute_diameter() * RIFprec(T) * B0,2).upper());
					continue;
				
				theta1 = self.theta_c_j(cIndex,j1,prec=prec);
				theta2 = self.theta_c_j(cIndex,j2,prec=prec);
				lamda = (delta2*theta1-delta1*theta2)/delta1;
				if verbose:
					print "lamda =",lamda;
				#print "Absolute diameter of log_abs_gamma_c_k[1,0] =",self.log_abs_gamma_c_k(0,0,prec).absolute_diameter();
				#print "Absolute diameter of Alpha[0,0] =",self.Alpha[prec][0,0].absolute_diameter();
				#print "Absolute diameter of delta1 =",delta1.absolute_diameter();
				#print "Absolute diameter of delta2 =",delta2.absolute_diameter();
				#print "Absolute diameter of theta1 =",theta1.absolute_diameter();
				#print "Absolute diameter of theta2 =",theta2.absolute_diameter();
				#print "Absolute diameter of lamda =",lamda.absolute_diameter();

				#Find good rational approximation to delta with denominator r bounded by T*B0:
				nTerms = 1;

				#In Sage version 6.5:
				cf = continued_fraction(delta.center());
				if debug:
					print "cf(delta):",cf;
					print "Bound for r:",T*B0;

				#r = NaN;
				r = cf.denominator(0); #which should be 1
				
				while True:
					#The following still worked in Sage 6.4.1, but not in Sage 6.5:
					#cf = continued_fraction(delta.center(),nterms=nTerms);   #bits=neededPrec+3*safetyBits);
					
					#if len(cf) == 0:
					#	print "Continued fractions where trivial. Redo with higher precision..."
					#	safetyBits *= 2;
					#	continue;

					#In version 6.4.1:
					#r = cf.denominator();
					#In version 6.5:
					rNew = cf.denominator(nTerms);
					if debug:
						print "rNew =",rNew;
					if rNew == r and nTerms > 1: #note that it may happen that nTerms = 1 and r=rNew=1, in which case we want to continue with larger nTerms!)
						break;
					if not rNew <= T*B0:
						break;
					r = rNew;
					nTerms += 1;
					continue;				

				try:
					if debug:
						print "r =",r;
					rTimesDelta_round = (r*delta).unique_round();
					if verbose:
						print "Rounding of r*delta worked";
					rTimesLamda_round = (r*lamda).unique_round();
					if verbose:
						print "Rounding of r*lambda worked";
					roundingFailed = False;
				except ValueError:
					roundingFailed = True;
				if roundingFailed:
					if verbose:
						print "Rounding failed. Redo with higher precision..."
					safetyBits *= 2;
					continue;
				if debug:
					print "Check whether r*delta is close enough to an integer:";
					print "(rTimesDelta_round - r*delta).abs():",(rTimesDelta_round - r*delta).abs();
					print "1/(RIFprec(T)*B0):",1/(RIFprec(T)*B0);
				if not (rTimesDelta_round - r*delta).abs() < 1/(RIFprec(T)*B0):
					if verbose:
						print "r*delta not close enough to an integer. Redo with higher precision...";
					safetyBits *= 2;
					continue;

				if verbose:
					print "Precision is enough!";
				break; #Precision is enough!
				

			rTimesLamda_distToZZ = (rTimesLamda_round-r*lamda).abs();
			if not (rTimesLamda_round-r*lamda).abs() >= 2/T:
				if verbose:
					print "Unlucky: r*lamda too close to an integer. Redo with larger T...";
				T *= 10;
				continue;

			if verbose:
				print "T is also large enough!";
			break; #T is large enough!

		#print "debug123",type(delta);
		#print "debug123",type(self.Theta[prec]);
		#print "debug123",type(T);
		#print "debug123",type(B0);
		if verbose:
			print rTimesLamda_distToZZ;
		q_BakerBound = RIFstd((RIFstd(3.2)*(1+delta.abs())*self.Theta[prec]*T*B0)/(rTimesLamda_distToZZ-1/T))^(-p);
		#print "debug 234";

		q_BakerBound = min(q_BakerBound,RIFstd(self.Theta[prec])^(-p),RIFstd(2^(-p)));
		#print "debug 23424";

		#Obtain improved B0:
		Xi = log(1/q_BakerBound);
		#print "debug 234wer24";
		B1 = delta1.abs()*Xi + theta1.abs() + 3.2;
		#print "debug 234sadf24";
		B1 = RIFstd(min(B0,B1).upper());
		#print "debug asdf";
			
		if verbose:
			print "Reduced bounds:";
			print "B1 =",B1;
			print "q_BakerBound =",q_BakerBound;

		return B1,q_BakerBound;

	def init_for_higher_precision(self,prec):
		'''
		Needs to be called when you want to use numerical variables
		of X_ns^+(p) with a new precision prec.
		'''
		self.init_Eta_and_Alpha(prec=prec);
		self.init_variables_around_bakerBound(prec=prec);

	def showMemoryUsage(self):

		print "Memory usage of X:";
		for var_name,variable in vars(self).iteritems():
			if hasattr(variable,"__sizeof__"):
				print var_name,"has size",variable.__sizeof__();

	def reduction_from_B0old_to_B0new(self):
		#The first and second arxiv versions of the paper use slightly
		#different initial bounds B0.
		#The computations up to p=97 were done using the old bound.
		#So it remains to reduce from B0new to B0 (where B0=B0old).
		
		X = self;
		p = X.p;
		Sigma = X.Sigma;
		#Take one sigma per cusp:
		Sigma0 = [Sigma[i*p] for i in range(ZZ((p-1)/2))];
		
		prec = self.stdPrecision;
		B0old = X.B0_old1[prec]; #from first arxiv version
		B0new = X.B0[prec]; #from second arxiv version
		if not B0new < B0old:
			print "Need to reduce B0new at each cusp (with any sigma) until it's below B0:"
			for sigma in Sigma0:
				print "sigma:",sigma.list();
				B1new = B0new;
				while not B1new < B0old:
					B1new,q1new = X.reduction_of_BakerBound_B0_and_qc_at_sigmaC(sigma,B0=B1new);
					#raise Exception("TODO...");
		print "done";

########################################################################
### Some calculus that was done in order to write the program: #########
########################################################################

def calculationForGammaForOptimalSystemOfRepresentatives():
	'''
	This is not needed in the program.
	(Is was used in writing the program.)
	'''
	var("k1 k2 a b c d s t at ct p ");
	kx = s*at + c*k1;
	ky = t*at - a*k1;
	kz = s*ct + c*k2;
	kw = t*ct - a*k2;
	null = p*(kx*kw - ky*kz)+kx+kw;
	print null.expand();

	#One has to choose k1 and k2 such that null becomes zero.
	#This makes the determinant of gamma equal to one.
	#This turns out to be equivalent to f1*k1 + f2*k2 = k
	#with f1, f2 and k as defined in the program above.
	return null

#calculationForGammaForOptimalSystemOfRepresentatives();
	
def calculationForLogEstimates():
	var("x y t phi alpha beta r ");
	#f(t) = abs(log(1+phi*t));
	phi = exp(I*beta);
	#1+phi*t = r*exp(I*alpha) = x + I*y
	x = 1 + t*cos(beta);
	y = t*sin(beta);
	r = sqrt(x^2+y^2);
	alpha = arctan(y/x);
	#log(1+phi*t) = log(r) + I*alpha, thus
	#f := |log(1+phi*t)| = sqrt((log(r))^2 + alpha^2):
	f = sqrt((log(r))^2 + alpha^2);
	#print f;
	#show(latex((f^2).derivative(t).expand()));

	#plot3d(lambda x,y:abs(log(x+I*y)),(0,2),(-1,1))
	P = plot3d(lambda t,phi:abs(log(1+t*exp(I*phi))),(0,1),(0,2*pi))
	P.show();
	return;

#calculationForLogEstimates();

def test_dualEllipsoids1():
	'''
	The following finds an ellipsoid that contains two ellipsoids that are centered at zero:
	'''
	Q1 = matrix(QQ,2,2,[2,1,1,3]);
	Q2 = matrix(QQ,2,2,[5,2,2,1]);
	Qall = (Q1.inverse()+Q2.inverse()).inverse();
	print Q1.det(),Q2.det(),Qall.det();
	P = contour_plot(vector([x,y])*Q1*vector([x,y]),(x,-5,5),(y,-5,5),contours=[0,1],fill=False);
	P += contour_plot(vector([x,y])*Q2*vector([x,y]),(x,-5,5),(y,-5,5),contours=[0,1],fill=False);
	P += contour_plot(vector([x,y])*Qall*vector([x,y]),(x,-5,5),(y,-5,5),contours=[0,1],fill=False);
	P.show();

def test_dualEllipsoids2():
	'''
	The following is complete non-sense:
	'''
	Q = matrix(QQ,2,2,[2,1,1,3]);
	a = vector([3,1]);
	print Q*a, a*Q, a*Q*a;
	Qtilda1 = block_matrix(QQ,2,2,[Q,(Q*a).column(),(a*Q).row(),matrix(QQ,1,1,[a*Q*a])]);
	Qtilda2 = block_matrix(QQ,2,2,[Q,(-Q*a).column(),(-a*Q).row(),matrix(QQ,1,1,[a*Q*a])]);
	QtildaAll = (Qtilda1.adjoint()+Qtilda2.adjoint()).adjoint();
	Qall = QtildaAll.submatrix(0,0,2,2);
	print Qall;
	#print Q1.det(),Q2.det(),Qall.det();
	P = contour_plot((vector([x,y])-a)*Q*(vector([x,y])-a),(x,-5,5),(y,-5,5),contours=[0,1],fill=False);
	P += contour_plot((vector([x,y])+a)*Q*(vector([x,y])+a),(x,-5,5),(y,-5,5),contours=[0,1],fill=False);
	P += contour_plot(vector([x,y])*Qall*vector([x,y]),(x,-5,5),(y,-5,5),contours=[0,1],fill=False);
	P.show();

def test_dualEllipsoids3():
	'''
	The following takes considers translates E+a and E-a of some
	 ellipsoid centered at the origin, and computes an ellipsoid
	 that contains both translates.
	'''
	
	Q = matrix(QQ,2,2,[2,1,1,3]);
	a = vector([3,1]);
	#Q = matrix(QQ,2,2,[5,2,2,3]);
	#a = vector([0,1]);
	#Q = matrix(QQ,2,2,[5,2,2,3]);
	#a = vector([0,3]);
	#proj = block_matrix(QQ,1,2,[identity_matrix(2),a.column()]);
	#Qtilda = block_diagonal_matrix(Q,identity_matrix(1))/2;
	#QtildaDual = Qtilda.inverse();
	#QallDual = proj * QtildaDual * proj.transpose();
	#Qall = QallDual.inverse();
	#print Qall;
	Qall = smallEllipsoidContaining_AplusQ_and_AminusQ(Q,a);
	#print Q1.det(),Q2.det(),Qall.det();
	P = contour_plot((vector([x,y])-a)*Q*(vector([x,y])-a),(x,-5,5),(y,-5,5),contours=[0,1],fill=False);
	P += contour_plot((vector([x,y])+a)*Q*(vector([x,y])+a),(x,-5,5),(y,-5,5),contours=[0,1],fill=False);
	P += contour_plot(vector([x,y])*Qall*vector([x,y]),(x,-5,5),(y,-5,5),contours=[0,1],fill=False);
	P.show();

#test_dualEllipsoids1();
#test_dualEllipsoids2();
#test_dualEllipsoids3();

########################################################################
### Extra search over some small set of j's: #.##########################
########################################################################

def test_j_and_p(j,p):
	'''
	Checks the image type of the mod-p Galois representation
	of ONE elliptic curve with j-invariant j.
	'''
	E = EllipticCurve_from_j(QQ(j));
	#print "j =",j,":",E, "-> is minimal =",E.is_minimal();
	r = E.galois_representation();
	imageType = r.image_type(p);
	output = "";
	if not imageType.startswith("The image is all of GL_2(F_"):
		output += "j = "+str(j);
		if j != -0:
			output += " = " + str(factor(j));
		output += ", hasCM = "+str(E.has_cm());
		output += ", "+imageType;
		output += "\n";
	return output;
	#sys.stdout.flush();
			
def testEllipticCurvesWithSmallJ(p):
	print "Elliptic curves with j in [0,1728] with image not all of GL_2(F_p) when p =",p,":";
	for j in range(0,1729):
		E = EllipticCurve_from_j(QQ(j));
		#print "j =",j,":",E, "-> is minimal =",E.is_minimal();
		r = E.galois_representation();
		imageType = r.image_type(p);
		if not imageType.startswith("The image is all of GL_2(F_"):
			print "j =",j,
			print ", hasCM =",E.has_cm(),
			print ", ",imageType;
		sys.stdout.flush();

#for p in prime_range(1,100):
#	testEllipticCurvesWithSmallJ(p);

########################################################################
### Fincke-Pohst: ######################################################
########################################################################

def QuadraticFormToFPForm_rif(A,base_field=None):
    '''
    We compute a symmetric matrix Q such that for any vector x,
     x^t*A*x = \sum_{i=1}^n q_{ii} * (x_i + \sum_{j=i+1}^n q_{ij}*x_j)^2
    Everything is done over the base_field.
     (If None, then take the base ring of A, which should be a field).
    '''
    
    n = A.nrows();
    if base_field==None:
		base_field = A.base_ring();
    Q = matrix(base_field,n,n);
    
    #Following Fincke-Pohst's paper:
    #Q = copy(A);
    #Q = Q.change_ring(QQ);
    #for i in range(n):
    #    for j in range(i+1,n):
    #        Q[j,i] = Q[i,j]; #used just as a temporary buffer
    #        Q[i,j] = Q[i,j] / Q[i,i];
    #for i in range(n):
    #    for k in range(i+1,n):
    #        for l in range(k,n):
    #            Q[k,l] = Q[k,l] - Q[k,i]*Q[i,l];
    #for i in range(n):
    #    for j in range(i):
    #        Q[i,j] = 0; #the lower diagonal is redundant now
            
    for i in range(n):
        for k in range(i+1):
            s = 0;
            for j in range(k):
                s = s + Q[j,i]*Q[j,k]*Q[j,j];
            if i==k:
                Q[i,i] = A[i,i]-s;
            else:
                Q[k,i] = (A[k,i]-s)/Q[k,k];
            
    return Q;            
   
def My_FinckePohst_ViaGramMatrix_rif(A,boundForNormSquared,center=0,solutionList=None,maxNumSolutions=None,finalTransformation=None,callback=None,callbackArgs=None,lllReduce=True,breakSymmetry=True):
	'''
	Input: A - a positive definite matrix with RIF coefficients
	       boundForNormSquared - a non-negative RIF
	       maxNumSolutions - a non-negative integer, or None if this should not be bounded.
	       solutionList - either None or [], depending on whether it shall be filled with solutions or not
	       center - a RIF vector, or 0 (which means the zero vector).
	       breakSymmetry - a boolean: If True and if center.base_ring()==ZZ, then the output will consider only one solution among each pair of antipodal solutions.
	       
	Output: If callback function appears as a parameter:
	       [notTooManySolutions,numSolutions], where
	       notTooManySolutions - true iff number of solutions is at most maxNumSolutions or maxNumSolutions==None
	       numSolutions - number of solutions that have been traversed
	Furthermore, for each solution:
	    if solutionList != None, then the solution will be appended to solutionList, and
	    if callback != None, then callback(solution,callbackArgs) will be called.

	In a previous version of this function:
	Output: If callback function does not appear as a parameter:
	       [bufferWasBigEnough,solutions], where
	       bufferWasBigEnough - true iff number of solutions is at most maxNumSolutions
	       solutions - a list whose entries are the integral vectors x!=center with
	              (x-center)^t * A * (x-center) <= boundForNormSquared
	'''

	def traverseShortVectors(Q,n,k,x,center,shortVectors,RemainingBoundForNormSquared,numSolutionsList,maxNumSolutions,xIsCenterSoFar,finalTransformation,callback,callbackArgs):
		'''
		Fills shortVectors or calls callback()
		Returns [bufferBigEnoughSoFar,newNumSolutionsSoFar]
		'''
		
		if k==-1:
			if xIsCenterSoFar:
				return True; #We don't consider the center point as a solution.
			if maxNumSolutions!=None and numSolutionsList[0]>=maxNumSolutions:
				return False; #As there are too many solutions.
			numSolutionsList[0] = numSolutionsList[0]+1;
			y = vector(x);
			if finalTransformation != None:
				y = finalTransformation * y;
			if shortVectors != None:
				shortVectors.append(y);
			if callback != None:
				callback(y,callbackArgs);
			return True;
			
		u = -center[k];
		for j in range(k+1,n):
			u = u + Q[k,j] * (x[j]-center[j]);
		if not contains_only_finite_rifs([u]):
			numSolutionsList[0] = NaN;
			return False;
		uCeil = u.upper().ceil();
		uFloor = u.lower().floor();

		xk0_up = -uFloor;
		xk0_down = -uCeil;

		#We check all x[k] between xk0_down and xk0_up no matter what.
		#On both sides outside that interval, t>RemainingBoundForNormSquared is a breaking condition,
		# that is, once this condition does not hold anymore,
		# we are free to stop searching for further x[k] in this direction.
		
		x[k] = xk0_down;
		#print "debug14",k,u,xk0_up,xk0_down;
		while True:
			t = Q[k,k] * (x[k] + u)^2;
			#print "t =",t;
			if t>RemainingBoundForNormSquared:
				if x[k] >= xk0_up: 
					break;
				x[k] += 1;
				continue;	
			#print k,x[k],t;
			if not traverseShortVectors(Q,n,k-1,x,center,shortVectors,RemainingBoundForNormSquared-t,numSolutionsList,maxNumSolutions,xIsCenterSoFar and x[k]==center[k], finalTransformation,callback,callbackArgs):
				return False; #too many solution found
			x[k] = x[k] + 1;            

		if not xIsCenterSoFar: #break antipodal symmetry: if x was so far center, then u is zero here, and we iterate only over x[k]>=0.
			x[k] = xk0_down-1;
			while True:
				t = Q[k,k] * (x[k] + u)^2;
				if t>RemainingBoundForNormSquared:
					break;
				if not traverseShortVectors(Q,n,k-1,x,center,shortVectors,RemainingBoundForNormSquared-t,numSolutionsList,maxNumSolutions,False, finalTransformation,callback,callbackArgs):
					return False; #too many solutions found
				x[k] = x[k] - 1;            
			
		return True; #prescribed bufferSize (maxNumSolutions) was enough
		
	global aLotOutput;
	#if aLotOutput:
	#    print "A =",A;    

	n = A.nrows();

	if center == 0:
		center = zero_vector(n);

	if center.base_ring() != ZZ: #We check this here, as after multiplication with U, center.base_ring() might become QQ! (even though U's base_ring is ZZ as well...)
		breakSymmetry = False;

	if lllReduce:
		U = LLL_gram_rif(A);
		A = U.transpose() * A * U;
		if finalTransformation != None:
			finalTransformation = finalTransformation * U;
		else:
			finalTransformation = U;
		if not center.is_zero():
			center = U.inverse() * center;

	Q = QuadraticFormToFPForm_rif(A);
	#print "center =",center;
	#print "Q =\n", Q;

	x = range(n);
	numSolutionsList = [0];
	bufferWasBigEnough = traverseShortVectors(Q,n,n-1,x,center,solutionList,boundForNormSquared,numSolutionsList,maxNumSolutions,breakSymmetry,finalTransformation,callback,callbackArgs);
		
	#return [bufferWasBigEnough,matrix(shortVectors).transpose()];
	return [bufferWasBigEnough,numSolutionsList[0]];

def normSquared(v,Q=None):
	'''
	INPUT
	* v - a vector
	* Q - a quadratic form, Q=None stands for the standard
	      Euclidean quadratic form (i.e. identity matrix).
	OUTPUT
	<v,v>_Q = v*Q*v
	'''
	if Q == None:
		return sum([x^2 for x in v]);
	else:
		v = vector(v);
		return v * Q * v;

########################################################################
### Ellipsoid basics: ##################################################
########################################################################

def smallestEllipsoidContaining_Q1_times_Q2(Q1,Q2):
	'''
	Given Ellipsoids E_Q1 = {x: x^t*Q1*x <= 1} and E_Q2,
	 compute the quadratic form Qprod giving rise to the
	 smallest volume ellipoid E_Qprod that contains
	 E_Q1 times E_Q2.
	'''
	
	n1 = Q1.nrows();
	n2 = Q2.nrows();
	Qprod = block_diagonal_matrix(n1*Q1,n2*Q2)/(n1+n2);
	return Qprod;

def smallEllipsoidContaining_AplusQ_and_AminusQ(Q,a):
	'''
	Given an Ellipsoid E_Q = {x: x^t*Q*x <= 1}, and a vector a,
	 consider the translates E+a and E-a.
	The output is a matrix Qboth such that the corresponding
	 ellipsoid E_Qboth contains both translates E+a and E-a.
	'''
	
	#print "Start computing ellipsoid containing E+a and E-a...",
	
	#We construct E_Q' as the image of an ellipsoid E_Qtilda of one dimension higher:
	#E_Qtilda is the smallest volume ellipsoid that contains E_Q \times [-1,+1].
	Qtilda = smallestEllipsoidContaining_Q1_times_Q2(Q,identity_matrix(1));
	#E_Qboth is now the projection of E_Qtilda via the map "projection":
	projection = block_matrix(a.base_ring(),1,2,[identity_matrix(len(a)),a.column()]);
	#I only see how to compute the projection in the following way:
	#I dualize E_Qtilda, then pull it back via an inclusion (the dual projection),
	# and then I dualize again:
	QtildaDual = Qtilda.inverse();
	QbothDual = projection * QtildaDual * projection.transpose();
	Qboth = QbothDual.inverse();

	#print "done.";
	
	return Qboth;

def volumeOfBall(dim,radius=1):
	'''
	Computes the volume of a 'dim'-dimensional ball of radius 'radius'
	in Euclidean space.
	'''
	return RIF(pi^(dim/2)*radius^dim/gamma(dim/2+1));

def volumeOfEllipsoid(Q,radiusSquared):
	'''
	Computes the volume of ellipsoid {x: x^t * Q * x <= radiusSquared}
	in Euclidean space.
	'''
	return RIF(volumeOfBall(Q.nrows(),sqrt(radiusSquared)) / sqrt(abs(det(Q))));

########################################################################
### Applying Fincke-Pohst to our setting: ##############################
########################################################################

def logU_asCoordinatewiseRif(q_rif,sigmaC,X,prec=None,verbose=False,coordinates="all",n_summands=1):
	'''
	We consider the curve $q\mapsto log of abs of (U^\sigma(q))_{\sigma\in Gal(K/Q)}$,
	where q ranges over the given interval q_rif.
	The output is an axis-parallel cube in R^d (i.e. a vector of rifs) that contains this curve segment.
	Note: if coordinates is not "all", the output will still be only a list, with the other coordinates just missing. Beware that this means a relabeling of the coordinates!
	
	Input:
	 q_rif: A real interval, which is the domain of the given curve segment.
	 sigmaC: sigmaC*D is the currently considered domain of X_ns^+(p).
	 rightAction_GH_on_Mp
	 O = Orbit0: A fixed orbit 
	'''
	
	#For P in \sigma_1(P) and \sigma_2 in G_H:
	#
	#u_O ^ {\sigma_2} (P) =
	#u_(O\sigma_2) (P) =
	#\prod_{a\in O\sigma_2} g_{\tilde a}^m (q) =
	#\prod_{a\in O\sigma_2\sigma_1} g_{\tilde a}^m (q\circ \sigma_1^{-1}).

	debug = False;
	#debug = True;
	#verbose = True;

	if prec == None:
		prec = X.stdPrecision;
	RIFprec = RealIntervalField(prec);
	CIFprec = ComplexIntervalField(prec);

	piPrec = CIFprec(pi);

	p = X.p;
	GalK_lift = X.GalK_lift;
	m = X.m; 

	qMin = RIFprec(q_rif.lower());
	qMax = RIFprec(q_rif.upper());
	tau1, tau2 = convert_realQparameter_to_tau(q_rif);
	tau = tau1.union(tau2);

	pi2Itau = 2*piPrec*I*tau; 
	if verbose:
		print "qMin =",qMin;
		print "qMax =",qMax;
		print "tau =",tau;
		print "exp(pi2Itau) =",exp(pi2Itau);
	abs_q = abs(q_rif);
	abs_q_up = RIFprec(abs_q.upper());
	if verbose:
		print "abs_q =",abs_q;
		print "abs_q_up =",abs_q_up;

	#sigmaLsigmaC = sigmaL * sigmaC;

	#Number of terms in the expansion of u_o (that's the n in (23):
	# (Make this larger in case you think the error term makes the ellipsoid too large.)
	#n_summands = 1;
	
	d = X.d; # = len(GalK_lift)
	cIndex = X.Sigma.index(sigmaC);
	
	if coordinates == "all":
		coordinates = range(d);

	log_Uk = [RIFprec(0) for k in coordinates];

	for index in range(len(coordinates)):
		kIndex = coordinates[index];
		#print "kIndex =",kIndex;
		k = GalK_lift[kIndex];
		sigma_k = X.lift_GalKIndex_to_G[kIndex]

		#The following is not needed anymore, as it is already in the memory:
		ord_c_Uk = RIFprec(0);
		log_abs_gamma_c = RIFprec(0);

		#a_list = []; #debugging...
		#aInOrbit = None;
		
		for a0_in_Mp in X.Orbit0:
			#a_in_Mp = X.rightAction_GH_on_Mp.action(sigmaC*sigma_k,a0_in_Mp); #Wrong order of sigmas!!!
			a_in_Mp = X.rightAction_GH_on_Mp.action(sigma_k*sigmaC,a0_in_Mp);
			a1,a2 = lift_of_Mp_to_ZZ2overP(a_in_Mp,p);

			#a_list.append(a_in_Mp);
			#aInOrbit = a_in_Mp;
			#print "debug3 OR:",X.rightAction_GH_on_Mp.OrbitRepresentative(X.rightAction_GH_on_Mp.OrbitOfM(aInOrbit));
		

			#print "a_in_Mp 2 =",a_in_Mp;
			
			#The following is not needed anymore, as it is already in the memory:
			#Above Proposition 6.6:
			l_a = bernoulli_polynomial(a1,2)/2;
			ord_c_Uk += p*m*l_a;
			#According to section "Approximate Formulas" of the paper:
			if a1 != 0:
				log_abs_gamma_a = log(abs(exp(piPrec*I*a2*(a1-1))));
			else:
				log_abs_gamma_a = log(abs(exp(piPrec*I*a2*(a1-1))*(1-X.exp_2pi_I_a2[prec][a2])));
			log_abs_gamma_c += RIFprec(m * log_abs_gamma_a); #taking RIF here seems to be necessary for some reason I don't understand (otherwise it's treated sometimes as a symbolic expression).


			#According to section "Approximate Formulas" of the paper:
			sumLogs_Main_tau = RIFprec(0);
			for i in range(n_summands):
				if i+a1 > 0: #otherwise this term is taken care of in log_abs_gamma_a:
					sumLogs_Main_tau += log(abs(CIFprec(1-exp(pi2Itau*(i+a1))*X.exp_2pi_I_a2[prec][a2])));
				#else:
				#	#print "debugT:",log(abs(CIFprec(1)-qMin^(i+a1)*X.exp_2pi_I_a2[prec][a2]));
				#	#print "debugq:",qMin^(i+a1);
				sumLogs_Main_tau += log(abs(CIFprec(1-exp(pi2Itau*(i+1-a1))*X.exp_2pi_I_a2[prec][-a2])));

			sumLogs_Remainder = RIFprec((abs_q_up^(n_summands+a1)+abs_q_up^(n_summands+1-a1))/(1-abs_q_up)^2);

			sumLogs = sumLogs_Main_tau + sumLogs_Remainder.union(-sumLogs_Remainder);

			if verbose:
				print "log_abs_gamma_a =",log_abs_gamma_a;
				print "log_abs_gamma_c =",log_abs_gamma_c;
				print "sumLogs_Main_tau =",sumLogs_Main_tau;
				print "sumLogs_Remainder =",sumLogs_Remainder;
				print "sumLogs =",sumLogs;

			log_Uk[index] += m * sumLogs;
			if verbose:
				print "m * sumLogs =",m * sumLogs;

		log_Uk[index] += log_abs_gamma_c;
		log_Uk[index] += ord_c_Uk/p * log(abs_q);
		if verbose:
			print "ord_c_Uk/p * log(abs_q) =",ord_c_Uk/p * log(abs_q);

	if verbose:
		print "log_Uk[coordinates]:",log_Uk;
		
	return log_Uk;

def DlogU_asCoordinatewiseRif(q_rif,sigmaC,X,prec=None,verbose=False,coordinates="all"):
	'''
	We consider the curve $q\mapsto log of abs of (U^\sigma(q))_{\sigma\in Gal(K/Q)}$,
	where q ranges over the given interval q_rif.
	The output is an axis-parallel cube in R^d (i.e. a vector of rifs) that contains all DERIVATIVES along the curve segment.
	Note: if coordinates is not "all", the output will still be only a list, with the other coordinates just missing. Beware that this means a relabeling of the coordinates!

	Input:
	 q_rif: A real interval, which is the domain of the given curve segment.
	 sigmaC: sigmaC*D is the currently considered domain of X_ns^+(p).
	 rightAction_GH_on_Mp
	 O = Orbit0: A fixed orbit 
	'''
	
	#For P in \sigma_1(P) and \sigma_2 in G_H:
	#
	#u_O ^ {\sigma_2} (P) =
	#u_(O\sigma_2) (P) =
	#\prod_{a\in O\sigma_2} g_{\tilde a}^m (q) =
	#\prod_{a\in O\sigma_2\sigma_1} g_{\tilde a}^m (q\circ \sigma_1^{-1}).

	debug = False;

	if prec == None:
		prec = X.stdPrecision;
	RIFprec = RealIntervalField(prec);
	CIFprec = ComplexIntervalField(prec);

	piPrec = CIFprec(pi);

	p = X.p;
	GalK_lift = X.GalK_lift;
	m = X.m; 

	qMin = RIFprec(q_rif.lower());
	qMax = RIFprec(q_rif.upper());
	tau1, tau2 = convert_realQparameter_to_tau(q_rif);
	tau = tau1.union(tau2);

	pi2Itau = 2*piPrec*I*tau; 
	if verbose:
		print "qMin =",qMin;
		print "qMax =",qMax;
		print "tau =",tau;
		print "exp(pi2Itau) =",exp(pi2Itau);
	abs_q = abs(q_rif);
	abs_q_up = RIFprec(abs_q.upper());

	#sigmaLsigmaC = sigmaL * sigmaC;

	#Number of terms in the expansion of u_o (that's the n in (23):
	# (Make this larger in case you think the error term makes the ellipsoid too large.)
	n_summands = 1;
	
	d = X.d; # = len(GalK_lift)
	cIndex = X.Sigma.index(sigmaC);
	
	if coordinates == "all":
		coordinates = range(d);

	Dlog_Uk = [RIFprec(0) for k in coordinates];

	for index in range(len(coordinates)):
		kIndex = coordinates[index];
		#print "kIndex =",kIndex;
		k = GalK_lift[kIndex];
		sigma_k = X.lift_GalKIndex_to_G[kIndex]

		#The following is not needed anymore, as it is already in the memory:
		ord_c_Uk = 0;
		log_abs_gamma_c = RIFprec(0);

		#a_list = []; #debugging...
		#aInOrbit = None;
		
		for a0_in_Mp in X.Orbit0:
			#a_in_Mp = X.rightAction_GH_on_Mp.action(sigmaC*sigma_k,a0_in_Mp); #Wrong order of sigmas!!!
			a_in_Mp = X.rightAction_GH_on_Mp.action(sigma_k*sigmaC,a0_in_Mp);
			a1,a2 = lift_of_Mp_to_ZZ2overP(a_in_Mp,p);

			#The following is not needed anymore, as it is already in the memory:
			#Above Proposition 6.6:
			l_a = bernoulli_polynomial(a1,2)/2;
			ord_c_Uk += p*m*l_a;

			sumDLogs_Main_tau = RIFprec(0);
			for i in range(n_summands):
				if i+a1 > 0: #otherwise this term is taken care of in log_abs_gamma_a:
					sumDLogs_Main_tau += CIFprec(-exp(pi2Itau*(i-1+a1))*(i+a1)*X.exp_2pi_I_a2[prec][a2]) / CIFprec(1-exp(pi2Itau*(i+a1))*X.exp_2pi_I_a2[prec][a2]);
				sumDLogs_Main_tau += CIFprec(-exp(pi2Itau*(i-a1))*(i+1-a1)*X.exp_2pi_I_a2[prec][-a2]) / CIFprec(1-exp(pi2Itau*(i+1-a1))*X.exp_2pi_I_a2[prec][-a2]);

			sumDLogs_Remainder = RIFprec(1/(1-abs_q)^2);
			for i in range(1,n_summands):
				sumDLogs_Remainder -= RIFprec(i * abs_q^(i-1));
			sumDLogs_Remainder *= RIFprec(1/(1-abs_q^(n_summands+a1)) + 1/(1-abs_q^(n_summands+1-a1)));

			sumDLogs = sumDLogs_Main_tau + sumDLogs_Remainder.union(-sumDLogs_Remainder);

			if verbose:
				print "sumLogs_Remainder =",sumDLogs_Remainder;

			Dlog_Uk[index] += m * sumDLogs;

		Dlog_Uk[index] += ord_c_Uk/p * 1/q_rif;
		
	return [x.real() for x in Dlog_Uk];
	
def ellipsoid_around_logU(q_rif,sigmaC,X,prec=None,verbose=False,putMidpointErrorIntoEllipsoid=True):
	'''
	We consider the curve $q\mapsto log of abs of (U^\sigma(q))_{\sigma\in Gal(K/Q)}$,
	where q ranges over the given interval q_rif.
	The output is an ellipsoid in R^d that contains this curve segment.
	
	Input:
	 q_rif: A real interval, which is the domain of the given curve segment.
	 sigmaC: sigmaC*D is the currently considered domain of X_ns^+(p).
	 rightAction_GH_on_Mp
	 O = Orbit0: A fixed orbit 
	'''
	
	#For P in \sigma_1(P) and \sigma_2 in G_H:
	#
	#u_O ^ {\sigma_2} (P) =
	#u_(O\sigma_2) (P) =
	#\prod_{a\in O\sigma_2} g_{\tilde a}^m (q) =
	#\prod_{a\in O\sigma_2\sigma_1} g_{\tilde a}^m (q\circ \sigma_1^{-1}).

	debug = False;

	if prec == None:
		prec = X.stdPrecision;
	RIFprec = RealIntervalField(prec);
	CIFprec = ComplexIntervalField(prec);

	piPrec = CIFprec(pi);

	p = X.p;
	GalK_lift = X.GalK_lift;
	m = X.m; 

	qMin = RIFprec(q_rif.lower());
	qMax = RIFprec(q_rif.upper());
	tau1, tau2 = convert_realQparameter_to_tau(q_rif);
	pi2Itau1 = 2*piPrec*I*tau1; 
	pi2Itau2 = 2*piPrec*I*tau2; 
	tau = tau1.union(tau2);
	pi2Itau = 2*piPrec*I*tau; 
	if verbose:
		print "qMin =",qMin;
		print "qMax =",qMax;
		print "tau1 =",tau1;
		print "tau2 =",tau2;
		print "exp(pi2Itau1) =",exp(pi2Itau1);
		print "exp(pi2Itau2) =",exp(pi2Itau2);
	abs_qMin = abs(qMin);
	abs_qMax = abs(qMax);
	abs_q = abs(q_rif);
	abs_q_up = RIFprec(abs_q.upper());

	#sigmaLsigmaC = sigmaL * sigmaC;

	#Number of terms in the expansion of u_o (that's the n in (23):
	# (Make this larger in case you think the error term makes the ellipsoid too large.)
	n_summands = 1;
	
	d = X.d; # = len(GalK_lift)
	cIndex = X.Sigma.index(sigmaC);

	log_Uk_qMin = [RIFprec(0) for k in range(d)];
	log_Uk_qMax = [RIFprec(0) for k in range(d)];

	log_Uk_error = [RIFprec(0) for k in range(d)];
	

	for kIndex in range(d):
		#print "kIndex =",kIndex;
		k = GalK_lift[kIndex];
		sigma_k = X.lift_GalKIndex_to_G[kIndex]

		#The following is not needed anymore, as it is already in the memory:
		ord_c_Uk = 0;
		log_abs_gamma_c = RIFprec(0);

		#a_list = []; #debugging...
		#aInOrbit = None;
		
		for a0_in_Mp in X.Orbit0:
			#a_in_Mp = X.rightAction_GH_on_Mp.action(sigmaC*sigma_k,a0_in_Mp); #Wrong order of sigmas!!!
			a_in_Mp = X.rightAction_GH_on_Mp.action(sigma_k*sigmaC,a0_in_Mp);
			a1,a2 = lift_of_Mp_to_ZZ2overP(a_in_Mp,p);

			#a_list.append(a_in_Mp);
			#aInOrbit = a_in_Mp;
			#print "debug3 OR:",X.rightAction_GH_on_Mp.OrbitRepresentative(X.rightAction_GH_on_Mp.OrbitOfM(aInOrbit));
		

			#print "a_in_Mp 2 =",a_in_Mp;
			
			#The following is not needed anymore, as it is already in the memory:
			#Above Proposition 6.6:
			l_a = bernoulli_polynomial(a1,2)/2;
			ord_c_Uk += p*m*l_a;
			#According to section "Approximate Formulas" of the paper:
			if a1 != 0:
				log_abs_gamma_a = log(abs(exp(piPrec*I*a2*(a1-1))));
			else:
				log_abs_gamma_a = log(abs(exp(piPrec*I*a2*(a1-1))*(1-X.exp_2pi_I_a2[prec][a2])));
			log_abs_gamma_c += RIFprec(m * log_abs_gamma_a); #taking RIF here seems to be necessary for some reason I don't understand (otherwise it's treated sometimes as a symbolic expression).



			#According to section "Approximate Formulas" of the paper:
			sumLogs_Main_tau = RIFprec(0);
			for i in range(n_summands):
				if i+a1 > 0: #otherwise this term is taken care of in log_abs_gamma_a:
					sumLogs_Main_tau += log(abs(CIFprec(1-exp(pi2Itau*(i+a1))*X.exp_2pi_I_a2[prec][a2])));
				sumLogs_Main_tau += log(abs(CIFprec(1-exp(pi2Itau*(i+1-a1))*X.exp_2pi_I_a2[prec][-a2])));

			##TODO: Improve this remainder term (this is (22) for m=0 applied twice!)
			sumLogs_Remainder = RIFprec((abs_q_up^(n_summands+a1)+abs_q_up^(n_summands+1-a1))/(1-abs_q_up)^2);

			sumLogs = sumLogs_Main_tau + sumLogs_Remainder.union(-sumLogs_Remainder);

			if verbose:
				print "sumLogs_Remainder =",sumLogs_Remainder;
				print "sumLogs =",sumLogs;				

			log_Uk_qMin[kIndex] += m * sumLogs;
			log_Uk_qMax[kIndex] += m * sumLogs;
			#log_Uk_error[kIndex] += m * (sumLogs_Main_error + sumLogs_Remainder);

		#a_list.sort();
		#print "debug a_list:",a_list;
		#print "debug OR:",X.rightAction_GH_on_Mp.OrbitRepresentative(X.rightAction_GH_on_Mp.OrbitOfM(aInOrbit));
		
		

		#log_abs_gamma_c = X.log_abs_gamma_c_k(cIndex,kIndex,prec=prec);
		#ord_c_Uk = X.ord_c_Uk(cIndex,kIndex);
		#print "Test log_abs_gamma_c: ",log_abs_gamma_c,"vs.",X.log_abs_gamma_c_k(cIndex,kIndex,prec=prec);
		#print "Test ord_c_Uk: ",ord_c_Uk,"vs.",X.ord_c_Uk(cIndex,kIndex);
		
		log_Uk_qMin[kIndex] += log_abs_gamma_c;
		log_Uk_qMax[kIndex] += log_abs_gamma_c;

		log_Uk_qMin[kIndex] += ord_c_Uk/p * log(abs_qMin);
		log_Uk_qMax[kIndex] += ord_c_Uk/p * log(abs_qMax);


		log_Uk_error[kIndex] = max(myRIF_radius(log_Uk_qMin[kIndex]),myRIF_radius(log_Uk_qMax[kIndex]));
		log_Uk_qMin[kIndex] = myRIF_center(log_Uk_qMin[kIndex]);
		log_Uk_qMax[kIndex] = myRIF_center(log_Uk_qMax[kIndex]);
		

	#We now computed a neighborhood of the curve segement.
	#It is the line segment between the vector log_Uk_qMin to the vector log_Uk_qMax, however
	# as these two vectors have only RIF-coordinates, it's actually a Minkowski sum of a proper line segment and a cube given by the errors.

	if verbose:
		print "Segment from:",log_Uk_qMin;
		print "          to:",log_Uk_qMax;
		print "Error box:",log_Uk_error;

	#Next we compute an ellipsoid that contains this polytope.
	#Of course we must take care for the numerical errors that occurred.

	#The midpoint of the final ellipsoid:
	midpoint_rif = [RIFprec((log_Uk_qMin[k]+log_Uk_qMax[k])/2) for k in range(d)];
	#print "midpoint_rif:",midpoint_rif;

	#Vector by which the error ellipsoid has to be translated:
	a1 = [log_Uk_qMax[k]-midpoint_rif[k] for k in range(d)];
	a2 = [midpoint_rif[k]-log_Uk_qMin[k] for k in range(d)];
	a_rif = [a1[k].union(a2[k]) for k in range(d)];
	a = vector([a_rif[k] for k in range(d)]);
	a_error = [myRIF_radius(a_rif[k]) for k in range(d)];

	if putMidpointErrorIntoEllipsoid:
		midpoint = vector([midpoint_rif[k].center() for k in range(d)]);
		midpoint_error = [myRIF_radius(midpoint_rif[k]) for k in range(d)];
		#print "midpoint_error:",midpoint_error;
		error = [RIFprec((log_Uk_error[k] + a_error[k] + midpoint_error[k]).upper()) for k in range(d)];
	else:
		midpoint = vector(midpoint_rif);
		error = [RIFprec((log_Uk_error[k] + a_error[k]).upper()) for k in range(d)];

	#The matrix whose corresponding ellipsoid bounds the error term:
	Qerror = diagonal_matrix([1/error[k]/d for k in range(d)]);
	Q_rif = smallEllipsoidContaining_AplusQ_and_AminusQ(Qerror,a);

	#The ellipsoid is now given by quadratic form Q and midpoint m.
	#The only problem is that the entries of Q are real intervals.

	if verbose:
		print "Ellipsoid Q_rif:\n",Q_rif
		print "Det(Q_rif) =",det(Q_rif);

	return (Q_rif, midpoint);

def ellipsoid_around_logU_inBcoordinates(q_rif,sigmaC,X,prec=None,verbose=False,putMidpointErrorIntoEllipsoid=True):
	'''
	Same as ellipsoid_around_logU() (see its documentation string),
	except that this ellipoid is transformed into b-coordinates via
	the matrix Alpha = Eta.inverse().
	'''

	if prec == None:
		prec = X.stdPrecision;
	RIFprec = RealIntervalField(prec);
	Q_rif, midpoint = ellipsoid_around_logU(q_rif,sigmaC,X,prec=prec,putMidpointErrorIntoEllipsoid=putMidpointErrorIntoEllipsoid);

	midpoint_B = X.Alpha[prec] * vector([RIFprec(x) for x in midpoint]);
	Q_B = X.Eta[prec].transpose() * Q_rif * X.Eta[prec];
	if verbose:
		print "Volume of ellipsoid Q_B:",volumeOfEllipsoid(Q_B,1);
	return Q_B,midpoint_B;
	
def findCandidatesB_inQinterval_at_sigmaC(q_rif,sigmaC,X,candidatesB,maxCandidates,prec=None,verbose=False):
	'''
	Find all candidate vectors B that possibly come from integral points
	 at sigma_C * D (D being the standard fundamental domain for X(1)),
	 with q_c parameter lying in the interval q_rif.
	'''
	
	if prec == None:
		prec = X.stdPrecision;
	RIFprec = RealIntervalField(prec);
	q_rif = RIFprec(q_rif);

	Q_B, midpoint_B = ellipsoid_around_logU_inBcoordinates(q_rif,sigmaC,X,prec=prec);
	entries = Q_B.list() + midpoint_B.list();
	if not contains_only_finite_rifs(entries):
		if verbose:
			print "Ellipsoid is not well-defined (some entries are infinite or NaN).";
		return False,NaN;
	lambda_low,lambda_up = lowerAndUpperBoundsForSmallestAndLargestEigenvalueOfSelfadjointRealMatrix(Q_B);
	if not lambda_low>0:
		if verbose:
			print "Matrix defining the ellipsoid is possibly not positive definite, lambda_low =",lambda_low;
		return False,NaN;		
	#if Infinity in Q_B.list() or NaN in Q_B.list():
	#	print "Ellipsoid defining Q_B contains Infinity or NaN.";
	#	return False,NaN;
	#print "Ellipsoid Q_B:\n",Q_B
	if verbose:
		print "volume of ellipsoid Q_B:",volumeOfEllipsoid(Q_B,1).upper().numerical_approx(digits=10);
	normSqBound = RIFprec(1);
	[bufferWasBigEnough,numCandidates] = My_FinckePohst_ViaGramMatrix_rif(Q_B,normSqBound,center=midpoint_B,solutionList=candidatesB,maxNumSolutions=maxCandidates,finalTransformation=None,callback=None,callbackArgs=None,lllReduce=True,breakSymmetry=False);
	return bufferWasBigEnough, numCandidates;

########################################################################
### Tests: #############################################################
########################################################################

def g_a_function(X,a_in_Mp,tau_re,tau_im,prec,n_terms=5):
	'''
	tau_re is a fixed RR or RIF.
	tau_im is the variable here!
	The function has coefficients in RIF and CIF.
	It includes already the error term!
	(A large value of n_terms will make the product expansion longer,
	 but this may also increase the numerical error for the RIFs and CIFs.)
	'''
	
	#TODO!!!!
	
	RIFprec = RealIntervalField(prec);
	CIFprec = ComplexIntervalField(prec);

	#q = CIFprec(exp(2*pi*tau));
	p = X.p;
	a1,a2 = lift_of_Mp_to_ZZ2overP(a_in_Mp,p);

	tau(tau_im) = tau_re + I*tau_im;

	ga(tau_im) = CIFprec(-1);
	#l_a = bernoulli_polynomial(a1,2)/2;

	#print "debug"

	l_a = X.ell_a[a_in_Mp];
	#print "l_a =",l_a;
	#The following should be the correct choice for q^la:
	ga *= exp(CIFprec(2*pi*l_a*I)*tau);


	ga *= CIFprec(exp(pi*I*a2*(a1-1)));

	for n in range(n_terms):
		ga *= 1-exp(CIFprec(2*pi*I*(n+a1))*tau + CIF(2*pi*I*a2));
		ga *= 1-exp(CIFprec(2*pi*I*(n+1-a1))*tau - CIF(2*pi*I*a2));

	print "TODO: Error from taking only finite product is not included yet!";

	return ga;

def g_a(X,a_in_Mp,tau,n_terms=5):
	#tau is in CIF
	prec = tau.prec();
	RIFprec = RealIntervalField(prec);
	CIFprec = ComplexIntervalField(prec);


	#q = CIFprec(exp(2*pi*tau));
	p = X.p;
	a1,a2 = lift_of_Mp_to_ZZ2overP(a_in_Mp,p);

	ga = CIFprec(-1);
	l_a = bernoulli_polynomial(a1,2)/2;
	#print "l_a =",l_a;
	#The following should be the correct choice for q^la:
	ga *= CIFprec(exp(2*pi*I*tau*l_a));

	ga *= CIFprec(exp(pi*I*a2*(a1-1)));

	for n in range(n_terms):
		ga *= CIFprec(1-exp((2*pi*I*tau)*(n+a1) + 2*pi*I*a2));
		ga *= CIFprec(1-exp((2*pi*I*tau)*(n+1-a1) - 2*pi*I*a2));

	return ga;

def vectorU_at_tau(X,sigmaC,tau):
	prec = tau.prec();
	RIFprec = RealIntervalField(prec);
	CIFprec = ComplexIntervalField(prec);

	result = [CIFprec(1) for kIndex in range(X.d)];

	#print "sigmaC =",sigmaC,", tau =",tau;

	for kIndex in range(X.d):
		#print "kIndex =",kIndex;
		k = X.GalK_lift[kIndex];
		sigma_k = X.lift_GalKIndex_to_G[kIndex];

		for a0_in_Mp in X.Orbit0:
			#a_in_Mp = X.rightAction_GH_on_Mp.action(sigmaC*sigma_k,a0_in_Mp); #Wrong order of sigmas!
			a_in_Mp = X.rightAction_GH_on_Mp.action(sigma_k*sigmaC,a0_in_Mp);
			ga = g_a(X,a_in_Mp,tau);
			#print "a_in_Mp =",a_in_Mp,", ga =",ga;
			result[kIndex] *= CIF(ga ^ X.m);
			
	return result;

def test_program_at_certain_js(X):

	def test_at_j(j,X):
		print "Test j =",j;
		prec = X.stdPrecision;
		RIFprec = RealIntervalField(prec);
		CIFprec = ComplexIntervalField(prec);
		E = EllipticCurve_from_j(j);
		print "E =",E;
		print "E has cm:",E.has_cm();
		rep = E.galois_representation();
		print "image type at p:",rep.image_type(X.p);
		tau = myCCtoCIF(E.period_lattice().tau(prec));
		#print "tau =",tau;
		q_cif = CIFprec(exp(2*pi*i*tau));
		#print "q_cif =",q_cif;
		q_rif = q_cif.real();
		print "q_rif =",q_rif;

		#t = var("t"); #the imaginary part of tau.
		#tau_re = tau.real();
		#tau_im = tau.imag();
		#a_in_Mp = X.Mp[0];
		#ga_fct = g_a_function(X,a_in_Mp,tau_re,t,prec,n_terms=2);
		##ga_fct = abs(ga_fct);
		#print "ga_fct:",ga_fct;
		#print ga_fct.derivative(t)(tau_im);

		#return;

		for sigmaC in X.Sigma:
			indexSigma = X.Sigma.index(sigmaC);
			#print "sigmaC:",sigmaC.list(),"at index",indexSigma;


			if False: #The following uses a fresh implementation for debugging purposes:
				vectorU = vectorU_at_tau(X,sigmaC,tau);
				#print "vectorU =",vectorU;
				vLogAbsU = [log(abs(x)) for x in vectorU];
				#print "vLogAbsU =",vLogAbsU;
				vLogAbsU_B = X.Alpha[prec] * vector([x for x in vLogAbsU]);
				print "vLogAbsU_B =",vLogAbsU_B;

				#prodUsigma = prod(vectorU);
				#print "prod U^sigma =", prodUsigma;
				#print "|prod U^sigma| =", abs(prodUsigma);
				#print "arg(prod U^sigma)/(2pi) =", CIFprec(prodUsigma.arg()/(2*pi));
						
			if False:
				Q_B,midpoint_B = ellipsoid_around_logU_inBcoordinates(q_rif,sigmaC,X,prec=prec,verbose=False,putMidpointErrorIntoEllipsoid=False);
				print "midpoint_B at sigmaC:",midpoint_B;
				print "volume of error ellipsoid:",volumeOfEllipsoid(Q_B,1);
				candidatesB = [];
				bufferWasBigEnough, numCandidates = findCandidatesB_inQinterval_at_sigmaC(q_rif,sigmaC,X,candidatesB,1,prec=prec,verbose=False);
				if numCandidates >= 1:
					print "Integral point is quite possibly here! With B =",candidatesB;
				else:
					print "Integral point corresponding to j is NOT in this sheet.";

			if True:
				logU = logU_asCoordinatewiseRif(q_rif,sigmaC,X,prec=None,n_summands=1);
				B = X.Alpha[X.stdPrecision] * vector(logU);
				#print "B =",B;
				if all([integralPoints_in_realInterval(b) != [] for b in B]):
					print "B =",B,"could be integral for sigmaC =",sigmaC,"with index",indexSigma;

	print "Test certain js:";

	for kIndex in range(X.d):
		#print "kIndex =",kIndex;
		k = X.GalK_lift[kIndex];
		sigma_k = X.lift_GalKIndex_to_G[kIndex];
		k_in_Zp = X.elementOfGalQzeta_AsElementIn_ZmodPstar(k);
		print "sigma_k =",sigma_k.list(),", det(sigma_k) =",det(sigma_k),", kIndex =",kIndex,", k_in_Zp =",k_in_Zp;

		#for a0_in_Mp in X.Orbit0:
		#	#a_in_Mp = X.rightAction_GH_on_Mp.action(sigmaC*sigma_k,a0_in_Mp); #Wrong order!!!
		#	a_in_Mp = X.rightAction_GH_on_Mp.action(sigma_k*sigmaC,a0_in_Mp);
		#	ga = g_a(X,a_in_Mp,tau);
		#	#print "a_in_Mp =",a_in_Mp,", ga =",ga;
		#	result[kIndex] *= CIF(ga ^ X.m);
			
	#raise Exception("Test finished!");


	CM = cm_j_invariants_and_orders(QQ);
	for (discriminant,conductor,j) in CM:
		D = Mod(discriminant,X.p);
		if D.is_zero():
			continue;
		if D.is_square():
			continue;
		print " ================================ ";
		print "j =",j,"(discr ="+str(discriminant)+") should give a point on X_ns^+(p)!";
		test_at_j(j,X);
		
	
	#j163 = -2^18*3^3*5^3*23^3*29^3;
	#test_at_j(j163,X);
	#j4 = 2^6*3^3;
	#test_at_j(j4,X);
	#j3 = 0;
	#test_at_j(j3,X);

	#For p = 7:
	#j4 = 2^6*3^3;
	#test_at_j(j4,X);
	#j8 = 2^6*5^3;
	#test_at_j(j8,X);
	raise Exception("Test finished!");

def testProductFormula(X):
	print "Test prod_{sigma in G/G_H} U^sigma:";
	prec = X.stdPrecision;
	RIFprec = RealIntervalField(prec);
	CIFprec = ComplexIntervalField(prec);
	tau = CIFprec(pi+I); #Take a random tau!!! (as it should work for any tau)
	product = CIFprec(1);
	exponent = 12*X.p;
	exponent = X.m;
	for a in X.Mp:
		ga = g_a(X,a,tau);
		product *= ga^exponent;
	print "Product of g_a^(12p) over all a in Mp:",product;
	print "p^(12p):",X.p^(CIFprec(exponent));
	print "The last two lines should give the same result."
	#Test worked!
	raise Exception("Test finished!");

def test_group_actions(X):

	print "Orbit0:",X.Orbit0,"has length",len(X.Orbit0);
	print "Length of Mp:",len(X.Mp);

	list_of_a = [];

	sigmaC = X.Sigma[1];
	
	for kIndex in range(X.d):
		#print "kIndex =",kIndex;
		k = X.GalK_lift[kIndex];
		sigma_k = X.lift_GalKIndex_to_G[kIndex];

		for a0_in_Mp in X.Orbit0:
			#a_in_Mp = X.rightAction_GH_on_Mp.action(sigmaC*sigma_k,a0_in_Mp); #This is wrong (was previously implemeneted like this...)!
			a_in_Mp = X.rightAction_GH_on_Mp.action(sigma_k*sigmaC,a0_in_Mp);
			list_of_a.append(a_in_Mp);
			
			#ga = g_a(X,a_in_Mp,tau);
			#print "a_in_Mp =",a_in_Mp,", ga =",ga;
			#result[kIndex] *= CIF(ga ^ X.m);

	print "length of list_of_a:",len(list_of_a);
	list_of_a.sort();
	print "list_of_a:",list_of_a;

	raise Exception("Test finished!");
	
def testRunningTimes():
	runningTimePer24CPUs = dict([(7,338),(11,420),(13,1397),(17,1389),(19,5579),(23,4734),(29,12332),(31,37204),(37,78528),(41,52444),(43,147720),(47,102184),(53,216312),(59,151221*2),(61,108212*8),(67,82131*16),(71,72397*16),(73,0),(79,0),(83,0),(89,0),(97,0)]);
	for p,t in runningTimePer24CPUs.iteritems():
		print p, t^(1/4.0)*2.1

########################################################################
### Main program: ######################################################
########################################################################

#@parallel(p_iter=parallelIterator,ncpus=numCPUs)
@parallel(p_iter=parallelIterator,ncpus=numCPUs)
def find_all_jCandidates_at_sigmaC_or_check_them_in_given_set(X,sigmaC,q_MaxNeg_sieve,q_MaxPos_sieve,js_to_check,verbose=False):
	'''
	There are two types of procedures unified in this method due to parallelization:
	 EITHER sigmaC is given, and we do Baker bound reduction and Fincke-Pohsting,
	 OR js_to_check is given (and sigmaC=None), and we directly check the image type of all j in this set.
	'''
	
	global pathData;
	
	debug = True;

	if sigmaC == None:
		#We are now to check all j in js_to_check:

		####################################################################
		### Extra search over 0 <= j <= 1728 and j in js_to_check: #########
		####################################################################

		if verbose:
			print "Extra search:"
		output_js = "";
		for j in js_to_check:
			#print "Check j =",j,
			#print ".",;
			if verbose:
				if j % 1000 == 0:
					print j,;
			output_j = test_j_and_p(j,X.p);
			if verbose:
				if output_j != "":
					print output_j,;
			if output_j.count("non-split") >= 1:
				output_js += output_j;

		if debug and len(js_to_check)!=0:
			filename = "extraSearch_maxJs"+str(max(js_to_check));
			if not pathData.endswith("/"):
				pathData += "/";
			path = pathData +"parts_"+str(X.p)+"/"
			if X.numParts > 1:
				path += "part_"+str(X.part)+"_of_"+str(X.numParts)+"/";
			filename = path + filename;
			if not os.path.exists(os.path.dirname(filename)):
				os.makedirs(os.path.dirname(filename))
			print "Write debug message in file",filename,"...";
			#save(output_js,filename);
			#Furthermore save the solutions to a text file:
			out = file(filename+'.txt','w')
			#out.write("###\n");
			out.write("Finished for js_to_check = "+str(js_to_check)+".");
			out.write("output_js:"+str(output_js));
			out.close();
			print "debug message written to file.";
				
		return output_js;

	#We are to do the Baker reduction and Fincke-Pohsting:	

	if verbose:
		print "=======================";
		print "sigmaC =",sigmaC.list(),", number",X.Sigma.index(sigmaC),"out of",len(X.Sigma);

	t0_sigmaC = cputime();

	####################################################################
	### Reduction of Baker bound: ######################################
	####################################################################

	if verbose:
		print "Reduce Baker bound:"
	#q_BakerBound = 10^(-100); #A stupid guess for what Baker might yield. Needs to be replaced by the actual Baker bound!!!

	B0 = X.B0[X.stdPrecision];
	q_BakerBound = RIF(0);
	#print "B0:",B0;
	#print "log q_BakerBound:",log(q_BakerBound);
	while True:
		B0_new, q_BakerBound_new = X.reduction_of_BakerBound_B0_and_qc_at_sigmaC(sigmaC,B0);
		q_BakerBound = max(q_BakerBound,q_BakerBound_new);
		if not B0_new < 0.99*B0:
			B0 = B0_new;
			break;
		B0 = B0_new;
		#print "B0:",B0;
		#print "log q_BakerBound:",log(q_BakerBound);

	#raise Exception("Done with Baker bound reduction.");
	
	#That is we need to search along the two intervals
	# [q_MaxNeg_sieve,-q_BakerBound] and [q_BakerBound,q_MaxPos_sieve].
	# (Plus the additonal j in the range in a small range
	#  depending on q_MaxNeg_sieve and q_MaxPos_sieve, which is not done in parallel.)


	####################################################################
	### Do the Fincke-Pohst: ###########################################
	####################################################################

	prec = X.stdPrecision;
	precIncrement = 200;
	RIFprec = RealIntervalField(prec);

	q_posInitialBound = q_BakerBound;
	q_negInitialBound = -q_BakerBound;

	if verbose:
		print "q_BakerBound:",q_BakerBound;

	#The following is done in a loop in order to be able to rerun it if prec is not enough:
	while True:
		if verbose:
			print "Start Fincke-Pohst with precision",prec,"...";

		js_to_check_at_sigmaC = set([]);

		#Parameters for the sieve:
		#max_j_interval_length = 10; #Determines the point at which we turn from Fincke-Pohsting to extra search.
		max_j_interval_length = 1; #Determines the point at which we turn from Fincke-Pohsting to extra search.

		#if sigmaC != matrix(X.ZmodP,2,2,[1,-2,0,1]):
		#	continue;

		numTimesPrecisionTooLow = 0;
		maxTimesPrecisionTooLow = 10; #In case this many times precision was not enough, we restart computation completely. (Otherwise there might be a dead-lock.)

		remainingIntervals = [];
		#print "q_BakerBound =",q_BakerBound;
		#print "q_MaxNeg_sieve =",q_MaxNeg_sieve;
		#print "q_MaxPos_sieve =",q_MaxPos_sieve;
		remainingIntervals.append(RIFprec(q_posInitialBound.lower(),q_MaxPos_sieve));
		remainingIntervals.append(RIFprec(q_MaxNeg_sieve,q_negInitialBound.upper()));
		#print "remainingIntervals =",remainingIntervals;

		#The following loop works through all remaining intervals
		# thereby possibly splitting some remaining intervals into smaller pieces:
		while True:
			remainingIntervals.sort(key = lambda x: -x.abs().lower());
			if remainingIntervals == []:
				finishedFinckePohstAtThisSigmaC = True;
				break;
			if verbose:
				print "Remaining intervals: "+str(len(remainingIntervals))+".",;
			#for q in remainingIntervals:
			#	print myRIFtoString(q);
			#raise Exception();
		
			q_rif = remainingIntervals.pop();
			
			#print "Now try to improve on q =",(q_rif.lower().sign()*q_rif.abs().lower()).numerical_approx(digits=10);
			
			if verbose:
				print "Check interval q_rif = "+rifSignToString(q_rif)+myRIFtoString(q_rif.abs());
			
			candidatesB = [];
			maxCandidates = 1; #We allow one candidate, because this cannot be avoided.
			bufferWasBigEnough, numCandidates = findCandidatesB_inQinterval_at_sigmaC(q_rif,sigmaC,X,candidatesB,maxCandidates=maxCandidates,prec=prec);

			if debug:
				if bufferWasBigEnough:
					if candidatesB != []:
						print "candidatesB =",candidatesB;

			if bufferWasBigEnough and bool(numCandidates == 0):
				#No candidates for the current interval q_rif, which is good!
				#print "No candidate for current interval!";
				continue;
				
			if bufferWasBigEnough and numCandidates == 1:
				B = candidatesB[0];
				if (q_rif.upper()/q_rif.lower()).log().abs() < 0.1:
					print "Only one candidate remaining in a short interval, so try to identify corresponding js:";
					
					prec_bisection = X.stdPrecision;
					RIFprec_bisection = RealIntervalField(prec_bisection);
					CIFprec_bisection = ComplexIntervalField(prec_bisection);

					#Q,midpoint = ellipsoid_around_logU(q_rif,sigmaC,X,prec=prec_bisection,verbose=False,putMidpointErrorIntoEllipsoid=True);
					#print "midpoint:",midpoint;
					print "logU:",logU_asCoordinatewiseRif(q_rif,sigmaC,X,prec=prec_bisection,verbose=False,coordinates="all");

					DlogU = DlogU_asCoordinatewiseRif(q_rif,sigmaC,X);
					print "DlogU:",DlogU;
					monotone_coordinates = [k for k in range(X.d) if DlogU[k].contains_zero()==False];
					if monotone_coordinates != []:
						k = monotone_coordinates[0];
						#In at least the coordinate k, DlogU is strictly monotone along the interval q_rif.
						print "DlogU is strictly monotone along the interval q_rif in coordinate k =",k;

						sign = DlogU[k].upper().sign();
						print "sign =",sign;

						a = RIFprec_bisection(q_rif.lower());
						b = RIFprec_bisection(q_rif.upper());

						if False and debug:
							print "Estimated derivative DlogU:";
							estimatedDlogU = (vector(logU_asCoordinatewiseRif(b,sigmaC,X,prec=prec_bisection))-vector(logU_asCoordinatewiseRif(a,sigmaC,X,prec=prec_bisection))) / (b-a);
							print estimatedDlogU;

						fa = sign * logU_asCoordinatewiseRif(a,sigmaC,X,prec=prec_bisection,coordinates=[k])[0];
						fb = sign * logU_asCoordinatewiseRif(b,sigmaC,X,prec=prec_bisection,coordinates=[k])[0];
						if fa > fb:
							raise Exception("Error detected: DlogU[k] has wrong sign, thus DlogU or logU must be implemented with at least one error!");

						f_at_B = sign * (X.Eta[X.stdPrecision] * vector(B))[k];
						print "fa =",fa;
						print "fb =",fb;
						print "f_at_B =",f_at_B;
						if f_at_B < fa or fb < f_at_B:
							print "f_at_B is not attained at the interval q_rif";
							continue;

						n_summands = 2;

						while True: #loop over precision

							if not X.Eta.has_key(prec_bisection):
								X.init_for_higher_precision(prec=prec_bisection);
							Eta = X.Eta[prec_bisection];
							
							logU_at_B = Eta * vector(B);
							print "logU_at_B:",logU_at_B;
							logUk_at_B = logU_at_B[k];
							print "logUk_at_B:",logUk_at_B;
							f_at_B = sign * logUk_at_B;

							print "a =",a;
							print "b =",b;
							print "f_at_B =",f_at_B;

							fa = sign*logU_asCoordinatewiseRif(a,sigmaC,X,prec=prec_bisection,n_summands=n_summands,coordinates=[k])[0];
							fb = sign*logU_asCoordinatewiseRif(b,sigmaC,X,prec=prec_bisection,n_summands=n_summands,coordinates=[k])[0];
							f_at_B = sign * (X.Eta[prec_bisection] * vector(B))[k];
							print "fa =",fa;
							print "fb =",fb;
							print "f_at_B =",f_at_B;
						
							
							fmid_hit_fatB_already = False;

							for loop_index in range(1000): #don't make a loop that goes on forever...
								mid = (a+b)/2;
								if (not (a<mid)) or not (mid<b):
									break;
								print "mid =",mid;
								fmid = sign*logU_asCoordinatewiseRif(mid,sigmaC,X,n_summands=n_summands,prec=prec_bisection,coordinates=[k])[0];
								print "fmid =",fmid;

								if fmid < f_at_B:
									a = RIFprec_bisection(mid.upper());
									#a_optimal = a;
								elif fmid > f_at_B:
									b = RIFprec_bisection(mid.lower());
									#b_optimal = b;
								elif fmid_hit_fatB_already:
									break;
								else:
									fmid_hit_fatB_already = True;
									
							print "a =",a;
							print "b =",b;

							j_rif = X.j_invariant_in_qInterval(RIFprec_bisection(a,b));
							print "j_rif =",j_rif;

							if j_rif.absolute_diameter() <= 0.001:
								print "j_rif is precise enough, can stop approximation now!";
								break;
							else:
								print "j_rif is too unprecise! Need to continue approximation with higher precision.";
								prec_bisection += precIncrement;
								print "New precision during approximation:",prec_bisection;
								RIFprec_bisection = RealIntervalField(prec_bisection);
								CIFprec_bisection = ComplexIntervalField(prec_bisection);
								a = RIFprec_bisection(a);
								b = RIFprec_bisection(b);
								n_summands += 2;
								print "Now take n_summands =",n_summands;
								continue;

						js_candidates = integralPoints_in_realInterval(j_rif);
						print "js_candidates =",js_candidates;

						#Test whether all coordinates of logU are (possibly) the ones we expect for candidate B:
						for j in js_candidates:
							print "j =",j,":";
							E = EllipticCurve_from_j(j);
							print "E =",E;
							tau = myCCtoCIF(E.period_lattice().tau(prec=prec_bisection));
							print "tau =",tau;
							q = RIFprec_bisection(CIFprec_bisection(exp(2*I*pi*tau)).real());
							print "q =",q;
							logU_at_j = logU_asCoordinatewiseRif(q,sigmaC,X,n_summands=2,prec=prec_bisection);
							print "logU_at_j =",logU_at_j;
							#Check whether at all coordinates, the intervals logU_at_B and logU_at_j overlap:
							if all([not (logU_at_B[k] != logU_at_j[k]) for k in range(X.d)]):
								print "All coordinates for j =",j,"match, so this j is a strong candidate!";
								js_to_check_at_sigmaC.add(j);
						
						#print "New j's to check:",integralPoints_in_realInterval(j_rif);
						#js_to_check_at_sigmaC.update(set(integralPoints_in_realInterval(j_rif)));

						continue;

						#print "DlogU:",DlogU_asCoordinatewiseRif(q_rif,sigmaC,X,prec=prec_bisection,verbose=False,coordinates="all");
				
				
			if bufferWasBigEnough and bool(numCandidates <= 1):
				#if (1/q_rif).diameter() <= 1:
				j_rif = X.j_invariant_in_qInterval(q_rif);
				if j_rif.absolute_diameter() <= max_j_interval_length:
					if verbose or debug:
						print "TODO: Found a small q-interval",q_rif,"which gives rise to at most",max_j_interval_length,"candidates, so check all relevant j's directly!";
						print "j_rif =",myRIFtoString(j_rif);
						print "candidatesB =",candidatesB
					#Sanity check:
					if len(integralPoints_in_realInterval(j_rif))>max_j_interval_length+1:
						print integralPoints_in_realInterval(j_rif);
						raise Exception("More js than expected!");
					js_to_check_at_sigmaC.update(set(integralPoints_in_realInterval(j_rif)));
					continue;
			if bufferWasBigEnough == False and numCandidates == NaN:
				if verbose:
					print "Precision was not high enough."
				numTimesPrecisionTooLow += 1;
				if numTimesPrecisionTooLow >= maxTimesPrecisionTooLow:
					msg = "Precision was too low already "+str(numTimesPrecisionTooLow)+" times.";
					msg += "\n It seems reasonable to stop trying at this point, so we do that.";
					msg += "\n We rerun the program at this sigma_c with higher precision.";
					#raise PrecisionError(msg);
					if verbose:
						print msg;

					finishedFinckePohstAtThisSigmaC = False;
					prec += precIncrement;
					RIFprec = RealIntervalField(prec);
					X.init_for_higher_precision(prec);

					if False and debug:
						print "Remaining intervals:";
						for remainingInterval in remainingIntervals:
							print myRIFtoString(remainingInterval);
						print "Current interval q_rif:",myRIFtoString(q_rif);

					#Update intial bounds by what was computed already:
					q_posInitialBound = RIFprec(q_MaxPos_sieve);
					q_negInitialBound = RIFprec(q_MaxNeg_sieve);
					for remainingInterval in remainingIntervals + [q_rif]:
						if False and debug:
							print "Remaining interval:",myRIFtoString(remainingInterval);
						if remainingInterval > 0:
							q_posInitialBound = q_posInitialBound.min(RIFprec(remainingInterval));
						elif remainingInterval < 0:
							q_negInitialBound = q_negInitialBound.max(RIFprec(remainingInterval));
						else:
							raise Exception("Some remainingInterval contains 0!");

					#In case of rounding errors, let's make sure we start at least from Baker's bound:
					q_posInitialBound = q_posInitialBound.max(RIFprec(q_BakerBound));
					q_negInitialBound = q_negInitialBound.min(-RIFprec(q_BakerBound));

					#Just to avoid confusion:
					q_posInitialBound = RIFprec(q_posInitialBound.lower());
					q_negInitialBound = RIFprec(q_negInitialBound.upper());					

					if debug:
						print "New q_posInitialBound:",myRIFtoString(q_posInitialBound);
						print "New q_negInitialBound:",myRIFtoString(q_negInitialBound);
					
					break;
					
			#Split q_rif into two smaller pieces.
			#print "Split interval into smaller pieces and retry.";
			q_split = q_rif.center().sign()*(1/2*log(RIFprec(q_rif.abs().lower())) + 1/2*log(RIFprec(q_rif.abs().upper()))).exp().center();
			#print log(RIFprec(q_rif.lower()));
			#print log(RIFprec(q_rif.upper()));
			#print (1/2*log(RIF(q_rif.lower())) + 1/2*log(RIF(q_rif.upper())));

			#print "q_split =",q_split.numerical_approx(digits=10);

			q_rif_low = RIFprec(q_rif.lower(),q_split);
			q_rif_up = RIFprec(q_split,q_rif.upper());
			remainingIntervals.append(q_rif_low);
			remainingIntervals.append(q_rif_up);
	
			#print "bufferWasBigEnough =",bufferWasBigEnough;
			#print "numCandidates =",numCandidates;
			#print "candidates =",candidatesB;

		if finishedFinckePohstAtThisSigmaC:
			break;
		else:
			continue; #With higher precision;

	if verbose or debug:	
		print "js_to_check at sigma_C:",js_to_check_at_sigmaC;
		print "Time taken for this sigma_c:",cputime(t0_sigmaC);

	if debug:
		filename = "sigma"+str(X.Sigma.index(sigmaC));
		if not pathData.endswith("/"):
			pathData += "/";
		path = pathData +"parts_"+str(X.p)+"/"
		if X.numParts > 1:
			path += "part_"+str(X.part)+"_of_"+str(X.numParts)+"/";
		filename = path + filename;
		if not os.path.exists(os.path.dirname(filename)):
			os.makedirs(os.path.dirname(filename))
		print "Write debug message in file",filename,"...";
		#save(output_js,filename);
		#Furthermore save the solutions to a text file:
		out = file(filename+'.txt','w')
		#out.write("###\n");
		out.write("Finished for sigmaC = "+str(sigmaC)+".");
		out.write("js_to_check_at_sigmaC:"+str(js_to_check_at_sigmaC));
		out.close();
		print "debug message written to file.";

	return js_to_check_at_sigmaC;
	
def integralPoints_on_XnsPlus_P(p,d=None,saveToFile = True,part=0,numParts=1,extraSearch=2^12):
	'''
	Computes the integral points on the modular curve X_{ns}^+(p).
	
	INPUT:
	 - A rational prime p >= 7, which determines the modular curve X_{ns}^+(p).
	 - A divisor d of (p-1)/2. (Default: d=None, meaning that d=
	   (p-1)/2 will be chosen automatically.) The algorithm will work
	   over the degree d subfield of the real cyclotomic field Q(Re(zeta_p)).
	 - A positive integer 'extraSearch' that determines the upper bound for j-invariants
	   that the program tests directly in the extra Search.
	   The larger p, the larger 'extraSearch' should be chosen, as for small j and large p,
	   the sieve becomes slower and slower, where as the extra search runs in almost
	   constant time for some fixed j.
	 - 0 <= part < numParts: The program cuts the modular curve into 'numParts' parts,
	   and only computes integral points in the part with index 'part'.
	   (Default: part=0, numParts=1, thus all integral points are computed.)
	   
	OUTPUT:
	  A string that lists the j-invariants of the non-cusp integral points of X_{ns}^+(p).
	'''
	
	global pathData;
	
	print "p =",p;

	print "Run part",part,"of",numParts;

	t00 = walltime();

	if d == None:
		d = ZZ((p-1)/2);

	#if d <= 2:
	#	return ValueError("d must be at least 3, otherwise we won't have a Baker bound for B0. Therefore also we must have p>=7.");

	####################################################################
	### Initialize curve and many related objects: #####################
	####################################################################

	X = X_ns_plus(p,d=d);

	X.part = part;
	X.numParts = numParts;

	#X.showMemoryUsage();

	#print "Debugging:";
	#test_group_actions(X);
	#testProductFormula(X);
	#test_program_at_certain_js(X);
	#find_all_jCandidates_at_sigmaC_or_check_them_in_given_set(X,sigmaC=X.Sigma[14],q_MaxNeg_sieve=-2^-12,q_MaxPos_sieve=2^-12,js_to_check=None,verbose=True);
	#raise Exception("End of debugging.");

	####################################################################
	### For each sigmaC, reduce Baker bound and start Fincke-Pohst: ####
	####################################################################

	q_MaxNeg = RIF(-exp(-pi*sqrt(3))).upper();
	q_MaxPos = RIF(exp(-2*pi)).upper();
	#q_MaxNeg_sieve = -2^-12;
	#q_MaxPos_sieve = 2^-12;
	q_MaxNeg_sieve = -1/extraSearch; #2^-12;
	q_MaxPos_sieve = 1/extraSearch; #2^-12;

	#We have to consider all q in the interval [q_MaxNeg,q_MaxPos].
	#The q in (-q_BakerBound,+q_BakerBound) are excluded by Baker method.
	#The q in [q_MaxNeg_sieve,-q_BakerBound] and [q_BakerBound,q_MaxPos_sieve]
	# are dealth with in sieve (reduction and ellipsoids).
	#The q in [q_MaxNeg,q_MaxNeg_sieve] and [q_MaxPos_sieve,q_MaxPos]
	# are dealt with in extra search.

	js_to_check_apriori = set([]);
	js_to_check_apriori.update(set(integralPoints_in_realInterval(X.j_invariant_in_qInterval(RIF(q_MaxNeg,q_MaxNeg_sieve)))));
	js_to_check_apriori.update(set(integralPoints_in_realInterval(X.j_invariant_in_qInterval(RIF(q_MaxPos_sieve,q_MaxPos)))));
	js_to_check_apriori.update(set(range(1728+1)));

	js_to_check_remaining = set([]);

	output_js = "";

	verboseDuringParallelInstances = True;

	#print "Test sieve:"
	#print find_all_jCandidates_at_sigmaC_or_check_them_in_given_set(X,X.Sigma[0],q_MaxNeg_sieve,q_MaxPos_sieve,None,True);
	#raise Exception("Test finished!");

	#print "Test extra search:"
	#print find_all_jCandidates_at_sigmaC_or_check_them_in_given_set(X,None,None,None,js_to_check_apriori,True);
	#raise Exception("Test finished!");

	print "Start parallel computing...";
	parameters = [];
	paramJsApriori = (X,None,None,None,js_to_check_apriori,verboseDuringParallelInstances);

	for sigmaC in X.Sigma:
		indexSigma = X.Sigma.index(sigmaC);

		##Debug:
		#if indexSigma != 435: #80: #110: #14
		#	continue;

		if not mod(indexSigma-part,numParts).is_zero():
			continue;

		paramSigmaC = (X,sigmaC,q_MaxNeg_sieve,q_MaxPos_sieve,None,True);
		parameters.append(paramSigmaC);

	#Let the last part of the program do the "extra search" for small j's:
	if mod(part+1,numParts).is_zero():
		parameters.append(paramJsApriori);

	gen = find_all_jCandidates_at_sigmaC_or_check_them_in_given_set(parameters);
	for x in gen:
		sigmaC = x[0][0][1];
		output = x[1];
		print "x:",x;

		##Debug:
		#for parameter in parameters:
		#	sigmaC = parameter[1];
		#	output = find_all_jCandidates_at_sigmaC_or_check_them_in_given_set(*parameter);
	
		print "Parallel instance returned:",output;
		print sigmaC;
		if sigmaC == None:
			#js have been checked:
			output_js += output; #x[1];
			print "js_to_check_apriori where checked!";
			print "output_js:",output_js;
		else:
			#Sieving was done:
			js_to_check_for_sigmaC = output;
			print "sigmaC =",sigmaC.list(),"was checked.";
			print "js_to_check_for_sigmaC:",js_to_check_for_sigmaC;
			js_to_check_for_sigmaC.difference_update(js_to_check_apriori);
			js_to_check_remaining = js_to_check_remaining.union(js_to_check_for_sigmaC);

	#Check remaining js which came as candidates from some sigmaC:
	print "js_to_check_remaining:",js_to_check_remaining;
	output_js_remaining = find_all_jCandidates_at_sigmaC_or_check_them_in_given_set(X,None,None,None,js_to_check_remaining,verboseDuringParallelInstances);
	print "output_js_remaining:",output_js_remaining;
	output_js += output_js_remaining;

	print "==================================";
	print "Finished for p = "+str(p)+"!";
	print "Set of remaining candidates j:";
	print output_js;

	totalTime = walltime(t00);
	print "Total time:",totalTime;

	##Debug:
	#return;

	if saveToFile:
		filename = "nonsplit_"+str(X.p);
		if numParts > 1:
			filename += "_part_"+str(part)+"_of_"+str(numParts);
		if not pathData.endswith("/"):
			pathData += "/";
		filename = pathData + filename;
		if not os.path.exists(os.path.dirname(filename)):
			os.makedirs(os.path.dirname(filename))
		#save(output_js,filename);
		#Furthermore save the solutions to a text file:
		out = file(filename+'.txt','w')
		#out.write("###\n");
		out.write("# List of j-invariants of possible integral points on X_ns^+(p).\n");
		out.write("# Computing this list took "+str(ceil(totalTime))+" seconds on "+str(numCPUs)+" cpus.\n");
		out.write("# Authors: Aurelien Bajolet, Yuri Bilu, Benjamin Matschke, 2016.\n");
		out.write("# License: Creative commons 3.0 by-nc.\n");
		out.write("#\n");
		out.write(output_js);
		out.close();
	
	return True;

########################################################################
### Further tests: #####################################################
########################################################################

#integralPoints_on_XnsPlus_P(5,d=2); #First reduction fails, as expected.
#integralPoints_on_XnsPlus_P(7,d=3); #6.7 min per sigmaC
#integralPoints_on_XnsPlus_P(11); #2.5 min per sigmaC
#integralPoints_on_XnsPlus_P(11,part=0,numParts=80); #2.5 min per sigmaC
#integralPoints_on_XnsPlus_P(13); #5.5 min per sigmaC
#integralPoints_on_XnsPlus_P(17,d=4); #12 min per sigmaC
#integralPoints_on_XnsPlus_P(17); #4 min per sigmaC (increased precision at the beginning automatically to 200)
#integralPoints_on_XnsPlus_P(19); #11 min per sigmaC (increased precision at the beginning automatically to 200 resp. 300)
#integralPoints_on_XnsPlus_P(23); #5.5 min per sigmaC (increased precision at the beginning automatically to 200)
#integralPoints_on_XnsPlus_P(29);#11 min per sigmaC (increased precision at the beginning automatically to 400)
#integralPoints_on_XnsPlus_P(31);#30 min per sigmaC (increased precision at the beginning automatically to 400)
#integralPoints_on_XnsPlus_P(37);#45 min per sigmaC 
#integralPoints_on_XnsPlus_P(41);#25 min per sigmaC 
#integralPoints_on_XnsPlus_P(41);#65 min per sigmaC 

#Tests for different d's:
#integralPoints_on_XnsPlus_P(29); #481 sec for first sigma
#integralPoints_on_XnsPlus_P(29,d=7); #972 sec for first sigma
#integralPoints_on_XnsPlus_P(41); #1021 sec for first sigma
#integralPoints_on_XnsPlus_P(41,d=10); #1659 sec for first sigma
#integralPoints_on_XnsPlus_P(97); #68564 sec for first sigma
#integralPoints_on_XnsPlus_P(97,d=12); #55448 sec for first sigma
#integralPoints_on_XnsPlus_P(97,d=16); #14795 sec for first sigma
#and p=97 on plafrim: d=48 75543 sec, 

def reduceB0OldToB0New():
	#In the arxiv versions 1 and 2 of the paper we use different initial
	#bounds B0 (or rather for q).
	#The new version reduces the amount of technicalities in the paper. 
	#In the computations for p<100 we used only the old bound.
	#In order to make all computations valid also wrt. the second arxiv version,
	#we reduce the new bound to the old bound whenever necessary:
	for p in prime_range(10,100):
		print "================================== p:",p;
		X = X_ns_plus(p);
		X.reduction_from_B0old_to_B0new();

#reduceB0OldToB0New(); #Succeeds with all reductions.

'''
delta's:
11 5
13 3
17 4
19 3
23 11
29 4
31 3
37 3
41 4
43 3
47 23
53 4
59 29
61 3
67 3
71 5
73 3
79 3
83 41
89 4
97 3
'''

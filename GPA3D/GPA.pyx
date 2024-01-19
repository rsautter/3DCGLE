import numpy
from libc.math cimport pow, fabs, sqrt, M_PI, sin, cos,tan,floor
from math import radians, atan2,factorial
from scipy.spatial import Delaunay as Delaunay
from scipy.spatial import QhullError
import matplotlib.pyplot as plt

from cpython cimport bool
cimport numpy
cimport cython

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cdef class GPA:
	cdef public double[:,:] mat,gradient_dx, gradient_dy
	cdef public double cx, cy
	cdef public int rows, cols
	
	cdef public double[:,:] phases, mods
	cdef public int[:,:] symmetricalP, asymmetricalP, unknownP
	cdef public object triangulation_points,triangles
	cdef public double maxGrad, tol
	cdef public object cvet

	cdef public int n_edges, n_points
	cdef public double G1, G2, G3, G1_Classic
	cdef public object G4

	#@profile
	def __cinit__(self, double tol=0.03):
		# setting matrix,and calculating the gradient field
		self.tol = tol

		# percentual Ga proprieties
		self.symmetricalP = numpy.array([[]],dtype=numpy.int32)
		self.asymmetricalP = numpy.array([[]],dtype=numpy.int32)
		self.unknownP = numpy.array([[]],dtype=numpy.int32)
		self.triangulation_points = []

	@cython.boundscheck(False)
	@cython.wraparound(False)
	@cython.nonecheck(False)
	@cython.cdivision(True)
	cpdef void setPosition(self, double cx, double cy):
		self.cx = cx
		self.cy = cy

	@cython.boundscheck(False)
	@cython.wraparound(False)
	@cython.nonecheck(False)
	@cython.cdivision(True)
	cdef void _setGradients(self):
		cdef double[:,:] gx, gy
		gy, gx = numpy.gradient(self.mat)
		#gy, gx = self.gradient(self.mat)
		#initialization
		self.gradient_dx = gx
		self.gradient_dy = gy

		self._setMaxGrad()
		self._setModulusPhase()
		
	@cython.boundscheck(False)
	@cython.wraparound(False)
	@cython.nonecheck(False)
	@cython.cdivision(True)
	cdef void _setMaxGrad(self):
		cdef int i,j,w,h
		self.maxGrad = -1.0
		w, h = self.cols, self.rows 
		for i in range(w):
			for j in range(h):
				if(self.maxGrad<0.0) or (sqrt(pow(self.gradient_dy[j, i],2.0)+pow(self.gradient_dx[j, i],2.0))>self.maxGrad):
					self.maxGrad = sqrt(pow(self.gradient_dy[j, i],2.0)+pow(self.gradient_dx[j, i],2.0))
		if self.maxGrad < 1e-5:
			self.maxGrad = 1.0
		
	@cython.boundscheck(False)
	@cython.wraparound(False)
	@cython.nonecheck(False)
	@cython.cdivision(True)
	cdef void _setModulusPhase(self):
		cdef int w, h, i, j
		cdef double[:,:] gx, gy
		
		gx = self.gradient_dx
		gy = self.gradient_dy
		w, h = self.cols, self.rows 
		
		self.phases = numpy.array([[atan2(self.gradient_dy[j, i],self.gradient_dx[j, i]) if atan2(gy[j, i],gx[j, i])>0 else atan2(gy[j, i],gx[j, i])+2.0*M_PI
									 for i in range(w) ] for j in range(h)],dtype=numpy.float64)
		self.mods = numpy.array([[self.getMod(self.gradient_dx[j, i], self.gradient_dy[j, i]) for i in range(w) ] for j in range(h)],dtype=numpy.float64)
		
	@cython.boundscheck(False)
	@cython.wraparound(False)
	@cython.nonecheck(False)
	@cython.cdivision(True)
	cdef double getMod(self,double x, double y):
		return sqrt(pow(x,2.0)+pow(y,2.0))/self.maxGrad
		
	@cython.boundscheck(False)
	@cython.wraparound(False)
	@cython.nonecheck(False)
	@cython.cdivision(True)
	cpdef char* version(self):
		return "GPA - 3.5"
	
	@cython.boundscheck(False)
	@cython.wraparound(False)
	@cython.nonecheck(False)
	@cython.cdivision(True)
	cdef void _update_asymmetric_mat(self,double[:] index_dist,double[:,:] dists,double tol, double ptol):
		cdef int ind, lx, px, py, px2, py2, i, j
		cdef int[:] x, y
		
		self.symmetricalP = numpy.zeros((len(dists),len(dists[0])),dtype=numpy.int32)
		self.asymmetricalP = numpy.zeros((len(dists),len(dists[0])),dtype=numpy.int32)
		self.unknownP = numpy.zeros((len(dists),len(dists[0])),dtype=numpy.int32)
		
		# distances loop
		for ind in range(0, len(index_dist)):
			x2, y2 =[], []
			for py in range(self.rows):
				for px in range(self.cols):
					if (fabs(dists[py, px]-index_dist[ind]) <= fabs(ptol)):
						x2.append(px)
						y2.append(py)
			x, y = numpy.array(x2,dtype=numpy.int32), numpy.array(y2,dtype=numpy.int32)
			lx = len(x)

			# compare each point in the same distance
			for i in range(lx):
				px, py = x[i], y[i]
				
				if self.mods[py,px]<= tol:
					self.unknownP[py,px] = 1
					continue
				
				for j in range(lx):
					px2, py2 = x[j], y[j]
					if self.mods[py2,px2]<= tol:
						continue			
					dx = self.gradient_dx[py, px]+self.gradient_dx[py2, px2]
					dy = self.gradient_dy[py, px]+self.gradient_dy[py2, px2]
					# Vetores sao simetricos
					if self.getMod(dx,dy)<= tol:
						self.symmetricalP[py,px] = 1
						#print("{",py, px,"}{",py2, px2,"} -> ", self.gradient_dx[py, px], self.gradient_dy[py, px], self.gradient_dx[py2, px2], self.gradient_dy[py2, px2])
						# se outro for simetrico ele vai se marcar
						break
						
		# Caso nao seja desconhecido ou simetrico, ele eh asimmetrico
		for py in range(self.rows):
			for px in range(self.cols):
				if self.symmetricalP[py,px] == 0 and self.unknownP[py,px] ==0:
					self.asymmetricalP[py,px] = 1


	@cython.boundscheck(False)
	@cython.wraparound(False)
	@cython.nonecheck(False)
	@cython.cdivision(True)
	def _G1_Classic(self,str symm):
		cdef int w, h, i, j
		cdef int[:,:] targetMat
		cdef float ratio
		self.triangulation_points = []
		d1 = self.symmetricalP


		if symm == 'S':# Symmetrical matrix 
			targetMat = self.symmetricalP
		elif symm == 'A':# Asymmetrical matrix 
			targetMat = self.asymmetricalP
		elif symm == 'F': # Full Matrix, including unknown vectors
			targetMat = numpy.ones((self.symmetricalP.shape[0],self.symmetricalP.shape[1]),dtype=numpy.int32)
		else:
			raise Exception("Unknown analysis type (should be S,A or K), got: "+symm)
		
		for i in range(self.rows):
			for j in range(self.cols):
					if targetMat[i,j] > 0:
						self.triangulation_points.append([j+0.5*self.gradient_dx[i, j], i+0.5*self.gradient_dy[i, j]])		
		self.triangulation_points = numpy.array(self.triangulation_points)
		if len(self.triangulation_points>1):
			self.triangulation_points = numpy.unique(self.triangulation_points,axis=0)
		self.n_points = len(self.triangulation_points)
		if self.n_points < 3:
			self.n_edges = 0
			self.G1_Classic = 0.0
		else:
			try:
				self.triangles = Delaunay(self.triangulation_points)
				neigh = self.triangles.vertex_neighbor_vertices
				self.n_edges = len(neigh[1])/2
				self.G1_Classic = (float(self.n_edges)-float(self.n_points))/float(self.n_points)
			except QhullError:
				self.n_edges = 0
				self.G1_Classic = 0.0
				self.n_points = 0
		if self.G1_Classic < 0.0:
			self.G1_Classic = 0.0

	@cython.boundscheck(False)
	@cython.wraparound(False)
	@cython.nonecheck(False)
	@cython.cdivision(True)
	def _getDistancesTriang(self,points,simplices):
		ds = []
		for p in simplices:
			p1 = points[p[0]]
			p2 = points[p[1]]
			p3 = points[p[2]]
			ds.append(numpy.sqrt(numpy.sum((p1-p2)**2)))
			ds.append(numpy.sqrt(numpy.sum((p2-p3)**2)))
			ds.append(numpy.sqrt(numpy.sum((p3-p1)**2)))
		return ds	
	
	@cython.boundscheck(False)
	@cython.wraparound(False)
	@cython.nonecheck(False)
	@cython.cdivision(True)
	def _G1(self,str symm):
		cdef int w, h, i, j
		cdef int[:,:] targetMat
		cdef float ratio
		self.triangulation_points = []
		d1 = self.symmetricalP


		if symm == 'S':# Symmetrical matrix 
			targetMat = self.symmetricalP
		elif symm == 'A':# Asymmetrical matrix 
			targetMat = self.asymmetricalP
		elif symm == 'F': # Full Matrix, including unknown vectors
			targetMat = numpy.ones((self.symmetricalP.shape[0],self.symmetricalP.shape[1]),dtype=numpy.int32)
		else:
			raise Exception("Unknown analysis type (should be S,A or K), got: "+symm)
		
		for i in range(self.rows):
			for j in range(self.cols):
					if targetMat[i,j] > 0:
						self.triangulation_points.append([j+0.5*self.gradient_dx[i, j], i+0.5*self.gradient_dy[i, j]])		
		self.triangulation_points = numpy.array(self.triangulation_points)
		if len(self.triangulation_points>1):
			self.triangulation_points = numpy.unique(self.triangulation_points,axis=0)
		self.n_points = len(self.triangulation_points)
		if self.n_points < 3:
			self.n_edges = 0
			self.G1 = 0.0
		else:
			try:
				self.triangles = Delaunay(self.triangulation_points)
				neigh = self.triangles.vertex_neighbor_vertices
				self.n_edges = len(neigh[1])/2
				ds = self._getDistancesTriang(self.triangulation_points,self.triangles.simplices)
				ds = numpy.sort(ds)/numpy.max(ds)
				#print(numpy.average(ds[len(ds)//2:]),numpy.average(ds[:len(ds)//2]))
				self.G1 = (numpy.average(ds[len(ds)//2:])-numpy.average(ds[:len(ds)//2]))/numpy.max(ds)
			except QhullError:
				self.n_edges = 0
				self.G1 = 0.0
				self.n_points = 0
		if self.G1 < 0.0:
			self.G1 = 0.0

	@cython.boundscheck(False)
	@cython.wraparound(False)
	@cython.nonecheck(False)
	@cython.cdivision(True)
	cdef void _G2(self,str symm):
		cdef int i,j
		cdef double somax, somay, phase, alinhamento, mod, smod, maxEntropy
		cdef int[:,:] targetMat,opositeMat
		cdef double[:,:] probabilityMat
		somax = 0.0
		somay = 0.0
		smod = 0.0
		
		if symm == 'S':# Symmetrical matrix
			targetMat = self.symmetricalP
			opositeMat = self.asymmetricalP
		elif symm == 'A':# Asymmetrical matrix
			targetMat = self.asymmetricalP
			opositeMat = self.symmetricalP
		elif symm == 'F':# Full Matrix, including unknown vectors
			targetMat = numpy.ones((self.symmetricalP.shape[0],self.symmetricalP.shape[1]),dtype=numpy.int32)
			opositeMat = numpy.zeros((self.symmetricalP.shape[0],self.symmetricalP.shape[1]),dtype=numpy.int32)
		else:
			raise Exception("Unknown analysis type (should be S,A or F), got: "+symm)
		
		if numpy.sum(targetMat)<1:
			self.G2 = 0.0
			return
		
		alinhamento = 0.0
		
		if symm != 'S':
			for i in range(self.rows):
				for j in range(self.cols):
					if targetMat[i,j] == 1:
						somax += self.gradient_dx[i,j]/self.maxGrad
						somay += self.gradient_dy[i,j]/self.maxGrad
						smod += self.mods[i,j]
			if smod <= 0.0:
				alinhamento = 0.0
			else:
				alinhamento = sqrt(pow(somax,2.0)+pow(somay,2.0))/(2*smod)
			if numpy.sum(opositeMat)+numpy.sum(targetMat)> 0:
				self.G2 = (float(numpy.sum(targetMat))/float(numpy.sum(opositeMat)+numpy.sum(targetMat)) )*(1.0-alinhamento)
			else: 
				self.G2 = 0.0
		else:
			probabilityMat = self.mods*numpy.array(targetMat,dtype=numpy.float64)
			probabilityMat = probabilityMat/numpy.sum(probabilityMat)
			maxEntropy = numpy.log(numpy.float64(numpy.sum(targetMat)))
			for i in range(self.rows):
				for j in range(self.cols):
					if targetMat[i,j] == 1:
						alinhamento = alinhamento - probabilityMat[i,j]*numpy.log(probabilityMat[i,j])/maxEntropy
			self.G2 = alinhamento

	@cython.boundscheck(False)
	@cython.wraparound(False)
	@cython.nonecheck(False)
	@cython.cdivision(True)
	cdef double distAngle(self,double a1,double a2):
		return (cos(a1)*cos(a2)+sin(a1)*sin(a2)+1)/2

	@cython.boundscheck(False)
	@cython.wraparound(False)
	@cython.nonecheck(False)
	@cython.cdivision(True)
	cdef void _G3(self,str symm):
		cdef int x1, y1, x2, y2, i, j, div
		cdef double sumPhases, alinhamento,nterms,angle
		cdef int[:,:] targetMat,opositeMat
		cdef int[:,:] targetList

		if symm == 'S':# Symmetrical matrix 
			targetMat = self.symmetricalP
			opositeMat = self.asymmetricalP
		elif symm == 'A':# Asymmetrical matrix 
			targetMat = self.asymmetricalP
			opositeMat = self.symmetricalP
		elif symm == 'F': # Full Matrix, including unknown vectors
			targetMat = numpy.ones((self.symmetricalP.shape[0],self.symmetricalP.shape[1]),dtype=numpy.int32)
			opositeMat = numpy.zeros((self.symmetricalP.shape[0],self.symmetricalP.shape[1]),dtype=numpy.int32)
		else:
			raise Exception("Unknown analysis type (should be S,A or F), got: "+symm)
		
	
		targetList = numpy.zeros((numpy.sum(targetMat),2),dtype=numpy.int32)
		
		i = 0
		for ty in range(self.rows):
			for tx in range(self.cols):
				if targetMat[ty,tx]>0:
					targetList[i,0] = ty
					targetList[i,1] = tx
					i = i+1
		
		sumPhases = 0.0
		nterms = 0.0
		alinhamento = 0.0
		for i in range(len(targetList)):
			x1, y1 = targetList[i,0],targetList[i,1] 
			y2, x2  = x1-int(self.cx), y1-int(self.cy)
			angle = atan2(y2,x2) if atan2(y2,x2)>0 else atan2(y2,x2)+2.0*M_PI
			sumPhases += self.distAngle(self.phases[x1,y1],angle)
			nterms = nterms + 1.0
		if nterms>0.0:
			alinhamento = sumPhases / nterms
		else:
			alinhamento = 0.0
		if numpy.sum(opositeMat)+numpy.sum(targetMat)> 0:
			self.G3 = ((float(numpy.sum(targetMat))/float(numpy.sum(opositeMat)+numpy.sum(targetMat)) ) + alinhamento)/2
		else: 
			self.G3 = 0.0
			
	@cython.boundscheck(False)
	@cython.wraparound(False)
	@cython.nonecheck(False)
	@cython.cdivision(True)
	def _G4(self,str symm):
		cdef int w, h, i, j
		cdef int[:,:] targetMat
		cdef object z
		cdef float sumZ

		if symm == 'S':# Symmetrical matrix 
			targetMat = self.symmetricalP
		elif symm == 'A':# Asymmetrical matrix 
			targetMat = self.asymmetricalP
		elif symm == 'F': # Full Matrix, including unknown vectors
			targetMat = numpy.ones((self.symmetricalP.shape[0],self.symmetricalP.shape[1]),dtype=numpy.int32)
		else:
			raise Exception("Unknown analysis type (should be S,A or F ), got: "+symm)
		
		self.G4 = 0.0+0.0j
		sumZ = 0.0
		
		for i in range(self.rows):
			for j in range(self.cols):
				if targetMat[i,j] > 0:
					if self.mods[i,j] > 1e-6:
						sumZ += numpy.abs(self.mods[i,j]*numpy.exp(1j*self.phases[i,j]))
		if sumZ<1e-5:
			return self.G4 
			
		for i in range(self.rows):
			for j in range(self.cols):
				if targetMat[i,j] > 0:
					if self.mods[i,j] > 1e-6:
						z = self.mods[i,j]*numpy.exp(1j*self.phases[i,j])/sumZ
						self.G4 = self.G4 - z*numpy.log(z)
		return self.G4

	@cython.boundscheck(False)
	@cython.wraparound(False)
	@cython.nonecheck(False)
	@cython.cdivision(True)
	def __call__(self,double[:,:] mat=None, double[:,:] gx=None,double[:,:] gy=None,list moment=["G2"],str symmetrycalGrad='A'):
		if (mat is None) and (gx is None) and (gy is None):
			raise Exception("Matrix or gradient must be stated!")
		if ((gx is None) and not(gy is None)) or (not(gx is None) and (gy is None)):
			raise Exception("Gradient must have 2 components (gx and gy)")
		if not(mat is None) and not(gx is None):
			raise Exception("Matrix or gradient must be stated, not both")
		
		if not(mat is None):
			return self._eval(mat,moment,symmetrycalGrad)
		else:
			return self._evalGradient(gx,gy,moment,symmetrycalGrad)

	@cython.boundscheck(False)
	@cython.wraparound(False)
	@cython.nonecheck(False)
	@cython.cdivision(True)
	def _eval(self,double[:,:] mat,list moment=["G2"],str symmetrycalGrad='A'):
		cdef int[:] i
		cdef int x, y
		cdef double minimo, maximo
		cdef dict retorno
		
		self.mat = mat
		self.cols = len(self.mat[0])
		self.rows = len(self.mat)
		self.setPosition(float(self.rows-1)/2.0,float(self.cols-1)/2.0)
		self._setGradients()
		

		cdef numpy.ndarray dists = numpy.array([[sqrt(pow(float(x)-self.cx, 2.0)+pow(float(y)-self.cy, 2.0)) \
												  for x in range(self.cols)] for y in range(self.rows)])
		
		minimo, maximo = numpy.min(dists),numpy.max(dists)
		sequence = numpy.arange(minimo,maximo,0.705).astype(dtype=numpy.float64)
		cdef numpy.ndarray uniq = numpy.array([minimo for minimo in  sequence])
		
		# removes the symmetry in gradient_asymmetric_dx and gradient_asymmetric_dy:
		self._update_asymmetric_mat(uniq.astype(dtype=numpy.float64), dists.astype(dtype=numpy.float64), self.tol, numpy.float64(1.41))
		
		#gradient moments:
		retorno = {}
		for gmoment in moment:
			if("G4" == gmoment):
				self._G4(symmetrycalGrad)
				retorno["G4"] = self.G4
			if("G3" == gmoment):
				self._G3(symmetrycalGrad)
				retorno["G3"] = self.G3
			if("G2" == gmoment):
				self._G2(symmetrycalGrad)
				retorno["G2"] = self.G2
			if("G1" == gmoment):
				self._G1(symmetrycalGrad)
				retorno["G1"] = self.G1
			if("G1C" == gmoment):
				self._G1_Classic(symmetrycalGrad)
				retorno["G1C"] = self.G1_Classic
		return retorno
		
	def getAsymmetricalMask(self):
		return numpy.array(self.asymmetricalP)
	
	def getSymmetricalMask(self):
		return numpy.array(self.symmetricalP)
	
	def getUnknownMask(self):
		return numpy.array(self.unknownP)	
	
	def getDx(self):
		return numpy.array(self.gradient_dx)

	def getDy(self):
		return numpy.array(self.gradient_dy)
	
	@cython.boundscheck(False)
	@cython.wraparound(False)
	@cython.nonecheck(False)
	@cython.cdivision(True)
	def _evalGradient(self,double[:,:] gradient_dx, double[:,:] gradient_dy, list moment=["G2"],str symmetrycalGrad='A'):
		cdef int[:] i
		cdef int x, y
		cdef double minimo, maximo
		cdef dict retorno
		
		self.cols = len(gradient_dx[0])
		self.rows = len(gradient_dx)
		
		self.gradient_dx = gradient_dx
		self.gradient_dy = gradient_dy
		
		self._setMaxGrad()
		self._setModulusPhase()
		
		self.setPosition(float(self.rows-1)/2.0,float(self.cols-1)/2.0)
		

		cdef numpy.ndarray dists = numpy.array([[sqrt(pow(float(x)-self.cx, 2.0)+pow(float(y)-self.cy, 2.0)) \
												  for x in range(self.cols)] for y in range(self.rows)])
		
		minimo, maximo = numpy.min(dists),numpy.max(dists)
		sequence = numpy.arange(minimo,maximo,0.705).astype(dtype=numpy.float64)
		cdef numpy.ndarray uniq = numpy.array([minimo for minimo in  sequence])
		
		# removes the symmetry in gradient_asymmetric_dx and gradient_asymmetric_dy:
		self._update_asymmetric_mat(uniq.astype(dtype=numpy.float64), dists.astype(dtype=numpy.float64), self.tol, numpy.float64(1.41))
		
		#gradient moments:
		retorno = {}
		for gmoment in moment:
			if("G4" == gmoment):
				self._G4(symmetrycalGrad)
				retorno["G4"] = self.G4
			if("G3" == gmoment):
				self._G3(symmetrycalGrad)
				retorno["G3"] = self.G3
			if("G2" == gmoment):
				self._G2(symmetrycalGrad)
				retorno["G2"] = self.G2
			if("G1" == gmoment):
				self._G1(symmetrycalGrad)
				retorno["G1"] = self.G1
			if("G1C" == gmoment):
				self._G1_Classic(symmetrycalGrad)
				retorno["G1C"] = self.G1_Classic
		return retorno
   
	@cython.boundscheck(False)
	@cython.wraparound(False)
	@cython.nonecheck(False)
	@cython.cdivision(True)
	cdef tuple gradient(self,double[:,:] mat):
		cdef double[:,:] dx, dy
		cdef double divx, divy
		cdef int i, j,w,h,i1,j1,i2,j2
		w, h = len(mat[0]),len(mat)
		dx = numpy.array([[0.0 for i in range(w) ] for j in range(h)],dtype=numpy.float64)
		dy = numpy.array([[0.0 for i in range(w) ] for j in range(h)],dtype=numpy.float64)
		for i in range(w):
			for j in range(h):
				divy =  2.0 if (j<len(mat)-1 and j>0) else 1.0
				divx =  2.0 if (i<len(mat[j])-1 and i>0) else 1.0
				i1 = (i+1) if i<len(mat[j])-1 else i
				j1 = (j+1) if j<len(mat)-1 else j
				i2 = (i-1) if i>0 else i
				j2 = (j-1) if j>0 else j
				dy[j, i] = (mat[j1, i] - mat[j2, i])/divy
				dx[j, i] = (mat[j, i1] - mat[j, i2])/divx
		return dy,dx

	

	


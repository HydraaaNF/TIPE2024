import math, random, time, matplotlib.pyplot as plt, json

class Vect:
	def __init__(self, *coords):
		self.coordonnees = list(coords)

	def __getitem__(self, i):
		if i > self.dim():
			return 0
		return self.coordonnees[i]

	def __setitem__(self, key, value):
		self.coordonnees[key] = value

	def __mul__(self, v2):
		if isinstance(v2, int) or isinstance(v2, float):
			return Vect(*[i * v2 for i in self.coordonnees])
		elif isinstance(v2, Vect):
			return sum([self[i] * v2[i] for i in range(self.dim())])

	def __add__(self, v2):
		if isinstance(v2, Vect):
			return Vect(*[self[i] + v2[i] for i in range(max(self.dim(), v2.dim()))])

	def __sub__(self, v2):
		if isinstance(v2, Vect):
			return v2 * (-1) + self

	def __eq__(self, v2):
		return isinstance(v2, Vect) and all([abs(self[i] - v2[i]) <= 10**-5 for i in range(max(self.dim(), v2.dim()))])

	def __str__(self):
		return str(tuple(self.coordonnees))

	def dim(self):
		return len(self.coordonnees)

	def copy(self):
		return Vect(*[i for i in self.coordonnees])

	def norme(self):
		return (self * self) **0.5

	def normalise(self):
		n = self.norme()
		if n != 0:
			for i in range(self.dim()):
				self[i] = self[i] / n

	def normalised(self):
		clone = self.copy()
		clone.normalise()
		return clone

	@classmethod
	def vectoriel(cls, v1, v2, d=0):
		if d == 0:
			d = max(v1.dim(), v2.dim())
		return Vect(*[v1[(i+1)%d] * v2[(i+2)%d] - v1[(i+2)%d] * v2[(i+1)%d] for i in range(d)])

	@classmethod
	def get_lz(cls, A, B, P):
		return (B[0] - A[0]) * (P[1] - A[1]) - (B[1] - A[1]) * (P[0] - A[0])
	
	@classmethod
	def colli(cls, u, v):
		w = Vect.vectoriel(u, v)
		return w.norme() <=  10**-5
	
	@classmethod
	def intersect(cls, a, u, b, v):
		if cls.colli(u,v):
			return
		n1 = Vect(u[1], -u[0])
		n2 = Vect(v[1], -v[0])
		det = n1[0]*n2[1] - n2[0] * n1[1]
		c = a[0] * n1[0] + a[1] * n1[1]
		f = b[0] * n2[0] + b[1] * n2[1]
		x = (c * n2[1] - n1[1] * f) / det
		y = (f * n1[0] - c * n2[0]) / det
		return Vect(x, y)

	@classmethod
	def projeter(cls, a, b, c, d):
		"""
		Projection de CD sur AB
		"""
		m = Vect.intersect(a, b - a, c, d - c)
		if m == None:
			return
		l = abs((m - a).norme() + (m - b).norme() - (b - a).norme())
		di = (d - c) * (m - c)
		if abs((m - a).norme() + (m - b).norme() - (b - a).norme()) <= 10**-5 and (d - c) * (m - c) > 0:
			return m

	@classmethod
	def raycast(cls, a, v, segments, d=None, b=None):
		proj = []
		for i in range(len(segments)):
			p = Vect.projeter(segments[i][0], segments[i][1], a, a + v)
			if i!= b and p != None and (d==None or (p-a).norme() < d):
				proj.append((p, i))
		if proj != []:
			distance = lambda x : (x[0] - a).norme()
			proj.sort(key=distance)
			return proj[0]
		return (None, None)

	@classmethod
	def reflect(cls, a, b, c):
		"""
		Réflection d'un objetau point A sur le 'miroir' CD
		"""
		#print("Réflection...")
		#print("a b c", a, b, c)
		u = (b - c).normalised()
		#print("cb", u)
		v = u * (u * (b - a))
		#print("v", v)
		return a + v * 2

	@classmethod
	def proximite(cls, a, b , c, r):
		#print("==Vect proximite !!")
		"""e = Vect(68.75, 90.625)
								if (b - e).norme() <= r:
									print("extrémité à proximité")"""
		if (c - a).norme() <= r:
			#print("dep")
			return a
		t = Triangle(a, b, c)
		h = t.projeter_h(2)
		#print("h", h)
		re = (h - c).norme()
		#print("re", re)
		if re <= r:
			hi = (r**2 - re**2)**0.5
			i = h - (b-a).normalised() * hi
			#print((i - a).norme(), (b - a).norme())
			if (i - a).norme() <= (b - a).norme():
				#print("TTT")
				return i

class Chemin:
	def __init__(self, *points):
		self.points = list(points)

	def __getitem__(self, i):
		return self.points[i % self.taille()]

	def __str__(self):
		if self.taille() == 0:
			return "<C: vide>"
		t = "<C: "
		for i in range(self.taille() - 1):
			t += str(self[i]) + " / "
		t += str(self[-1]) + ">"
		return t

	def __add__(self, autre):
		if isinstance(autre, Vect):
			p = self.points + [autre]
			return Chemin(*p)
		elif isinstance(autre, Chemin):
			p = self.points + autre.points
			return Chemin(*p)

	def __mul__(self, a):
		if isinstance(a, int):
			c = Chemin()
			for i in range(a):
				c += self
		elif isinstance(a, float):
			assert 0 <= a and a <= 1
			c = Chemin(self[0])
			lmax = self.longueur()
			l = 0
			for i in range(self.taille() - 1):
				n = (self[i+1] - self[i]).norme()
				if l + n > lmax:
					fin = self[i] + (self[i+1] - self[i]) * (self[i+1] - self[i]).norme() * (lmax - l)
					c += fin
					return c
				else:
					l += (self[i+1] - self[i]).norme()
					c += self[i+1]
			return c

	def taille(self):
		return len(self.points)
	
	def longueur(self):
		if len(self.points) <= 1:
			return 0
		segments = [self[i+1] - self[i] for i in range(self.taille()-1)]
		l = sum([s.norme() for s in segments])
		return l

	def proximite(self, c, r):
		#print("calcul proximite")
		#print("chemin:", self, self.longueur())
		assert self.taille() >=1, "Chemin vide"
		l = 0
		for i in range(self.taille() - 1):
			p = Vect.proximite(self[i], self[i+1], c, r)
			if p != None:
				l += (p - self[i]).norme()
				#print("trouvé !", l)
				#print("p", p)
				#print("c", c)
				#print("last seg len", (p - self[i]).norme())
				return (l, True)
			else:
				l += (self[i+1] - self[i]).norme()
				#print("pas trouvé", l)
		return (l, False)

	def reversed(self):
		return Chemin(*[p for p in reversed(self.points)])

	def reverse(self):
		self.points.reverse()

	def adapt_from(self, point):
		if (self[0] - point).norme() > (self[-1] - point).norme():
			return self.reversed()
		return self

class Triangle:
	def __init__(self, *points):
		self.points = list(points)
		assert len(self.points) == 3, "ERREUR : Mauvais nombre de sommets"

	def __getitem__(self, i):
		return self.points[i%3]

	def __eq__(self, t2):
		return isinstance(t2, Triangle) and all([p in t2 for p in t2.points])

	def __contains__(self, v):
		return isinstance(v, Vect) and True in [p == v for p in self.points]

	def __str__(self):
		return "<T: " + str(self[0]) + " / " + str(self[1]) + " / " + str(self[2]) + ">"

	def dans(self, P):
		return Vect.get_lz(self[0], self[1], P) > 0 and Vect.get_lz(self[1], self[2], P) > 0 and Vect.get_lz(self[2], self[0], P) > 0

	def plat(self):
		for i in range(3):
			if self[i] == self[i+1]:
				return (self[i-1], self[i])

	def centre_g(self):
		p = self.plat()
		if p != None:
			return (p[0] + p[1]) * (1/2)
		d = 2*self[0][0]*self[1][1] - 2*self[0][0]*self[2][1] -2*self[1][0]*self[0][1] +2*self[1][0]*self[2][1] +2*self[2][0]*self[0][1] -2*self[2][0]*self[1][1]
		x = -((self[2][1] - self[1][1]) * self[0][0]**2 + (self[0][1] - self[2][1]) * self[1][0]**2 + (self[1][1] - self[0][1]) * self[2][0]**2 + (self[2][1] - self[1][1]) * self[0][1]**2 + (self[0][1] - self[2][1]) * self[1][1]**2 + (self[1][1] - self[0][1]) * self[2][1]**2)/d
		y = -((self[1][0] - self[2][0]) * self[0][0]**2 + (self[2][0] - self[0][0]) * self[1][0]**2 + (self[0][0] - self[1][0]) * self[2][0]**2 + (self[1][0] - self[2][0]) * self[0][1]**2 + (self[2][0] - self[0][0]) * self[1][1]**2 + (self[0][0] - self[1][0]) * self[2][1]**2)/d
		return Vect(x, y)

	def rayon_circonscrit(self):
		rayon = self.centre_g() - self[0]
		return rayon.norme()

	def angle_droit(self, i):
		return (self[i] - self[i+1]) * (self[i] - self[i-1]) <= 10**-5

	def rectangle(self):
		for i in range(3):
			if self.angle_droit(i):
				return i
		return None

	def projeter_h(self, i):
		a = self[i+1]
		b = self[i-1]
		c = self[i]
		#print("a", a, "b", b, "c", c)
		ac = self[i] - self[i+1]
		ab = self[i-1] - self[i+1]
		u = ab.normalised()
		#print("u", u)
		#print("ac", ac, "ab", ab)
		ah = ac * ab / ab.norme()
		#print("ah", ah)
		h = self[i+1] + u * ah
		return h

	def edge(self):
		n = 0
		maxi = 0
		for i in range(3):
			l = (self[i+1] - self[i+2]).norme()
			if l > maxi:
				maxi = l
				n = i
		return n

	def polya_rec(self, d):
		if d >= self.rayon_circonscrit():
			return Chemin(self.centre_g())
		i = self.edge()
		h = self.projeter_h(i)
		return Triangle(self[i], self[i+1], h).polya(d) + Triangle(self[i], h, self[i-1]).polya(d)

	def polya(self, d, c=None):
		pile = [self]
		terminaux = []
		while pile != []:
			t = pile.pop()
			if d >= t.rayon_circonscrit():
				terminaux.append(t)
			else:
				i = t.edge()
				h = t.projeter_h(i)
				t1, t2 = Triangle(t[i], t[i+1], h), Triangle(t[i], h, t[i-1])
				pile.append(t2)
				pile.append(t1)

		if c != None:
			for t in terminaux:
				if t.dans(c):
					#print("dedans", t, c)
					r = (t.centre_g() - c).norme()
					#print(t.centre_g())
					#print("r", r, "d", d)

		chemin = Chemin(*[t.centre_g() for t in terminaux])
		if c != None:
			for p in chemin.points:
				if (p-c).norme() <= d:
					pass
					#print("toujours dans chemin")
		return chemin

	
	def area(self):
		i = self.edge()
		base = (self[i+1] - self[i-1]).norme()
		h = self.projeter_h(i)
		hauteur = self[i] - h
		return base * hauteur.norme() / 2

class Polygon:
	def __init__(self, *points):
		self.points = list(points)
		self.triangles = Polygon.triangulation(self)

	def __str__(self):
		t = "<P: "
		for i in range(self.taille() - 1):
			t += str(self[i]) + " / "
		t += str(self[self.taille() - 1]) + ">"
		return t

	def __getitem__(self, i):
		if isinstance(i, int):
			return self.points[i % self.taille()]
		if i.step == None:
			return [self[j] for j in range(i.start, i.stop)]
		else:
			return [self[j] for j in range(i.start, i.step, i.stop)]
			
	def taille(self):
		return len(self.points)
	
	def rand_point(self):
		x_min = min([self[i][0] for i in range(self.taille())])
		x_max = max([self[i][0] for i in range(self.taille())])
		y_min = min([self[i][1] for i in range(self.taille())])
		y_max = max([self[i][1] for i in range(self.taille())])
		p = Vect(x_min + random.random() * (x_max - x_min), y_min + random.random() * (y_max + y_min))
		while not self.dans(p):
			p = Vect(x_min + random.random() * (x_max - x_min), y_min + random.random() * (y_max + y_min))
		return p

	def plus_basse_abscisse(self):
		m = self[0][0]
		j = 0
		for i in range(1, self.taille()):
			if m > self[i][0]:
				j = i
		return j

	def sample(self, start, end):
		points = []
		if start > end:
			for i in range(start, self.taille()):
				points.append(self[i].copy())
			for i in range(0, end+1):
				points.append(self[i].copy())
		else:
			for i in range(start, end+1):
				points.append(self[i].copy())
		p = Polygon(*points)
		return p

	def dans(self, point):
		return True in [triangle.dans(point) for triangle in self.triangles]

	def segments(self):
		return [(self[i], self[i+1]) for i in range(self.taille())]

	def trianguler(self):
		return Polygon.triangulation(self)

	@classmethod
	def triangulation(cls, polygon):
		if polygon.taille() <= 2:
			return []
		elif polygon.taille() == 3:
			return [Triangle(*polygon.points)]
		else:
			i = polygon.plus_basse_abscisse()
			t = Triangle(*polygon[i-1:i+2])
			mx_lz = -1
			m = 0
			for j in range(polygon.taille()):
				if not polygon[j] in t and t.dans(polygon[j]):
					lz = Vect.get_lz(polygon[i+1], polygon[i-1], polygon[j])
					if lz > mx_lz:
						m = j
						mx_lz = lz
			if mx_lz > 0:
				return Polygon.triangulation(polygon.sample(i, m)) + Polygon.triangulation(polygon.sample(m, i))
			else:
				return [t] + Polygon.triangulation(polygon.sample(i+1, i-1))

	@classmethod
	def random(cls, n, d):
		modules = [random.random()*d/2 for i in range(n)]
		thetas = [random.random()*math.pi*2 for i in range(n)]
		thetas.sort()
		points = [Vect(d/2 + math.cos(thetas[i]) * modules[i], d/2 + math.sin(thetas[i]) * modules[i]) for i in range(n)]
		return Polygon(*points)
	
	def area(self):
		return sum([t.area() for t in self.triangles])

"""polygone = [[0,0],[0.5,-1],[1.5,-0.2],[2,-0.5],[2,0],[1.5,1],[0.3,0,],[0.5,1]]
p = Polygon(Vect(0, 0), Vect(0.5, -1), Vect(1.5, -0.2), Vect(2, -0.5), Vect(2, 0), Vect(1.5, 1), Vect(0.3, 0), Vect(0.5, 1))
triangles = p.trianguler()
print(p)
print("triangles:")
for t in triangles:
	print(t)
"""
class Agent:
	def __init__(self, position):
		self.position = position
		self.distance = 0

	@classmethod
	def get_random_dir(self, longueur):
		angle = random.random() * math.pi * 2
		x = math.cos(angle) * longueur
		y = math.sin(angle) * longueur
		return Vect(x, y)

	def get_path(self, direction, segments):
		pos = self.position.copy()
		u = direction.copy()
		l = direction.norme()
		h, i = Vect.raycast(pos, u, segments, l)
		#print("h i", h, i)
		path = Chemin(pos)
		while h != None:
			#print("chemin", path)
			#print(pos, u, l, h, segments[i][0])
			path += h
			l -= (h - pos).norme()
			r = Vect.reflect(pos, h, segments[i][0])
			pos = h.copy()
			#print("r", r)
			u = (r - h).normalised() * l
			h, i = Vect.raycast(pos, u, segments, l, i)
		pos += u
		path += pos
		return path

	def get_random_path(self, longueur, segments):
		return self.get_path(self.get_random_dir(longueur), segments)

	def marche(self, chemin):
		self.position = chemin[-1]
			
class Configuration:

	def __init__(self, n_agents, n_sommets, portee, d, polygone=None, limit=999999):
		if polygone == None:
			self.polygone = Polygon.random(n_sommets, 100)
		else:
			self.polygone = polygone
		#print("pol", self.polygone)
		self.spawn = self.polygone.rand_point()
		self.agents = [Agent(self.spawn) for i in range(n_agents)]
		self.balise = self.polygone.rand_point()
		self.rayon = portee
		self.pas = d
		self.limit = limit
		self.aire = self.polygone.area()
		#print("depart", self.spawn, "balise", self.balise)
		print("aire", self.aire)

	def random_search(self):
		chemins = []
		segments = self.polygone.segments()
		distance = 0
		chemins = [0 for i in self.agents]

		trouve = False
		while not trouve:
			#print("=========\n\n==========")
			m = self.pas
			for i in range(len(self.agents)):
				chemins[i] = self.agents[i].get_random_path(self.pas, segments)
				#print("ag", i, "coords", self.agents[i].position)
				l, t = chemins[i].proximite(self.balise, self.rayon)
				if t:
					trouve = True

				self.agents[i].position = chemins[i][-1]
				m = min(l, m)
			
			distance += m * len(self.agents)
			assert distance <= self.limit, "limite dépassé"
		#print("Trouvé !!")
		return distance

	def curve_search(self):
		#print("depart", self.spawn, "balise", self.balise)
		print("début de recherche")
		curves = [triangle.polya(self.rayon, self.balise) for triangle in self.polygone.triangles]
		print("polya finie")
		chemins = [Chemin(self.agents[i].position) for i in range(len(self.agents))]
		while len(curves) > 0:
			i_min = 0
			for i in range(len(self.agents)):
				if chemins[i].longueur() <= chemins[i_min].longueur():
					i_min = i

			j_min = 0
			for j in range(len(curves)):
				c_min = curves[j_min].adapt_from(chemins[i_min][-1]) + chemins[i_min]
				chemin_courant = curves[j].adapt_from(chemins[i_min][-1]) + chemins[i_min]
				if c_min.longueur() > chemin_courant.longueur():
					j_min = j
			chemins[i_min] = chemins[i_min] + c_min
			curves.pop(j_min)

		print("répartition finie")
		#print("chemins")

		l_min = None
		for i in range(len(chemins)):
			#print(chemins[i])
			l, trouve = chemins[i].proximite(self.balise, self.rayon)
			#print("l, trouve", l, trouve)
			if trouve:
				if l_min == None:
					l_min = l
				else:
					l_min = min(l, l_min)

		distance = sum([min(c.longueur(), l_min) for c in chemins if c.taille() >= 1])
		return distance

	def gaz_search(self, a):
		segments = self.polygone.segments()
		directions = [self.agents[i].get_random_dir(self.pas) for i in range(len(self.agents))]
		chemins = [None for i in range(len(self.agents))]
		trouve = True in [ (self.agents[i].position - self.balise).norme() <= self.rayon for i in range(len(self.agents))]
		distance = 0

		while not trouve:
			for i in range(len(self.agents)):
				d = Vect(0, 0)
				for j in range(len(self.agents)):
					u = self.agents[i].position - self.agents[j].position
					if i != j and u.norme() > 10**-5:
						d += u.normalised() * (1 / (u.norme()**5))
				directions[i] += d * a
				directions[i].normalise()
				directions[i] = directions[i] * self.pas

			l_min = self.pas
			print("===\n\n===")
			for i in range(len(self.agents)):
				print(self.agents[i].position)
				#print("dir", directions[i].norme())
				chemins[i] = self.agents[i].get_path(directions[i], segments)
				#print("che", chemins[i].longueur())
				l, t = chemins[i].proximite(self.balise, self.rayon)
				#print("→", l, t)
				if t:
					trouve = t
					l_min = min(l, l_min)
				self.agents[i].marche(chemins[i])
			distance += len(self.agents) * l_min

		return distance



"""A = Vect(1, 2)
B = Vect(16, 2)
C = Vect(6, 14)
T = Triangle(A, B, C)
c = Polygon(Vect(0, 0), Vect(1, 0), Vect(1, 1), Vect(0, 1))
print(T.centre_g())
print(T.rayon_circonscrit())
chemin = T.polya(4)
print(T.area())
print(chemin.longueur())"""

pol = Polygon(Vect(0, 0), Vect(100, 0), Vect(100, 100), Vect(0, 100))
r = 4
n = 2
pas = 1
conf = Configuration(n, 10, r, pas, pol)
conf.spawn = Vect(61.098, 41.587)
conf.balise = Vect(67.9, 89.95)

def random_simulation(k, n_agents, n_sommets, bornes, pas, limit=999999, pol=None):
	X = []
	Y = []
	for i in range(k):
		print("===n°", i)
		r = random.uniform(*bornes)
		conf = Configuration(n_agents, n_sommets, r, pas, pol, limit)
		try:
			l = conf.random_search()
		except:
			pass
		else:
			print("l", l)
			X.append(conf.aire/(r**2))
			Y.append(l)
	return X, Y

def curve_simulation(k, n_agents, n_sommets, bornes, pas, pol=None):
	X, Y = [], []
	for i in range(k):
		print("===n°", i)
		r = random.uniform(*bornes)
		conf = Configuration(n_agents, n_sommets, r, pas, pol)
		try:
			l = conf.curve_search()
		except:
			pass
		else:
			print(l)
			#X.append(conf.aire/(r**2))
			X.append(conf.aire/(r**2))
			Y.append(l)
	return X, Y

def save(title, X, Y):
	d = {"x":X, "y":Y}
	f = open(title, "wt")
	f.write(json.dumps(d))
	f.close()

#X, Y = random_simulation(2000, 10, 10, (2, 5), 5, 200000, pol)
X, Y = curve_simulation(2000, 2, 5, (2, 5), 5, pol)
print("simulation finished")
save("test2", X, Y)
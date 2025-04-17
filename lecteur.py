import matplotlib.pyplot as plt, json

def replace(l, a, b):
	l2 = []
	for e in l:
		if e == a:
			l2.append(b)
		else:
			l2.append(e)
	return l2

def load(title):
	f = open(title, "rt")
	d = json.loads(f.read())
	f.close()
	X = d["x"]
	Y = d["y"]
	return X, Y

def ecartT(X):
	if len(X) in [0, 1]:
		return 0
	moy = sum(X)/len(X)
	e = (sum([(x - moy)**2 for x in X])/(len(X)-1))**0.5
	return e

X, Y = load("test2")
"""plt.plot(X, Y, '+')
plt.xlabel("rapport des surfaces")
plt.ylabel("distance parcourue")
plt.show()"""

bars = {}
x_max = max(X)
pas = 300
tranche = lambda y, pas : int((y - (y % pas))/pas)
b_max = tranche(x_max, pas)
for i in range(len(X)):
	t = tranche(X[i], pas)
	if not t in bars:
		bars[t] = []
	bars[t].append(Y[i])

moyennes = []

for k in range(b_max+1):
	if not k in bars or len(bars[k]) == 0:
		moyennes.append(0)
	else:
		s = sum(bars[k])
		moyennes.append(s/len(bars[k]))

names = [f"{i*pas} - {(i+1)*pas}" for i in range(b_max + 1)]

eT = []

for k in range(b_max+1):
	if not k in bars:
		eT.append(0)
	else:
		eT.append(ecartT(bars[k]))

plt.bar(names, moyennes, color = "#A0AAE4", edgecolor="red", linewidth=3, yerr=eT, ecolor = "green",capsize = 10)
plt.xlabel("rapport des surfaces")
plt.ylabel("distance parcourue")
plt.show()
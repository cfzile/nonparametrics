# <=
class GaltonStatistic:

    def __init__(self, n):
        self.n = n

    def getP(self, g):
        if g >= 0 and g <= self.n:
            return 1. / (self.n + 1)
        return 0

    def getF(self, g):
        if g < 0:
            return 0
        if g > self.n:
            return 1
        return g / (self.n + 1.)


# <= >=
class GaltonTestTS:
    def findC1(self):
        c = self.n
        while abs(1 - self.galton_statistic.getF(c) - self.alpha / 2.) > eps:
            c -= 1
        return c

    def __init__(self, x, y, alpha=0.05):
        self.n = len(x)
        x = sorted(x)
        y = sorted(y)
        self.galton_statistic = GaltonStatistic(self.n)
        self.g = sum([x[i] > y[i] for i in range(self.n)])
        self.alpha = alpha

    def p_value(self):
        return 2 * min(1 - self.galton_statistic.getF(self.g), self.galton_statistic.getF(self.g))

    def isRejected(self):
        self.c1 = self.findC1()
        self.c2 = self.n - self.c1
        return self.g >= self.c1 and self.g <= self.c2


# <=
class GaltonTestG:
    def findC1(self):
        c = self.n
        while abs(1 - self.galton_statistic.getF(c) - self.alpha) > eps and c >= 0:
            c -= 1.
        return c

    def __init__(self, x, y, alpha=0.05):
        self.n = len(x)
        x = sorted(x)
        y = sorted(y)
        self.galton_statistic = GaltonStatistic(self.n)
        self.g = sum([x[i] > y[i] for i in range(self.n)])
        self.alpha = alpha

    def p_value(self):
        return 1 - self.galton_statistic.getF(self.g)

    def isRejected(self):
        self.c1 = self.findC1()
        return self.g >= self.c1


# >=
class GaltonTestL:
    def findC1(self):
        c = self.n
        while abs(self.galton_statistic.getF(c) - self.alpha) > eps:
            c -= 1
        return c

    def __init__(self, x, y, alpha=0.05):
        self.n = len(x)
        x = sorted(x)
        y = sorted(y)
        self.galton_statistic = GaltonStatistic(self.n)
        self.g = sum([x[i] > y[i] for i in range(self.n)])
        self.alpha = alpha

    def p_value(self):
        return self.galton_statistic.getF(self.g + 1.)

    def isRejected(self):
        self.c1 = self.findC1()
        return self.g >= self.c1


class GaltonTest:
    def __init__(self, x, y, alternative='two-sided', alpha=0.05):
        if alternative == 'two-sided':
            self.test = GaltonTestTS(x, y, alpha)
        elif alternative == 'less':
            self.test = GaltonTestL(x, y, alpha)
        else:
            self.test = GaltonTestG(x, y, alpha)

    def p_value(self):
        return self.test.p_value()

    def isRejected(self):
        return self.test.p_value() <= self.test.alpha
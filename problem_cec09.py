import numpy as np
import abc

class Problem:
    def __init__(self, n, m):
        self._n = n
        self._m = m

    @abc.abstractmethod
    def is_feasible(self, x):
        """
        :return: True if the point x is feasible, False othwerwise
        """

    @abc.abstractmethod
    def evaluate(self, x):
        """
        :return: the objective functions evaluated in x
        """
    @property
    def n(self):
        return self._n

    @property
    def m(self):
        return self._m

    def fconstr_a(self,x):
        J = np.arange(len(x)-2)
        return  (3-2*x[J+1])*x[J+1] - x[J] - 2*x[J+2] + 1

    def fconstr_b(self,x):
        J = np.arange(len(x)-2)
        return  (3-2*x[J+1])*x[J+1] - x[J] - 2*x[J+2] + 2.5

    def fconstr_c(self,x):
        J = np.arange(len(x)-1)
        return x[J]**2 + x[J+1]**2 + x[J]*x[J+1] - 2*x[J] - 2*x[J+1] + 1

    def fconstr_d(self,x):
        J = np.arange(len(x)-1)
        return x[J]**2 + x[J+1]**2 + x[J]*x[J+1] - 1

    def fconstr_e(self,x):
        J = np.arange(len(x)-2)
        return (3-0.5*x[J+1])*x[J+1] - x[J] -2*x[J+2] +1

    def fconstr_f(self,x):
        J = np.arange(len(x)-2)
        return np.array([np.sum((3-0.5*x[J+1])*x[J+1] - x[J] -2*x[J+2] +1)])

    def fconstr_z(self,x):
        return np.array([-1.0])

class Prob1(Problem):
    @staticmethod
    def get_matlab_name():
        return "UF1"

    def __str__(self):
        return "Cec09_Prob1"

    def __repr__(self):
        return self.__str__()

    def __init__(self, n, ctype='z', constrained=False):
        assert n >= 3
        """
        :param n: the number of optimization variables
        :return:
        """

        super(Prob1, self).__init__(n, 2)
        self.constrained = constrained

        self.lb = -1 * np.ones(n)
        self.lb[0] = 0.0

        self.ub = np.ones(n)
        self.ctype = ctype;

    def eval_constr(self,x):
        if self.ctype == 'a':
            return self.fconstr_a(x)
        if self.ctype == 'b':
            return self.fconstr_b(x)
        if self.ctype == 'c':
            return self.fconstr_c(x)
        if self.ctype == 'd':
            return self.fconstr_d(x)
        if self.ctype == 'e':
            return self.fconstr_e(x)
        if self.ctype == 'f':
            return self.fconstr_f(x)
        if self.ctype == 'z':
            return self.fconstr_z(x)

    def is_feasible(self, x):
        if self.constrained:
            return (self.lb <= x).all() and (self.ub >= x).all()
        return True

    def evaluate1(self, x):
        n = self._n
        #print(self.eval_constr(x))
        #input()
        J1 = np.arange(2, n, 2)
        J2 = np.arange(1, n, 2)
        valpen = 0.1*np.maximum(0.,np.max(self.eval_constr(x)))
        obj_1 = x[0] + (2 / len(J1)) * sum(
            [(x[i] - np.sin(6 * np.pi * x[0] + (i + 1) * np.pi / n)) ** 2 for i in J1])
        obj_2 = 1 - np.sqrt(x[0]) + (2 / len(J2)) * sum(
            [(x[i] - np.sin(6 * np.pi * x[0] + (i + 1) * np.pi / n)) ** 2 for i in J2])

        return np.array([obj_1+valpen, obj_2+valpen])

    def evaluate(self, x):
        n = self._n
        #print(self.eval_constr(x))
        #input()
        J1 = np.arange(2, n, 2)
        J2 = np.arange(1, n, 2)
        obj_1 = x[0] + (2 / len(J1)) * sum(
            [(x[i] - np.sin(6 * np.pi * x[0] + (i + 1) * np.pi / n)) ** 2 for i in J1])
        obj_2 = 1 - np.sqrt(x[0]) + (2 / len(J2)) * sum(
            [(x[i] - np.sin(6 * np.pi * x[0] + (i + 1) * np.pi / n)) ** 2 for i in J2])

        return np.array([obj_1, obj_2])

class Prob2(Problem):
    @staticmethod
    def get_matlab_name():
        return "UF2"

    def __str__(self):
        return "Cec09_Prob2"

    def __repr__(self):
        return self.__str__()

    def __init__(self, n, ctype='z', constrained=False):
        assert n >= 3
        """
        :param n: the number of optimization variables
        :return:
        """

        super(Prob2, self).__init__(n, 2)
        self.constrained = constrained

        self.lb = -1 * np.ones(n)
        self.lb[0] = 0.0
        self.ub = np.ones(n)
        self.ctype = ctype;

    def eval_constr(self,x):
        if self.ctype == 'a':
            return self.fconstr_a(x)
        if self.ctype == 'b':
            return self.fconstr_b(x)
        if self.ctype == 'c':
            return self.fconstr_c(x)
        if self.ctype == 'd':
            return self.fconstr_d(x)
        if self.ctype == 'e':
            return self.fconstr_e(x)
        if self.ctype == 'f':
            return self.fconstr_f(x)
        if self.ctype == 'z':
            return self.fconstr_z(x)

    def is_feasible(self, x):
        if self.constrained:
            return (self.lb <= x).all() and (self.ub >= x).all()
        return True

    def evaluate(self, x):
        n = self._n
        J1 = np.arange(2, n, 2)
        J2 = np.arange(1, n, 2)

        y_odd = 2 * sum([(x[j] -
                                    (0.3 * x[0] ** 2 * np.cos(
                                        24 * np.pi * x[0] + 4 * (j + 1) * np.pi / n) + 0.6 * x[0])
                                    * np.cos(6 * np.pi * x[0] + (j + 1) * np.pi / n)) ** 2 for j in J1]) / len(J1)

        y_even = 2 * sum([(x[j] -
                                     (0.3 * x[0] ** 2 * np.cos(
                                         24 * np.pi * x[0] + 4 * (j + 1) * np.pi / n) + 0.6 * x[0])
                                     * np.sin(6 * np.pi * x[0] + (j + 1) * np.pi / n)) ** 2 for j in J2]) / len(J2)

        obj_1 = x[0] + y_odd
        obj_2 = 1 - np.sqrt(x[0]) + y_even
        return np.array([obj_1, obj_2])

class Prob3(Problem):
    @staticmethod
    def get_matlab_name():
        return "UF3"

    def __str__(self):
        return "Cec09_Prob3"

    def __repr__(self):
        return self.__str__()

    def __init__(self, n, ctype='z', constrained=False):
        assert n >= 3
        """
        :param n: the number of optimization variables
        :return:
        """
        super(Prob3, self).__init__(n, 2)

        self.constrained = constrained

        self.lb = np.zeros(n)
        self.ub = np.ones(n)
        self.ctype = ctype;

    def eval_constr(self,x):
        if self.ctype == 'a':
            return self.fconstr_a(x)
        if self.ctype == 'b':
            return self.fconstr_b(x)
        if self.ctype == 'c':
            return self.fconstr_c(x)
        if self.ctype == 'd':
            return self.fconstr_d(x)
        if self.ctype == 'e':
            return self.fconstr_e(x)
        if self.ctype == 'f':
            return self.fconstr_f(x)
        if self.ctype == 'z':
            return self.fconstr_z(x)

    def is_feasible(self, x):
        if self.constrained:
            return (self.lb <= x).all() and (self.ub >= x).all()
        return True

    def evaluate(self, x):
        n = self._n
        J1 = np.arange(2, n, 2)
        J2 = np.arange(1, n, 2)

        y = [x[j] - x[0] ** (0.5 * (1 + 3 * (j - 1) / (n - 2))) for j in np.arange(n)]
        y_1 = sum([y[j] ** 2 for j in J1])

        y_2 = np.prod([np.cos(20 * y[j] * np.pi / np.sqrt(j + 1.0)) for j in J1])

        y_odd = 2 * (4 * y_1 - 2 * y_2 + 2) / len(J1)

        obj_1 = x[0] + y_odd

        y_even = 2 * (
            4 * sum([y[j] ** 2 for j in J2])
            - 2 * np.prod([np.cos(20 * y[j] * np.pi / np.sqrt(j + 1.0)) for j in J2])
            + 2) / len(J2)

        obj_2 = 1 - np.sqrt(x[0]) + y_even
        return np.array([obj_1, obj_2])

class Prob4(Problem):
    @staticmethod
    def get_matlab_name():
        return "UF4"

    def __str__(self):
        return "Cec09_Prob4"

    def __repr__(self):
        return self.__str__()

    def __init__(self, n, ctype='z', constrained=False):
        assert n >= 3
        """
        :param n: the number of optimization variables
        :return:
        """
        super(Prob4, self).__init__(n, 2)

        self.constrained = constrained

        self.lb = -2 * np.ones(n)
        self.ub = 2 * np.ones(n)

        self.lb[0] = 0.0
        self.ub[0] = 1.0
        self.ctype = ctype;

    def eval_constr(self,x):
        if self.ctype == 'a':
            return self.fconstr_a(x)
        if self.ctype == 'b':
            return self.fconstr_b(x)
        if self.ctype == 'c':
            return self.fconstr_c(x)
        if self.ctype == 'd':
            return self.fconstr_d(x)
        if self.ctype == 'e':
            return self.fconstr_e(x)
        if self.ctype == 'f':
            return self.fconstr_f(x)
        if self.ctype == 'z':
            return self.fconstr_z(x)

    def is_feasible(self, x):
        if self.constrained:
            return (self.lb <= x).all() and (self.ub >= x).all()
        return True

    def evaluate(self, x):
        n = self._n
        J1 = np.arange(2, n, 2)
        J2 = np.arange(1, n, 2)

        y = [x[j] - np.sin(6 * np.pi * x[0] + (j + 1) * np.pi / n) for j in np.arange(n)]

        obj_1 = x[0] + 2 * sum([np.abs(y[j] / (1 + np.exp(2 * np.abs(y[j])))) for j in J1]) / len(J1)
        obj_2 = 1 - x[0] ** 2 + 2 * sum(
            [np.abs(y[j] / (1 + np.exp(2 * np.abs(y[j])))) for j in J2]) / len(J2)
        return np.array([obj_1, obj_2])

class Prob5(Problem):
    @staticmethod
    def get_matlab_name():
        return "UF5"

    def __str__(self):
        return "Cec09_Prob5"

    def __repr__(self):
        return self.__str__()

    def __init__(self, n, ctype='z', constrained=False):
        assert n >= 3
        """
        :param n: the number of optimization variables
        :return:
        """
        super(Prob5, self).__init__(n, 2)

        self.constrained = constrained

        self.lb = -1 * np.ones(n)
        self.ub = np.ones(n)
        self.lb[0] = 0.0
        self.ctype = ctype;

    def eval_constr(self,x):
        if self.ctype == 'a':
            return self.fconstr_a(x)
        if self.ctype == 'b':
            return self.fconstr_b(x)
        if self.ctype == 'c':
            return self.fconstr_c(x)
        if self.ctype == 'd':
            return self.fconstr_d(x)
        if self.ctype == 'e':
            return self.fconstr_e(x)
        if self.ctype == 'f':
            return self.fconstr_f(x)
        if self.ctype == 'z':
            return self.fconstr_z(x)

    def is_feasible(self, x):
        if self.constrained:
            return (self.lb <= x).all() and (self.ub >= x).all()

        return True

    def evaluate(self, x):
        n = self._n
        J1 = np.arange(2, n, 2)
        J2 = np.arange(1, n, 2)

        N = 10
        eps = 0.1

        def h(t):
            return 2 * t ** 2 - np.cos(4 * np.pi * t) + 1

        y = [x[j] - np.sin(6 * np.pi * x[0] + (j + 1) * np.pi / n) for j in np.arange(n)]

        obj_1 = x[0] + (1 / (2 * N) + eps) * np.abs(np.sin(2 * N * np.pi * x[0])) + 2 * sum(
            [h(y[j]) for j in J1]) / len(J1)
        obj_2 = 1 - x[0] + (1 / (2 * N) + eps) * np.abs(
            np.sin(2 * N * np.pi * x[0])) + 2 * sum([h(y[j]) for j in J2]) / len(J2)
        return np.array([obj_1, obj_2])

class Prob6(Problem):
    @staticmethod
    def get_matlab_name():
        return "UF6"

    def __str__(self):
        return "Cec09_Prob6"

    def __repr__(self):
        return self.__str__()

    def __init__(self, n, ctype='z', constrained=False):
        assert n >= 3
        """
        :param n: the number of optimization variables
        :return:
        """
        super(Prob6, self).__init__(n, 2)

        self.constrained = constrained

        self.lb = -1 * np.ones(n)
        self.ub = np.ones(n)
        self.lb[0] = 0.0
        self.ctype = ctype;

    def eval_constr(self,x):
        if self.ctype == 'a':
            return self.fconstr_a(x)
        if self.ctype == 'b':
            return self.fconstr_b(x)
        if self.ctype == 'c':
            return self.fconstr_c(x)
        if self.ctype == 'd':
            return self.fconstr_d(x)
        if self.ctype == 'e':
            return self.fconstr_e(x)
        if self.ctype == 'f':
            return self.fconstr_f(x)
        if self.ctype == 'z':
            return self.fconstr_z(x)

    def is_feasible(self, x):
        if self.constrained:
            return (self.lb <= x).all() and (self.ub >= x).all()
        return True

    def evaluate(self, x):
        n = self._n
        J1 = np.arange(2, n, 2)
        J2 = np.arange(1, n, 2)

        N = 2
        eps = 0.1
        y = [x[j] - np.sin(6 * np.pi * x[0] + (j + 1) * np.pi / self._n) for j in np.arange(n)]

        obj_1 = x[0] + (
            np.maximum(0.0, 2 * (1. / (2 * N) + eps) * np.sin(2 * N * np.pi * x[0]))
            + 2 * (
                4 * sum([y[j] ** 2 for j in J1]) - 2 * np.prod(
                    [np.cos(20 * y[j] * np.pi / np.sqrt(j + 1)) for j in J1])
                + 2) / len(J1))

        obj_2 = 1 - x[0] + (
            np.maximum(0.0, 2 * (1. / (2 * N) + eps) * np.sin(2 * N * np.pi * x[0]))
            + 2 * (
                4 * sum([y[j] ** 2 for j in J2]) - 2 * np.prod(
                    [np.cos(20 * y[j] * np.pi / np.sqrt(j + 1)) for j in J2])
                + 2) / len(J2))

        return np.array([obj_1, obj_2])

class Prob7(Problem):
    @staticmethod
    def get_matlab_name():
        return "UF7"

    def __str__(self):
        return "Cec09_Prob7"

    def __repr__(self):
        return self.__str__()

    def __init__(self, n, ctype='z', constrained=False):
        assert n >= 3
        """
        :param n: the number of optimization variables
        :return:
        """

        super(Prob7, self).__init__(n, 2)

        self.constrained = constrained

        self.lb = -1 * np.ones(n)
        self.ub = np.ones(n)
        self.lb[0] = 0.0
        self.ctype = ctype;

    def eval_constr(self,x):
        if self.ctype == 'a':
            return self.fconstr_a(x)
        if self.ctype == 'b':
            return self.fconstr_b(x)
        if self.ctype == 'c':
            return self.fconstr_c(x)
        if self.ctype == 'd':
            return self.fconstr_d(x)
        if self.ctype == 'e':
            return self.fconstr_e(x)
        if self.ctype == 'f':
            return self.fconstr_f(x)
        if self.ctype == 'z':
            return self.fconstr_z(x)

    def is_feasible(self, x):
        if self.constrained:
            return (self.lb <= x).all() and (self.ub >= x).all()
        return True

    def evaluate(self, x):
        n = self._n
        J1 = np.arange(2, n, 2)
        J2 = np.arange(1, n, 2)

        y = [x[j] - np.sin(6 * np.pi * x[0] + (j + 1) * np.pi / self._n) for j in np.arange(n)]
        obj_1 = x[0] ** (1 / 5) + (2 / len(J1)) * sum([y[j] ** 2 for j in J1])
        obj_2 = 1 - x[0] ** (1 / 5) + (2 / len(J2)) * sum([y[j] ** 2 for j in J2])
        return np.array([obj_1, obj_2])

class Prob8(Problem):
    @staticmethod
    def get_matlab_name():
        return "UF8"

    def __str__(self):
        return "Cec09_Prob8"

    def __repr__(self):
        return self.__str__()

    def __init__(self, n, ctype='z', constrained=False):
        assert n >= 5
        """
        :param n: the number of optimization variables
        :return:
        """

        super(Prob8, self).__init__(n, 3)

        self.constrained = constrained

        self.lb = -2 * np.ones(n)
        self.ub = 2 * np.ones(n)
        self.lb[0] = 0.0
        self.lb[1] = 0.0
        self.ub[0] = 1.0
        self.ub[1] = 1.0
        self.ctype = ctype;

    def eval_constr(self,x):
        if self.ctype == 'a':
            return self.fconstr_a(x)
        if self.ctype == 'b':
            return self.fconstr_b(x)
        if self.ctype == 'c':
            return self.fconstr_c(x)
        if self.ctype == 'd':
            return self.fconstr_d(x)
        if self.ctype == 'e':
            return self.fconstr_e(x)
        if self.ctype == 'f':
            return self.fconstr_f(x)
        if self.ctype == 'z':
            return self.fconstr_z(x)

    def is_feasible(self, x):
        if self.constrained:
            return (self.lb <= x).all() and (self.ub >= x).all()
        return True

    def evaluate(self, x):
        n = self._n
        J1 = np.arange(3, n, 3)
        J2 = np.arange(4, n, 3)
        J3 = np.arange(2, n, 3)

        obj_1 = np.cos(0.5 * np.pi * x[0]) * np.cos(0.5 * np.pi * x[1]) + (2 / len(J1)) * sum(
            [(x[j] - 2 * x[1] * np.sin(2 * np.pi * x[0] + (j + 1) * np.pi / n)) ** 2 for j in J1]
        )

        obj_2 = np.cos(0.5 * np.pi * x[0]) * np.sin(0.5 * np.pi * x[1]) + (2 / len(J2)) * sum(
            [(x[j] - 2 * x[1] * np.sin(2 * np.pi * x[0] + (j + 1) * np.pi / self._n)) ** 2 for j in J2]
        )
        obj_3 = np.sin(0.5 * np.pi * x[0]) + (2 / len(J3)) * sum(
            [(x[j] - 2 * x[1] * np.sin(2 * np.pi * x[0] + (j + 1) * np.pi / self._n)) ** 2 for j in J3]
        )
        return np.array([obj_1, obj_2, obj_3])

class Prob9(Problem):
    @staticmethod
    def get_matlab_name():
        return "UF9"

    def __str__(self):
        return "Cec09_Prob9"

    def __repr__(self):
        return self.__str__()

    def __init__(self, n, ctype='z', constrained=False):
        assert n >= 5
        """
        :param n: the number of optimization variables
        :return:
        """

        super(Prob9, self).__init__(n, 3)

        self.constrained = constrained

        self.lb = -2 * np.ones(n)
        self.ub = 2 * np.ones(n)
        self.lb[0] = 0.0
        self.lb[1] = 0.0
        self.ub[0] = 1.0
        self.ub[1] = 1.0
        self.ctype = ctype;

    def eval_constr(self,x):
        if self.ctype == 'a':
            return self.fconstr_a(x)
        if self.ctype == 'b':
            return self.fconstr_b(x)
        if self.ctype == 'c':
            return self.fconstr_c(x)
        if self.ctype == 'd':
            return self.fconstr_d(x)
        if self.ctype == 'e':
            return self.fconstr_e(x)
        if self.ctype == 'f':
            return self.fconstr_f(x)
        if self.ctype == 'z':
            return self.fconstr_z(x)

    def is_feasible(self, x):

        if self.constrained:
            return (self.lb <= x).all() and (self.ub >= x).all()

        return (-1.0e+3 <= x).all() and (1.0e+3 >= x).all()

    def evaluate(self, x):
        n = self._n
        J1 = np.arange(3, n, 3)
        J2 = np.arange(4, n, 3)
        J3 = np.arange(2, n, 3)
        eps = 0.1

        obj_1 = 0.5 * (np.maximum(0.0, (1 + eps) * (1 - 4 * (2 * x[0] - 1) ** 2))
                              + 2 * x[0]) * x[1] + (2 / len(J1)) * sum([(x[j] - 2 * x[1] * np.sin(
            2 * np.pi * x[0] + (j + 1) * np.pi / n)) ** 2 for j in J1])

        obj_2 = 0.5 * (np.maximum(0.0, (1 + eps) * (1 - 4 * (2 * x[0] - 1) ** 2))
                              - 2 * x[0] + 2) * x[1] + (2 / len(J2)) * sum([(x[j] - 2 * x[1] * np.sin(
            2 * np.pi * x[0] + (j + 1) * np.pi / n)) ** 2 for j in J2])

        obj_3 = 1 - x[1] + (2 / len(J3)) * sum([(x[j] - 2 * x[1] * np.sin(
            2 * np.pi * x[0] + (j + 1) * np.pi / n)) ** 2 for j in J3])

        return np.array([obj_1, obj_2, obj_3])

class Prob10(Problem):
    @staticmethod
    def get_matlab_name():
        return "UF10"

    def __str__(self):
        return "Cec09_Prob10"

    def __repr__(self):
        return self.__str__()

    def __init__(self, n, ctype='z', constrained=False):
        assert n >= 5
        """
        :param n: the number of optimization variables
        :return:
        """
        super(Prob10, self).__init__(n, 3)

        self.constrained = constrained

        self.lb = -2 * np.ones(n)
        self.ub = 2 * np.ones(n)
        self.lb[0] = 0.0
        self.lb[1] = 0.0
        self.ub[0] = 1.0
        self.ub[1] = 1.0
        self.ctype = ctype;

    def eval_constr(self,x):
        if self.ctype == 'a':
            return self.fconstr_a(x)
        if self.ctype == 'b':
            return self.fconstr_b(x)
        if self.ctype == 'c':
            return self.fconstr_c(x)
        if self.ctype == 'd':
            return self.fconstr_d(x)
        if self.ctype == 'e':
            return self.fconstr_e(x)
        if self.ctype == 'f':
            return self.fconstr_f(x)
        if self.ctype == 'z':
            return self.fconstr_z(x)

    def is_feasible(self, x):
        if self.constrained:
            return (self.lb <= x).all() and (self.ub >= x).all()
        return True

    def evaluate(self, x):
        n = self._n
        J1 = np.arange(3, n, 3)
        J2 = np.arange(4, n, 3)
        J3 = np.arange(2, n, 3)
        y = [x[j] - 2 * x[1] * np.sin(2 * np.pi * x[0] + (j + 1) * np.pi / self._n) for j in np.arange(n)]

        obj_1 = np.cos(0.5 * x[0] * np.pi) * np.cos(0.5 * x[1] * np.pi) + (2 / len(J1)) * sum(
            [4 * y[j] ** 2 - np.cos(8 * np.pi * y[j]) + 1 for j in J1])

        obj_2 = np.cos(0.5 * x[0] * np.pi) * np.sin(0.5 * x[1] * np.pi) + (2 / len(J2)) * sum(
            [4 * y[j] ** 2 - np.cos(8 * np.pi * y[j]) + 1 for j in J2])

        obj_3 = np.sin(0.5 * x[0] * np.pi) + (2 / len(J3)) * sum(
            [4 * y[j] ** 2 - np.cos(8 * np.pi * y[j]) + 1 for j in J3])

        return np.array([obj_1, obj_2, obj_3])

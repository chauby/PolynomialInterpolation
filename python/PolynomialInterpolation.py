#-*-coding:UTF-8 -*-
# ***************************************************************************
# Polynomial Interpolation
# Code Library
# ***************************************************************************
# Author: Chaobin
# Email:  chaubyZou@163.com
# Date: October 2020
# ***************************************************************************
# Language: Python
# Also available in: Matlab
# Required library: Numpy, Matplotlib
# ***************************************************************************

import numpy as np
import matplotlib.pyplot as plt

class LinearInterpolation():
    def __init__(self, name='Linear', q_via=None, t_via=None):
        """
        :param: name: string
            name of objective
        :param: q_via: N x 3 array
            given q array
        :param: t_via: N x 1 array
            given t array
        """
        super(self.__class__, self).__init__()
        self.name = name
        self.q_via = q_via
        self.t_via = t_via

        try:
            q_via.shape[1] != t_via.shape[0]
        except ValueError:
            print('The q_via and t_via must have a same length')

    def linear(self, q0, q1, t0, t1):
        """
        :param: q0: float
            the first data point
        :param: q1: float
            the second data point
        :param: t0: float
            the time of the first data point
        :param: t1: float
            the time of the second data point
        """
        try:
            abs(t0 - t1) < 1e-6
        except ValueError:
            print('t0 and t1 must be different')

        a0 = q0
        a1 = (q1 - q0)/(t1 - t0)
        return a0, a1

    def getPosition(self, t):
        """
        :param: t: float
            specified time
        :return: q: float
            output of the interpolation at time t
        """
        try:
            (t < self.t_via[0]) or (t > self.t_via[-1])
        except ValueError:
            print('The specific time error, time ranges error')

        j_array = np.where(self.t_via >= t) # find the index of t1
        j = j_array[0][0]
        if j == 0:
            i = 0
            j = 1
        else:
            i = j-1

        q = np.zeros((1, 3))

        # position
        q0 = self.q_via[i,0]
        t0 = self.t_via[i]
        q1 = self.q_via[j,0]
        t1 = self.t_via[j]
        a0, a1 = self.linear(q0, q1, t0, t1)
        q[0, 0] = a0 + a1*(t - t0)

        # velocity
        q[0, 1] = a1

        # acceleration
        q[0, 2] = 0 # for linear model, the acceleration is infinite, here we set to zero
        return q

        
class ParabolicInterpolation():
    def __init__(self, name='Parabolic', q_via=None, t_via=None):
        """
        :param: name: string
            name of objective
        :param: q_via: N x 3 array
            given q array
        :param: t_via: N x 1 array
            given t array
        """
        super(self.__class__, self).__init__()
        self.name = name
        self.q_via = q_via
        self.t_via = t_via

        try:
            q_via.shape[1] != t_via.shape[0]
        except ValueError:
            print('The q_via and t_via must have a same length')

    def parabolic(self, q0, q1, v0, v1, t0, t1, tf, qf):
        """
        :param: q0: float
            the first data point
        :param: q1: float
            the second data point
        :param: v0: float
            the velocity of the first data point
        :param: v1: float
            the velocity of the second data point
        :param: t0: float
            the time of the first data point
        :param: t1: float
            the time of the second data point
        :param: tf: float
            the time of the flex point
        :param: qf: float
            the position of the flex point
        """
        
        try:
            abs(t0 - t1) < 1e-6
        except ValueError:
                print('t0 and t1 must be different')

        try:
            ((tf <= t0) or (tf >= t1))
        except ValueError:
            print('tf must satisfy t0 < tf < t1')

        try:
            ((qf <= min(q0, q1)) or (qf >= max(q0, q1)))
        except ValueError:
            print('qf must satisfy min(q0, q1) < qf < max(q0, q1)')

        T = t1 - t0
        h = q1 - q0
        Ta = tf - t0
        Td = t1 - tf

        a0 = q0
        a1 = v0
        a2 = (2*h - v0*(T + Ta) - v1*Td)/(2*T*Ta)
        a3 = (2*q1*Ta + Td*(2*q0 + Ta*(v0 - v1)))/(2*T)
        a4 = (2*h - v0*Ta - v1*Td)/T
        a5 = -(2*h - v0*Ta - v1*(T+Td))/(2*T*Td)
        return a0, a1, a2, a3, a4, a5

    def getPosition(self, t):
        """
        :param: t: float
            specified time
        :return: q: float
            output of the interpolation at time t
        """
        try:
            (t < self.t_via[0]) or (t > self.t_via[-1])
        except ValueError:
            print('The specific time error, time ranges error')

        j_array = np.where(self.t_via >= t) # find the index of t1
        j = j_array[0][0]
        if j == 0:
            i = 0
            j = 1
        else:
            i = j-1

        q = np.zeros((1, 3))

        # get given position
        q0 = self.q_via[i,0]
        v0 = self.q_via[i,1]
        t0 = self.t_via[i]

        q1 = self.q_via[j,0]
        v1 = self.q_via[j,1]
        t1 = self.t_via[j]

        # symmetric acceleration
        tf = (t0 + t1)/2
        qf = (q0 + q1)/2

        # asymmetric acceleration, specify tf and qf by users
        # tf = ?
        # qf = ?

        a0, a1, a2, a3, a4, a5 = self.parabolic(q0, q1, v0, v1, t0, t1, tf, qf)

        if t <= tf:
            q[0, 0] = a0 + a1*(t - t0) + a2*(t-t0)**2
            q[0, 1] = a1 + 2*a2*(t - t0)
            q[0, 2] = 2*a2
        else:
            q[0, 0] = a3 + a4*(t - tf) + a5*(t-tf)**2
            q[0, 1] = a4 + 2*a5*(t - tf)
            q[0, 2] = 2*a5

        return q


class CubicInterpolation():
    def __init__(self, name='Cubic', q_via=None, t_via=None):
        """
        :param: name: string
            name of objective
        :param: q_via: N x 3 array
            given q array
        :param: t_via: N x 1 array
            given t array
        """
        super(self.__class__, self).__init__()
        self.name = name
        self.q_via = q_via
        self.t_via = t_via

        try:
            q_via.shape[1] != t_via.shape[0]
        except ValueError:
            print('The q_via and t_via must have a same length')

    def cubic(self, q0, q1, v0, v1, t0, t1):
        """
        :param: q0: float
            the first data point
        :param: q1: float
            the second data point
        :param: v0: float
            the velocity of the first data point
        :param: v1: float
            the velocity of the second data point
        :param: t0: float
            the time of the first data point
        :param: t1: float
            the time of the second data point
        """
        try:
            abs(t0 - t1) < 1e-6
        except ValueError:
                print('t0 and t1 must be different')


        T = t1 - t0
        h = q1 - q0

        a0 = q0
        a1 = v0
        a2 = (3*h - (2*v0 + v1)*T) / (T**2)
        a3 = (-2*h + (v0 + v1)*T) / (T**3)
        return a0, a1, a2, a3

    def getPosition(self, t):
        """
        :param: t: float
            specified time
        :return: q: float
            output of the interpolation at time t
        """
        try:
            (t < self.t_via[0]) or (t > self.t_via[-1])
        except ValueError:
            print('The specific time error, time ranges error')

        j_array = np.where(self.t_via >= t) # find the index of t1
        j = j_array[0][0]
        if j == 0:
            i = 0
            j = 1
        else:
            i = j-1

        q = np.zeros((1, 3))

        # get given position
        q0 = self.q_via[i,0]
        v0 = self.q_via[i,1]
        t0 = self.t_via[i]

        q1 = self.q_via[j,0]
        v1 = self.q_via[j,1]
        t1 = self.t_via[j]

        a0, a1, a2, a3 = self.cubic(q0, q1, v0, v1, t0, t1)

        q[0, 0] = a0 + a1*(t - t0) + a2*(t-t0)**2 + a3*(t - t0)**3 # position
        q[0, 1] = a1 + 2*a2*(t - t0) + 3*a3*(t - t0)**2 # velocity
        q[0, 2] = 2*a2 + 6*a3*(t - t0) # acceleration

        return q


class Polynomial5Interpolation():
    def __init__(self, name='Polynomial 5', q_via=None, t_via=None):
        """
        :param: name: string
            name of objective
        :param: q_via: N x 3 array
            given q array
        :param: t_via: N x 1 array
            given t array
        """
        super(self.__class__, self).__init__()
        self.name = name
        self.q_via = q_via
        self.t_via = t_via

        try:
            q_via.shape[1] != t_via.shape[0]
        except ValueError:
            print('The q_via and t_via must have a same length')

    def polynomial(self, q0, q1, v0, v1, acc0, acc1, t0, t1):
        """
        :param: q0: float
            the first data point
        :param: q1: float
            the second data point
        :param: v0: float
            the velocity of the first data point
        :param: v1: float
            the velocity of the second data point
        :param: acc0: float
            the acceleration of the first data point
        :param: acc1: float
            the acceleration of the second data point
        :param: t0: float
            the time of the first data point
        :param: t1: float
            the time of the second data point
        """
        try:
            abs(t0 - t1) < 1e-6
        except ValueError:
                print('t0 and t1 must be different')


        T = t1 - t0
        h = q1 - q0

        a0 = q0
        a1 = v0
        a2 = acc0/2
        a3 = (20*h - (8*v1 + 12*v0)*T - (3*acc0 - acc1)*T**2) / (2*T**3)
        a4 = (-30*h + (14*v1 + 16*v0)*T + (3*acc0 - 2*acc1)*T**2) / (2*T**4)
        a5 = (12*h - 6*(v1 + v0)*T + (acc1 - acc0)*T**2) / (2*T**5)
        return a0, a1, a2, a3, a4, a5

    def getPosition(self, t):
        """
        :param: t: float
            specified time
        :return: q: float
            output of the interpolation at time t
        """
        try:
            (t < self.t_via[0]) or (t > self.t_via[-1])
        except ValueError:
            print('The specific time error, time ranges error')

        j_array = np.where(self.t_via >= t) # find the index of t1
        j = j_array[0][0]
        if j == 0:
            i = 0
            j = 1
        else:
            i = j-1

        q = np.zeros((1, 3))

        # get given position
        q0 = self.q_via[i,0]
        v0 = self.q_via[i,1]
        acc0 = self.q_via[i,2]
        t0 = self.t_via[i]

        q1 = self.q_via[j,0]
        v1 = self.q_via[j,1]
        acc1 = self.q_via[j,2]
        t1 = self.t_via[j]

        a0, a1, a2, a3, a4, a5 = self.polynomial(q0, q1, v0, v1, acc0, acc1, t0, t1)

        q[0, 0] = a0 + a1*(t - t0) + a2*(t-t0)**2 + a3*(t - t0)**3 + a4*(t -  t0)**4 + a5*(t - t0)**5 # position
        q[0, 1] = a1 + 2*a2*(t - t0) + 3*a3*(t - t0)**2 + 4*a4*(t - t0)**3 + 5*a5*(t - t0)**4 # velocity
        q[0, 2] = 2*a2 + 6*a3*(t - t0) + 12*a4*(t- t0)**2 + 20*a5*(t - t0)**3 # acceleration

        return q



if __name__ == "__main__":
    q_given = np.array([[0, 1.6, 3.2, 2, 4, 0.2, 1.2],
                        [0, 1.0, 2.0, -2.0, -1.0, 0, 0],
                        [0, 1, 2, 3, 2, 1, 0]]).transpose()
    t_given = np.array([0, 1, 3, 4.5, 6, 8, 10]).transpose()

    # time for interpolation
    t = np.linspace(t_given[0], t_given[-1], 1000)

    #%% ************************ Linear interpolation *******************************
    linear_interpolation = LinearInterpolation('Linear', q_given, t_given)
    linear_trajectory = np.zeros((t.shape[0], 3)) # N x 3 array: position, velocity, acceleration

    for i in range(t.shape[0]):
        linear_trajectory[i,:] = linear_interpolation.getPosition(t[i])

    plt.figure(figsize=(8, 6))
    plt.subplot(3,1,1)
    plt.plot(t_given, q_given[:, 0], 'ro')
    plt.plot(t, linear_trajectory[:,0], 'k')
    plt.grid('on')
    plt.title('Linear interpolation')
    plt.xlabel('time')
    plt.ylabel('position')
    plt.xlim(t_given[0]-1, t_given[-1]+1)
    plt.ylim(min(q_given[:,0]) - 1, max(q_given[:,0]) + 1)

    plt.subplot(3,1,2)
    plt.plot(t_given, q_given[:, 1], 'ro')
    plt.plot(t, linear_trajectory[:,1], 'k')
    plt.grid('on')
    plt.xlabel('time')
    plt.ylabel('velocity')
    plt.xlim(t_given[0]-1, t_given[-1]+1)

    plt.subplot(3,1,3)
    plt.plot(t_given, q_given[:, 2], 'ro')
    plt.plot(t, linear_trajectory[:,2], 'k')
    plt.grid('on')
    plt.xlabel('time')
    plt.ylabel('acceleration')
    plt.xlim(t_given[0]-1, t_given[-1]+1)


    #%% ************************ Parabolic interpolation *******************************
    parabolic_interpolation = ParabolicInterpolation('Parabolic', q_given, t_given)
    parabolic_trajectory = np.zeros((t.shape[0], 3)) # N x 3 array: position, velocity, acceleration

    for i in range(t.shape[0]):
        parabolic_trajectory[i,:] = parabolic_interpolation.getPosition(t[i])

    plt.figure(figsize=(8, 6))
    plt.subplot(3,1,1)
    plt.plot(t_given, q_given[:, 0], 'ro')
    plt.plot(t, parabolic_trajectory[:,0], 'k')
    plt.grid('on')
    plt.title('Parabolic interpolation')
    plt.xlabel('time')
    plt.ylabel('position')
    plt.xlim(t_given[0]-1, t_given[-1]+1)
    plt.ylim(min(q_given[:,0]) - 1, max(q_given[:,0]) + 1)

    plt.subplot(3,1,2)
    plt.plot(t_given, q_given[:, 1], 'ro')
    plt.plot(t, parabolic_trajectory[:,1], 'k')
    plt.grid('on')
    plt.xlabel('time')
    plt.ylabel('velocity')
    plt.xlim(t_given[0]-1, t_given[-1]+1)

    plt.subplot(3,1,3)
    plt.plot(t_given, q_given[:, 2], 'ro')
    plt.plot(t, parabolic_trajectory[:,2], 'k')
    plt.grid('on')
    plt.xlabel('time')
    plt.ylabel('acceleration')
    plt.xlim(t_given[0]-1, t_given[-1]+1)


    #%% ************************ Cubic interpolation *******************************
    cubic_interpolation = CubicInterpolation('Cubic', q_given, t_given)
    cubic_trajectory = np.zeros((t.shape[0], 3)) # N x 3 array: position, velocity, acceleration

    for i in range(t.shape[0]):
        cubic_trajectory[i,:] = cubic_interpolation.getPosition(t[i])

    plt.figure(figsize=(8, 6))
    plt.subplot(3,1,1)
    plt.plot(t_given, q_given[:, 0], 'ro')
    plt.plot(t, cubic_trajectory[:,0], 'k')
    plt.grid('on')
    plt.title('Cubic interpolation')
    plt.xlabel('time')
    plt.ylabel('position')
    plt.xlim(t_given[0]-1, t_given[-1]+1)
    plt.ylim(min(q_given[:,0]) - 1, max(q_given[:,0]) + 1)

    plt.subplot(3,1,2)
    plt.plot(t_given, q_given[:, 1], 'ro')
    plt.plot(t, cubic_trajectory[:,1], 'k')
    plt.grid('on')
    plt.xlabel('time')
    plt.ylabel('velocity')
    plt.xlim(t_given[0]-1, t_given[-1]+1)

    plt.subplot(3,1,3)
    plt.plot(t_given, q_given[:, 2], 'ro')
    plt.plot(t, cubic_trajectory[:,2], 'k')
    plt.grid('on')
    plt.xlabel('time')
    plt.ylabel('acceleration')
    plt.xlim(t_given[0]-1, t_given[-1]+1)


    #%% *************** Polynomial of degree five interpolation *********************
    polynomial5_interpolation = Polynomial5Interpolation('Polynomial5', q_given, t_given)
    polynomial5_trajectory = np.zeros((t.shape[0], 3)) # N x 3 array: position, velocity, acceleration

    for i in range(t.shape[0]):
        polynomial5_trajectory[i,:] = polynomial5_interpolation.getPosition(t[i])

    plt.figure(figsize=(8, 6))
    plt.subplot(3,1,1)
    plt.plot(t_given, q_given[:, 0], 'ro')
    plt.plot(t, polynomial5_trajectory[:,0], 'k')
    plt.grid('on')
    plt.title('Polynomial of degree 5 interpolation')
    plt.xlabel('time')
    plt.ylabel('position')
    plt.xlim(t_given[0]-1, t_given[-1]+1)
    plt.ylim(min(q_given[:,0]) - 1, max(q_given[:,0]) + 1)

    plt.subplot(3,1,2)
    plt.plot(t_given, q_given[:, 1], 'ro')
    plt.plot(t, polynomial5_trajectory[:,1], 'k')
    plt.grid('on')
    plt.xlabel('time')
    plt.ylabel('velocity')
    plt.xlim(t_given[0]-1, t_given[-1]+1)

    plt.subplot(3,1,3)
    plt.plot(t_given, q_given[:, 2], 'ro')
    plt.plot(t, polynomial5_trajectory[:,2], 'k')
    plt.grid('on')
    plt.xlabel('time')
    plt.ylabel('acceleration')
    plt.xlim(t_given[0]-1, t_given[-1]+1)


    #%% ************************** Comparison ***************************
    plt.figure(figsize=(8, 6))
    plt.subplot(3,1,1)
    plt.plot(t_given, q_given[:, 0], 'ro', label='given pos')
    plt.plot(t, linear_trajectory[:,0], 'k', label='linear')
    plt.plot(t, parabolic_trajectory[:,0], 'b', label='parabolic')
    plt.plot(t, cubic_trajectory[:,0], 'g', label='cubic')
    plt.plot(t, polynomial5_trajectory[:,0], 'm', label='poly 5')
    plt.grid('on')
    plt.legend(loc='upper right')
    plt.title('Comparison')
    plt.xlabel('time')
    plt.ylabel('position')
    plt.xlim(t_given[0]-1, t_given[-1]+3)
    plt.ylim(min(q_given[:,0]) - 1, max(q_given[:,0]) + 1)

    plt.subplot(3,1,2)
    plt.plot(t_given, q_given[:, 1], 'ro', label='given vel')
    plt.plot(t, linear_trajectory[:,1], 'k', label='linear')
    plt.plot(t, parabolic_trajectory[:,1], 'b', label='parabolic')
    plt.plot(t, cubic_trajectory[:,1], 'g', label='cubic')
    plt.plot(t, polynomial5_trajectory[:,1], 'm', label='poly 5')
    plt.grid('on')
    plt.legend(loc='upper right')
    plt.xlabel('time')
    plt.ylabel('velocity')
    plt.xlim(t_given[0]-1, t_given[-1]+3)

    plt.subplot(3,1,3)
    plt.plot(t_given, q_given[:, 2], 'ro', label='given acc')
    plt.plot(t, linear_trajectory[:,2], 'k', label='linear')
    plt.plot(t, parabolic_trajectory[:,2], 'b', label='parabolic')
    plt.plot(t, cubic_trajectory[:,2], 'g', label='cubic')
    plt.plot(t, polynomial5_trajectory[:,2], 'm', label='poly 5')
    plt.grid('on')
    plt.legend(loc='upper right')
    plt.xlabel('time')
    plt.ylabel('acceleration')
    plt.xlim(t_given[0]-1, t_given[-1]+3)


    plt.show()

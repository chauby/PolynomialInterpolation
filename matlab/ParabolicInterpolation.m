% ***************************************************************************
% Parabolic Interpolation
% ***************************************************************************
% Author: Chaobin
% Email:  chaubyZou@163.com
% Date: October 2020
% ***************************************************************************
% Language: Matlab
% Also available in: Python
% Required library: None
% ***************************************************************************

classdef ParabolicInterpolation < handle
    properties
        name = 'parabolic interpolation';
        q_via = [];
        t_via = [];
    end

    methods
        % crate the objective
        % name: string
        % q_via: N x 3 array
        % t_via: N x 1 array
        function obj = ParabolicInterpolation(name, q_via, t_via)
            obj.name = name;
            obj.q_via = q_via;
            obj.t_via = t_via;
            if length(q_via(:,1)) ~= length(t_via)
                error('The q_via and t_via must have a same length');
            end
        end

        % parabolic interpolation with two data points
        % q0: the first data point
        % q1: the second data point
        % v0: the velocity of the first data point
        % v1: the velocity of the second data point
        % t0: the time of the first data point
        % t1: the time of the second data point
        % tf: the time of the flex point 
        % qf: the position of the flex point 
        % a0~a5: parameters
        function [a0, a1, a2, a3, a4, a5] = parabolic(obj, q0, q1, v0, v1, t0, t1, tf, qf)
            if abs(t0 - t1) < 1e-6
                error('t0 and t1 must be different');
            end

            if ((tf <= t0) || (tf >= t1))
                error('tf must satisfy t0 < tf < t1');
            end

            if ((qf <= min(q0, q1)) || (qf >= max(q0, q1)))
                error('qf must satisfy min(q0, q1) < qf < max(q0, q1)');
            end

            T = t1 - t0;
            h = q1 - q0;
            Ta = tf - t0;
            Td = t1 - tf;

            a0 = q0;
            a1 = v0;
            a2 = (2*h - v0*(T + Ta) - v1*Td)/(2*T*Ta);
            a3 = (2*q1*Ta + Td*(2*q0 + Ta*(v0 - v1)))/(2*T);
            a4 = (2*h - v0*Ta - v1*Td)/T;
            a5 = -(2*h - v0*Ta - v1*(T+Td))/(2*T*Td);
        end

        % parabolic interpolation for all data points
        % t: time
        % q: 1 x 3 array, output of the interpolation at time t
        function q = getPosition(obj, t)
            if (t < obj.t_via(1)) || (t > obj.t_via(end))
                error('The specific time error, time ranges error');
            end

            j = find(obj.t_via >= t, 1, 'first'); % find the index of t1
            if j == 1
                i = 1;
                j = 2;
            else
                i = j-1;
            end

            % get given position
            q0 = obj.q_via(i, 1);
            v0 = obj.q_via(i, 2);
            t0 = obj.t_via(i);

            q1 = obj.q_via(j, 1);
            v1 = obj.q_via(j, 2);
            t1 = obj.t_via(j);

            % symmetric acceleration
            tf = (t0 + t1)/2;
            qf = (q0 + q1)/2;

            % asymmetric acceleration, specify tf and qf by users
            % tf = ?
            % qf = ?

            [a0, a1, a2, a3, a4, a5] = obj.parabolic(q0, q1, v0, v1, t0, t1, tf, qf);

            if t <= tf
                q(1, 1) = a0 + a1*(t - t0) + a2*(t-t0)^2;
                q(1, 2) = a1 + 2*a2*(t - t0);
                q(1, 3) = 2*a2;
            else
                q(1, 1) = a3 + a4*(t - tf) + a5*(t-tf)^2;
                q(1, 2) = a4 + 2*a5*(t - tf);
                q(1, 3) = 2*a5;
            end
        end
    end
end

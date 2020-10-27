% ***************************************************************************
% Cubic Interpolation
% ***************************************************************************
% Author: Chaobin
% Email:  chaubyZou@163.com
% Date: October 2020
% ***************************************************************************
% Language: Matlab
% Also available in: Python
% Required library: None
% ***************************************************************************

classdef CubicInterpolation < handle
    properties
        name = 'cubic interpolation';
        q_via = [];
        t_via = [];
    end

    methods
        % crate the objective
        % name: string
        % q_via: N x 3 array
        % t_via: N x 1 array
        function obj = CubicInterpolation(name, q_via, t_via)
            obj.name = name;
            obj.q_via = q_via;
            obj.t_via = t_via;
            if length(q_via(:,1)) ~= length(t_via)
                error('The q_via and t_via must have a same length');
            end
        end

        % cubic interpolation with two data points
        % q0: the first data point
        % q1: the second data point
        % v0: the velocity of the first data point
        % v1: the velocity of the second data point
        % t0: the time of the first data point
        % t1: the time of the second data point
        % a0~a5: parameters
        function [a0, a1, a2, a3] = cubic(obj, q0, q1, v0, v1, t0, t1)
            if abs(t0 - t1) < 1e-6
                error('t0 and t1 must be different');
            end

            T = t1 - t0;
            h = q1 - q0;

            a0 = q0;
            a1 = v0;
            a2 = (3*h - (2*v0 + v1)*T) / (T^2);
            a3 = (-2*h + (v0 + v1)*T) / (T^3);
        end

        % cubic interpolation for all data points
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

            [a0, a1, a2, a3] = obj.cubic(q0, q1, v0, v1, t0, t1);

            q(1, 1) = a0 + a1*(t - t0) + a2*(t-t0)^2 + a3*(t - t0)^3; % position
            q(1, 2) = a1 + 2*a2*(t - t0) + 3*a3*(t - t0)^2; % velocity
            q(1, 3) = 2*a2 + 6*a3*(t - t0); % acceleration
        end
    end
end

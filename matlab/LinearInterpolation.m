% ***************************************************************************
% Linear Interpolation
% ***************************************************************************
% Author: Chaobin
% Email:  chaubyZou@163.com
% Date: October 2020
% ***************************************************************************
% Language: Matlab
% Also available in: Python
% Required library: None
% ***************************************************************************

classdef LinearInterpolation < handle
    properties
        name = 'linear interpolation';
        q_via = [];
        t_via = [];
    end

    methods
        % crate the objective
        % name: string
        % q_via: N x 3 array
        % t_via: N x 1 array
        function obj = LinearInterpolation(name, q_via, t_via)
            obj.name = name;
            obj.q_via = q_via;
            obj.t_via = t_via;
            if size(q_via(:,1)) ~= length(t_via)
                error('The q_via and t_via must have a same length');
            end
        end

        % linar interpolation with two data points
        % q0: the first data point
        % q1: the second data point
        % t0: the time of the first data point
        % t1: the time of the second data point
        % a0, a1: parameters
        function [a0, a1] = linear(obj, q0, q1, t0, t1)
            if abs(t0 - t1) < 1e-6
                error('t0 and t1 must be different');
            end
            a0 = q0;
            a1 = (q1 - q0)/(t1 - t0);
        end

        % linar interpolation for all data points
        % t: specified time
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

            % position
            q0 = obj.q_via(i,1);
            t0 = obj.t_via(i);
            q1 = obj.q_via(j,1);
            t1 = obj.t_via(j);
            [a0, a1] = obj.linear(q0, q1, t0, t1);
            q(1, 1) = a0 + a1*(t - t0);

            % velocity
            q(1, 2) = a1;

            % acceleration
            q(1, 3) = 0; % for linear model, the acceleration is infinite, here we set to zero
        end
    end
end

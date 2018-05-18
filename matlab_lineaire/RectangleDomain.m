classdef RectangleDomain
    % RectangleDomain: defines a rectangular domain object
    % A RectangleDomain object is caracterized by its lower left and
    % upper right points coordinates

    %% PROPERTIES DECLARATION

    properties
        pointA
        pointB
    end


    %% METHODS DEFINITION

    methods
        % Constructor
        function obj = RectangleDomain(varargin)
            if 2 == nargin
                obj.pointA.x = varargin{1}.x;
                obj.pointA.y = varargin{1}.y;
                obj.pointB.x = varargin{2}.x;
                obj.pointB.y = varargin{2}.y;
            elseif 4 == nargin
                obj.pointA.x = varargin{1};
                obj.pointA.y = varargin{2};
                obj.pointB.x = varargin{3};
                obj.pointB.y = varargin{4};
            else
                error('Wrong number of argument to create instance of "RectangleDomain" object')
            end
        end

        % Check if a point is contained in the domain given its coordinates
        function bool = contains_coord(obj, x, y)
            if    obj.pointA.x <= x && x <= obj.pointB.x ...
               && obj.pointA.y <= y && y <= obj.pointB.y
                bool = true;
            else
                bool = false;
            end
        end

        % Check if a point is contained in the domain given its value
        function bool = contains_point(obj, point)
            if    obj.pointA.x <= point.x && point.x <= obj.pointB.x ...
               && obj.pointA.y <= point.y && point.y <= obj.pointB.y
                bool = true;
            else
                bool = false;
            end
        end

        % Plot the rectangular domain with the indicated color
        function show(obj, color, z_pos)
            patch([obj.pointA.x, obj.pointB.x, obj.pointB.x, obj.pointA.x, obj.pointA.x], ...
                     [obj.pointA.y, obj.pointA.y, obj.pointB.y, obj.pointB.y, obj.pointA.y], ...
                     [z_pos, z_pos, z_pos, z_pos, z_pos], color);
        end
    end

end
classdef Rectangle < matlab.mixin.Copyable
    %RECTANGLE A rectangle that makes up beam cross sections

    properties
        
        x double
        y double
        w double
        h double

        xRefRect int16 % the rectangle number that this rectangle will base its x coordinate off of
        xRefFunc % function that defines how to calculate the x coordinate based off the reference rectangle
        isRelativeX % stores whether the rectangle's x coordinate will be specified based off a reference rectangle
        
        constX % holds the static the x coordinate of center of rectangle if there is no reference rectangle.
    end

    methods
        function obj = Rectangle(varargin)
            %RECTANGLE Construct a rectangle with left bottom point at
            %(x,y), width w and height h
            % If only two inputs are given the (x,y) coordinates are
            % assumed to be (0,0)
            if nargin == 2
                obj.x = 0;
                obj.y = 0;
                obj.w = varargin{1};
                obj.h = varargin{2};
            else
                obj.x = varargin{1};
                obj.y = varargin{2};
                obj.w = varargin{3};
                obj.h = varargin{4};
            end

            obj.constX = -1;

        end

        function movedRect = moveRect(obj, new_x, new_y)
            %METHOD1 Returns a value copy of this rectangle but moved to the new
            %specified coordinates
            movedRect = obj.copy();
            movedRect.x = new_x;
            movedRect.y = new_y;
            % disp("MOVED RECT: "); disp(movedRect);
        end

        function copiedRect = copyRect(obj)
            copiedRect = obj.copy();
        end

        function isRelX = isXRel(obj)
            %isXRel returns true if this rectangle's x coordinates are
            %build relative to another rectangle's coordinates
            if obj.constX == -1
                isRelX = true;
            else 
                isRelX = false;
            end
        end

        function printInfo(obj)
            disp("X: "); disp(obj.x);
            disp("y: "); disp(obj.y);
            disp("w: "); disp(obj.w);
            disp("h: "); disp(obj.h);
        end

        function area = getArea(obj)
            %getArea returns the area of this rectangle in mm^2
            area = obj.w * obj.h;
        end

        function centroidY = getCentroidY(obj)
            %getCentroid returns the centroid of this rectangle relative to
            %y = 0. 
            centroidY = obj.h / 2 + obj.y;
        end

        function rawI = getRawI(obj)
            %getI returns the raw I of rectangle. Does not include
            %contribution from parallel axis theorem
            rawI = (obj.w * obj.h^3) / 12; % Moment of inertia for a rectangle about its base
        end

        function area = getPartialArea(obj, yUpperBound)
            %getPartialArea returns the area of this rectangle that is below a
            %horizontal line yUpperBound
            if yUpperBound <= obj.y
                area = 0;
            else
                area = obj.w * (min(yUpperBound - obj.y, obj.h));
            end
        end

        function Q = getQ(obj, yUpperBound, ybar)
            %getPartialCentroid returns the Q value of area below y =
            %yUpperBound, where Q = A*d
            %d is the distance from centroid of cut rectangle to ybar
            %yUpperBound dictates how much of the rectangle matters in
            %centroid and area calculation
            if obj.y < yUpperBound
                cutRect = obj.copyRect();
                cutUpperBound = min(cutRect.y + cutRect.h, yUpperBound);
                cutRect.h = cutUpperBound - cutRect.y; % Adjust height of the cut rectangle
                
                Q = cutRect.getArea() * (ybar - cutRect.getCentroidY());
            else
                Q = 0;
            end
        end

    end
end
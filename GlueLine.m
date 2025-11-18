classdef GlueLine < matlab.mixin.Copyable
    %GLUELINE A line of glue connecting two components in a cross section

    properties
        L
        x % x coordinate of left side of glue
        y % y coordinate of left side of glue

        refRect % number of the rectangle that this glue is attached to
        yRefFunc
        xRefFunc
        lRefFunc 

    end

    methods
        function obj = GlueLine(defaultLength)
            %GLUELINE Construct a glue line that is attached to a reference
            %rectangle
            obj.L = defaultLength;
        end
        
        function newGlue = newGlueLine(obj, x, y, L)
            %newGlueLine returns a value copy of this glue line moved to
            %(x, y), and length of L
            newGlue = obj.copy();
            newGlue.x = x;
            newGlue.y = y;
            newGlue.L = L;
        end

    end
end
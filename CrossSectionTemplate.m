classdef CrossSectionTemplate < handle
    %CROSSSECTIONTEMPLATE A template of how to build a cross section
    %    MUST build cross section templates from left to right, bottom to top. 
    %    Start with rectangles, then add on glues

    properties
        rectangles {} % holds default rectangle sizes
        layers {} % information about rectangles on each vertical layer. first index is layer number from the bottom, layers(index, j) is rectangle j that may determine the height of that layer
        varDistances % holds distances that will be variable over length of bridge. layers(rect, i) is the ith variable dimension parameter for rectangle rect. [x y w h]
        glueTabs {}
        bucklingRectsHori {} % each row's first array is the horizontal rectangles' numbers. That is followed by an array of all vertical rectangles that intersect with the horizontal ones
        bucklingRectsVerti % an array holding the rectangle numbers of all vertical rectangles to be considered for case 3 buckling 
    end

    methods
        function obj = CrossSectionTemplate(baseRectangle)
            %CROSSSECTIONTEMPLATE Construct a new template with the rectangle at
            %the bottom.
            %The first rectangle will default to centered position
            obj.rectangles = {};
            obj.rectangles{1} = (baseRectangle);

            obj.layers = {};
            obj.layers{1, 1} = baseRectangle;

            baseRectangle.constX = 0;
        end

        function attachRectangleBeside(obj, varargin)
            %attachRectangleBeside attach a rectangle on the same layer as the
            %last rectangle added. Call if you want to add rectangle beside
            %the last rectangle you added.

            % Can be called as attachRectangleBeside(rect, xRefRect,
            % xRefFunc) where xRefRect is the rectangle number of the
            % reference rectangle and xRefFunc is the function that
            % calculates rect's x coord based on xRefRect
            % xRefFunc TRY NOT TO HAVE ANY CONSTANTS IN IT. If you do have
            % constants, be careful that it won't cause unindended
            % behaviour as other dimensions vary along length of beam.

            % or can call as attachRectangleBeside(rect, constX)
            % constX: the static x coordinate of the center of the rectangle

            rect = varargin{1};
            if nargin == 4
                rect.xRefRect = varargin{2};
                rect.xRefFunc = varargin{3};
            elseif nargin == 3
                rect.constX = varargin{2};
            end
            obj.layers{end, end + 1} = rect; % add the rectangle to current layer

            obj.rectangles{end + 1} = rect;
        end

        function attachRectangleAbove(obj, varargin)
            %attachRectangleAbove attach a rectangle on the layer above the
            %last rectangle added.

            % Can be called as attachRectangleAbove(rect, xRefRect,
            % xRefFunc) where xRefRect is the rectangle number of the
            % reference rectangle and xRefFunc is the function that
            % calculates rect's x coord based on xRefRect
            % xRefFunc TRY NOT TO HAVE ANY CONSTANTS IN IT. If you do have
            % constants, be careful that it won't cause unindended
            % behaviour as other dimensions vary along length of beam.

            % or can call as attachRectangleAbove(rect, constX)
            % constX: the static x coordinate of the center of the rectangle

            rect = varargin{1};
            if nargin == 4
                rect.xRefRect = varargin{2};
                rect.xRefFunc = varargin{3};
            elseif nargin == 3
                rect.constX = varargin{2};
            end
            obj.layers{end + 1, 1} = rect; % add the rectangle to new layer

            obj.rectangles{end + 1} = rect;
        end

        function attachGlue(obj, glueLine, refRect, xRefFunc, yRefFunc, lRefFunc)
            %attachGlue attaches a glue line to a specific rectangle
            %refRect on the cross section. 
            %MUST PLACE GLUE FROM BOTTOM TO TOP, LEFT TO RIGHT
            glueLine.refRect = refRect;
            glueLine.xRefFunc = xRefFunc;
            glueLine.yRefFunc = yRefFunc;
            glueLine.lRefFunc = lRefFunc;
            obj.glueTabs{end + 1} = glueLine;
        end
        
        function setBucklingRectangleHori(obj, rectNums, refRectNums)
            %rectNums are the horizontal rectangles that must be connected
            %that buckling analysis will be run on. This is to account for
            %cases where two rectangles are stacked on top of each other,
            %where the thickness is doubled for cases 1 and 2.
            %refRects are the rectangles that divide the horizontal
            %rectangle into multiple cases
            obj.bucklingRectsHori{end+1, 1} = rectNums;
            obj.bucklingRectsHori{end, 2} = refRectNums;
        end

        function setBucklingRectangleVerti(obj, rectNums)
             obj.bucklingRectsVerti = rectNums;
        end

    end
end
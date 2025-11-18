classdef CrossSection < handle
    %CROSSSECTION Builds a cross section using a template with specified
    %parameters for variable lengths
    % Places all rectangles and glue lines in their correct positions on a
    % cartesian plane with (0, 0) defined at the lower middle of the cross section

    properties
        rectangles {} % will hold the rectangles at their absolute positions relative to cartesian plane
        glueTabs {}
        CST CrossSectionTemplate % cross section template
        layers
        bucklingRectsHori {} % each row's first array is the horizontal rectangles' numbers. That is followed by an array of all vertical rectangles that intersect with the horizontal ones
        bucklingRectsVerti % an array holding the rectangle numbers of all vertical rectangles to be considered for case 3 buckling 
    end

    methods
        function obj = CrossSection(crossSectionTemplate, params)
            %CROSSSECTION Builds a cross section with the associated cross
            % section template and specific rectangle dimension parameters
            %   Detailed explanation goes here
            obj.CST = crossSectionTemplate;
            obj.layers = crossSectionTemplate.layers;
            obj.rectangles = {};
            obj.glueTabs = {};
            obj.bucklingRectsHori = crossSectionTemplate.bucklingRectsHori;
            
            % build the cross section from template
            layerHeight = 0;
            absRectNum = 1; % stores the current rectangle number we are processing

            for layerNum = 1:height(obj.layers)
                maxRectHeight = 0;
                % disp("obj.layers(layerNum)")
                % disp(obj.layers(layerNum))
                for rectNum = 1:length(obj.layers(layerNum, :)) % rectNum is the index of rectangle in that layer
                    % Process each rectangle in the current layer
                    % disp("LayerNum")
                    % disp(layerNum)
                    % disp("Rect Num: ")
                    % disp(rectNum)
                    curRect = obj.layers{layerNum, rectNum};

                    if ~isempty(curRect) % since layers may contain empty elements this checks that a rectangle exists

                        % change width and height to specified ones so that
                        % cross section can vary linearly along depth
                        varW = params(absRectNum * 2); 
                        varH = params(absRectNum * 2 - 1);
                        if varW ~= -1
                            curRect.w = varW; 
                        end
    
                        if varH ~= -1
                            curRect.h = varH;
                        end

                        if curRect.isXRel()
                            % disp("curRect = "); disp(curRect);
                            % disp("curRect.xRefRect = "); disp(curRect.xRefRect);
                            % disp("class(curRect.xRefRect) = "); disp(class(curRect.xRefRect));
                            % disp("curRect.xRefFunc result = "); disp(curRect.xRefFunc(curRect, curRect.xRefRect));
                            obj.rectangles{end + 1} = curRect.moveRect(curRect.xRefFunc(curRect, obj.rectangles{curRect.xRefRect}), layerHeight); % move the rectangle to the x coordinates specified by reference rectangle and y height determined by layer height
                        else
                            % disp("SET REC TO: "); disp(curRect.moveRect(curRect.constX - curRect.w, layerHeight))
                            obj.rectangles{end + 1} = curRect.moveRect(curRect.constX - curRect.w / 2, layerHeight); % move rectangle to its constant X defined coordinate, adjusting to center it 
                        end

                        if curRect.h > maxRectHeight
                            maxRectHeight = curRect.h; % update maxHeight of layer if this is the tallest rectangle
                        end
                        absRectNum = absRectNum + 1;
                    end

                end
                layerHeight = maxRectHeight + layerHeight;
            end

            for glueNum = 1:length(crossSectionTemplate.glueTabs)
                cur_gt = crossSectionTemplate.glueTabs{glueNum};
                cur_ref_rect = obj.rectangles{cur_gt.refRect};
                obj.glueTabs{end+1} = cur_gt.newGlueLine(cur_gt.xRefFunc(cur_gt, cur_ref_rect), ...
                                                          cur_gt.yRefFunc(cur_gt, cur_ref_rect), ...
                                                          cur_gt.lRefFunc(cur_gt, cur_ref_rect));
            end

        end

        function rects = getRects(obj)
            rects = obj.rectangles;
        end

        function glues = getGlues(obj)
            glues = obj.glueTabs;
        end

        function totalArea = getTotalArea(obj)
            totalArea = 0;
            for i = 1:length(obj.rectangles)
                totalArea = totalArea + obj.rectangles{i}.w * obj.rectangles{i}.h; % calculate area for each rectangle
            end
        end

        function smallestWidth = smallestWidthAtY(obj, y)
            % smallestWidthAtY imagines a horizontal line at y and finds the
            % smallest merged width of all rectangles intersected by that line.
            % rectangles touching or overlapping horizontally are treated as one block.
            % rectangles array must hold rectangles from left to right
            % order for this to work. This also assumes that rectangles do
            % not overlap
            
            rects = {};
            smallestWidth = inf;
        
            for i = 1:length(obj.rectangles)
                cur_rect = obj.rectangles{i};
                if cur_rect.y <= y && (cur_rect.y + cur_rect.h) >= y
                    rects{end + 1} = cur_rect;
                end
            end
            
            cur_rect = rects{1};
            cur_width = cur_rect.w;

            for i = 2:length(rects)
                next_rect = rects{i};
                if cur_rect.x + cur_rect.w >= next_rect.x % rectangles touch
                    cur_width = cur_width + next_rect.w;
                else
                    smallestWidth = min(smallestWidth, cur_width); % rectangles do not touch, update smallest width
                    cur_width = next_rect.w;
                end
                cur_rect = next_rect;
            end

            smallestWidth = min(smallestWidth, cur_width); % account for case where last rectangle is smallest 


        end

        function total = totalWidthAtY(obj, y)
            % totalWidthAtY imagines a horizontal line at y and finds the
            % total width of all rectangles intersected by that line.
            % rectangles are assumed to not overlap
            
            total = 0;

            for i = 1:length(obj.rectangles)
                cur_rect = obj.rectangles{i};
                if cur_rect.y <= y && (cur_rect.y + cur_rect.h) >= y
                    total = total + cur_rect.w;
                end
            end
        end
         
        function glueLayers = getGlueLayers(obj)
            %getGlueLayers returns a cell matrix with each row being a
            %group of glue tabs that lie on the same horizontal line.
            %iterating through rows and columns going from left to right,
            %top to bottom, will yield glue tabs from bottom to top, left
            %to right.
            glueLayers = {obj.glueTabs{1}};
            for i = 2:length(obj.glueTabs)
                if obj.glueTabs{i}.y == glueLayers{i-1}.y
                    glueLayers{end, end + 1} = obj.glueTabs{i};
                else 
                    glueLayers{end + 1, 1} = obj.glueTabs{i};
                end
            end
        end
        
        function SCritMin = calcSCritCase1(obj, E, mu)
            %calcSCritCase1 calculates the critical stress for thin plate
            %buckling case 1. E is the Young's modulus of the material, mu
            %is Poisson's ratio
            tempRects = {};
            SCritMin = inf;

            for i = 1:height(obj.bucklingRectsHori)
                for j = 1:2 % bucklingRectsHori{i,j} will be the array of either horizontal or vertical rectangle numbers
                    tempRects{i, j} = arrayfun(@(rectNum) obj.rectangles{rectNum}, obj.bucklingRectsHori{i, j}); %for each array of numbers convert to array of rectangles
                end
            end
            
            for row = 1:height(obj.bucklingRectsHori)
                horiRects = tempRects(row, 1);
                t = CrossSection.minThickness(horiRects{1});
                vertiRects = tempRects(row, 2);
                vertiRectsX = arrayfun(@(rect) rect.x, vertiRects{1});
                for i = 2:length(vertiRectsX) % iterate through each interval that is located between vertical rectangles
                    b = vertiRectsX(i) - vertiRectsX(i-1);
                    SCrit = 4*pi^2*E/12/(1-mu^2)*(t/b)^2;
                    SCritMin = min(SCritMin, SCrit);
                end
            end
        end

        function SCritMin = calcSCritCase2(obj, E, mu)
            %calcSCritCase1 calculates the critical stress for thin plate
            %buckling case 2. E is the Young's modulus of the material, mu
            %is Poisson's ratio
            tempRects = {};
            SCritMin = inf;

            for i = 1:height(obj.bucklingRectsHori)
                for j = 1:2 % bucklingRectsHori{i,j} will be the array of either horizontal or vertical rectangle numbers
                    tempRects{i, j} = arrayfun(@(rectNum) obj.rectangles{rectNum}, obj.bucklingRectsHori{i, j}); %for each array of numbers convert to array of rectangles
                end
            end
            
            for row = 1:height(obj.bucklingRectsHori)
                horiRects = tempRects(row, 1); %this is a cell array that wraps an array of the rectangles.
                horiRectsX = arrayfun(@(rect) rect.x, horiRects{1});
                horiRectsRightX = arrayfun(@(rect) rect.x + rect.w, horiRects{1});
                t = CrossSection.minThickness(horiRects{1});

                vertiRects = tempRects(row, 2);
                vertiRectsX = arrayfun(@(rect) rect.x, vertiRects{1});
                vertiRectsW = arrayfun(@(rect) rect.w, vertiRects{1});

                leftRectRightBoundary = vertiRectsX(1); % the right x boundary for the left flange of horizontal beam
                rightRectLeftBoundary = vertiRectsX(end) + vertiRectsW(end); % the left x boundary for the right flange of horizontal beam
                leftB = leftRectRightBoundary - min(horiRectsX); % have to find the x coordinate most to the left to find largest b value
                rightB = max(horiRectsRightX) - rightRectLeftBoundary;

                SCritLeft = .425*pi^2*E/12/(1-mu^2)*(t/leftB)^2;
                SCritRight = .425*pi^2*E/12/(1-mu^2)*(t/rightB)^2;

                SCritMin = min([SCritMin SCritLeft SCritRight]);
            end
        end

        function SCritMin = calcSCritCase3(obj, E, mu, centroidY)
            %calcSCritCase3 calculates the critical stress for thin plate
            %buckling case 3. E is the Young's modulus of the material, mu
            %is Poisson's ratio
            %centroidY is the ybar location for this cross section
            SCritMin = inf;

            for rectNum = 1:length(obj.rectangles)
                curRect = obj.rectangles{rectNum};
                if curRect.h > curRect.w % the rectangle is vertical if the height is greater than width, we will then run case 3 analysis on it
                    t = curRect.w;
                    rectTopY = curRect.y + curRect.h;
                    bottomOfB = max(centroidY, curRect.y); % the lower y value of b calculation
                    b = rectTopY - bottomOfB;
                    SCrit = 6*pi^2*E/12/(1-mu^2)*(t/b)^2;
                    SCritMin = min(SCrit, SCritMin);
                end
            end

            
        end
        
    end


    methods (Static)
        function m = minThickness(rects)
        % minThickness returns the minimum height of 
        % rects is an array of Rectangles
        % returns the smallest total thickness across all x.
        
            N = numel(rects);
        
            % convert to [x1 x2 h]
            intervals = zeros(N,3);
            for i = 1:N
                intervals(i,1) = rects(i).x;
                intervals(i,2) = rects(i).x + rects(i).w;
                intervals(i,3) = rects(i).h;
            end
        
            % all possible change points along x
            edges = unique([intervals(:,1); intervals(:,2)]);
        
            % initialize minimum with a large value
            m = inf;
        
            % check thickness in each subinterval
            for k = 1:length(edges)-1
                a = edges(k);
                b = edges(k+1);
                mid = (a + b) / 2;  % any interior point works
        
                % compute total thickness at this midpoint
                hsum = 0;
                for i = 1:N
                    if mid >= intervals(i,1) && mid <= intervals(i,2)
                        hsum = hsum + intervals(i,3);
                    end
                end
        
                % update minimum
                if hsum < m
                    m = hsum;
                end
            end


        end


    end
end
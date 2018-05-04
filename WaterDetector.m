classdef WaterDetector < handle
    
    properties (Access = private)
        image
        histeq
        info
        crop
        mask
        water
        index
        threshold
        lat
        lon
        shape
        ctr
    end
    
    methods
        %COSTRUCTOR
        function obj = WaterDetector()
            obj.image = [];
            obj.histeq = false;
            obj.crop.image = [];
            obj.crop.xMin = 0;
            obj.crop.xMax = 0;
            obj.crop.yMin = 0;
            obj.crop.yMax = 0;
            obj.mask = [];
            obj.info = {};
            obj.water = [];
            obj.index = '';
            obj.threshold = [];
            obj.lat = [];
            obj.lon = [];
            obj.shape = {};
            obj.ctr = {};
        end
        
        %INTERFACE
        function I = getImage(obj)
            if (obj.histeq)
                I = obj.image(:,:,1:3);
            else
                I(:,:,1) = adapthisteq(obj.image(:,:,1));
                I(:,:,2) = adapthisteq(obj.image(:,:,2));
                I(:,:,3) = adapthisteq(obj.image(:,:,3));
                
            end
        end
        
        function crop = getCrop(obj)
            crop = obj.crop;
            if (~obj.histeq && ~isempty(obj.crop.image))
                crop.image(:,:,1) = adapthisteq(obj.crop.image(:,:,1));
                crop.image(:,:,2) = adapthisteq(obj.crop.image(:,:,2));
                crop.image(:,:,3) = adapthisteq(obj.crop.image(:,:,3));
            end
        end
        
        function mask = getMask(obj)
            mask = obj.mask;
        end
        
        function water = getWater(obj)
            water = obj.water;
        end
        
        function shape = getShape(obj)
            shape = obj.shape;
        end
        
        function info = getInfo(obj)
            info = obj.info;
        end
        
        function ctr = getCTR(obj)
            ctr = obj.ctr;
        end
        
        function ctr = getWaterFromCTR(obj,descr)
            if (nargin > 1)
                ctr = {};
                k = 1;
                for i = 1 : size(obj.ctr,1)
                    if (ismember(obj.ctr(i).DESCR,descr))
                        ctr(k,1).Geometry = obj.ctr(i,1).Geometry;
                        ctr(k,1).BoundingBox = obj.ctr(i,1).BoundingBox;
                        ctr(k,1).X = obj.ctr(i,1).X;
                        ctr(k,1).Y = obj.ctr(i,1).Y;
                        ctr(k,1).LAYER = obj.ctr(i,1).LAYER;
                        ctr(k,1).DESCR = obj.ctr(i,1).DESCR;
                        ctr(k,1).id = obj.ctr(i,1).id;
                        ctr(k,1).Area = obj.ctr(i,1).Area;
                        ctr(k,1).found = obj.ctr(i,1).found;
                        k = k + 1;
                    end
                end
            end
        end
        
        %METHOD
        function readImage(obj,filename,eq)
            if (nargin > 1)
                obj.image = geotiffread(filename);
                
                temp = obj.image(:,:,1);
                obj.image(:,:,1) = obj.image(:,:,3);
                obj.image(:,:,3) = temp;
                if (exist('eq','var') && eq)
                    obj.histeq = true;
                    obj.image(:,:,1) = adapthisteq(obj.image(:,:,1));
                    obj.image(:,:,2) = adapthisteq(obj.image(:,:,2));
                    obj.image(:,:,3) = adapthisteq(obj.image(:,:,3));
                    obj.image(:,:,4) = adapthisteq(obj.image(:,:,4));
                end
                obj.info = geotiffinfo(filename);
                obj.lat = [];
                obj.lon = [];
            end
        end
        
        function crop = cropImage(obj)
            if (~isempty(obj.image))
                if (isempty(obj.mask))
                    [~, obj.crop.rect] = imcrop(obj.getImage());
                    close;
                else
                    I=obj.getImage();
                    r = I(:,:,1);
                    g = I(:,:,2);
                    b = I(:,:,3);
                    r(obj.mask) = 0;
                    g(obj.mask) = 0;
                    b(obj.mask) = 65535;
                    [~, obj.crop.rect] = imcrop(cat(3,r,g,b));
                    close;
                end
                if (~isempty(obj.crop.rect))
                    obj.crop.xMin = round(obj.crop.rect(1));
                    width = round(obj.crop.rect(3));
                    height = round(obj.crop.rect(4));
                    obj.crop.yMin = round(obj.crop.rect(2));
                    obj.crop.xMax = obj.crop.xMin + (width - 1);
                    obj.crop.yMax = obj.crop.yMin + (height - 1);
                    obj.crop.image = obj.image(obj.crop.yMin:obj.crop.yMax,...
                        obj.crop.xMin:obj.crop.xMax,:);
                    obj.water = [];
                    obj.index = '';
                    obj.threshold = [];
                    obj.lat = [];
                    obj.lon = [];
                    crop = obj.crop;
                end
            end
        end
        
        function removeCrop(obj)
            obj.crop.image = [];
            obj.crop.rect = [];
            obj.crop.xMin = 0;
            obj.crop.xMax = 0;
            obj.crop.yMin = 0;
            obj.crop.xMax = 0;
            if (~isempty(obj.index) || ~isempty(obj.threshold))
                obj.calcNdvi(obj.index, obj.threshold);
            end
            obj.water = [];
            obj.index = '';
            obj.threshold = [];
            obj.lat = [];
            obj.lon = [];
        end
        
        function Imm = calcNdvi(obj, index, threshold)
            if (~isempty(obj.image) && nargin >1)
                if (~isempty(obj.crop.image))
                    I = obj.crop.image;
                else
                    I = obj.image;
                end
                I = double(I);
                index = lower(index);
                tMin = -1;
                tMax = 1;
                switch index
                    case 'ndvi'
                        Imm= (I(:,:,4) - I(:,:,1)) ./ (I(:,:,4) + I(:,:,1));
                        tMin = -0.7;
                        tMax = -0.3;
                    case 'ndwi2'
                        Imm = (I(:,:,2) - I(:,:,4)) ./ (I(:,:,2) + I(:,:,4));
                        tMin = 0.7;
                    case 'ndwi2m'
                        Imm = (I(:,:,4) - I(:,:,2)) ./ (I(:,:,4) + I(:,:,2));
                        tMax = -0.7;
                    otherwise
                        return;
                end
                obj.index = index;
                obj.threshold = threshold;
                obj.water = false(size(Imm));
                if (nargin > 2)
                    if (size(threshold,2) == 1)
                        if (threshold>=-1 && threshold<=1)
                            tMin = threshold;
                        end
                    elseif (size(threshold,2) > 1)
                        if (threshold(1)>=-1 && threshold(1)<=1)
                            tMin = threshold(1);
                        end
                        if (threshold(2)>=threshold(1) && threshold(2)<=1)
                            tMax = threshold(2);
                        end
                    end
                end
                obj.water = Imm > tMin & Imm < tMax;
                Imm(~obj.water) = 0;
            end
        end
        
        function image = showMask(obj)
            image = [];
            if (~isempty(obj.mask))
                if (isempty(obj.crop.image))
                    r = obj.image(:,:,1);
                    g = obj.image(:,:,2);
                    b = obj.image(:,:,3);
                    r(obj.mask) = 0;
                    g(obj.mask) = 0;
                    b(obj.mask) = 65535;
                    image = cat(3,r,g,b);
                else
                    r = obj.crop.image(:,:,1);
                    g = obj.crop.image(:,:,2);
                    b = obj.crop.image(:,:,3);
                    r(obj.mask(obj.crop.yMin:obj.crop.yMax,obj.crop.xMin:obj.crop.xMax)) = 0;
                    g(obj.mask(obj.crop.yMin:obj.crop.yMax,obj.crop.xMin:obj.crop.xMax)) = 0;
                    b(obj.mask(obj.crop.yMin:obj.crop.yMax,obj.crop.xMin:obj.crop.xMax)) = 65535;
                    image = cat(3,r,g,b);
                end
            end
        end
        
        function shape2mask(obj, shape)
            if (nargin > 1 && ~isempty(obj.image) && ~isempty(shape))
                if (isempty(obj.mask))
                    obj.mask = false(size(obj.image(:,:,1)));
                end
                [x,y] = pixcenters(obj.info);
                [X,Y] = meshgrid(x,y);
                for i=1:size(shape)
                    rx = shape(i).X(1:end-1);
                    ry = shape(i).Y(1:end-1);
                    obj.mask = obj.mask | inpolygon(X,Y,rx,ry);
                end
            end
        end
        
        function removeShapeMask(obj, shape)
            if (nargin > 1 && ~isempty(obj.mask) && ~isempty(shape))
                [x,y] = pixcenters(obj.info);
                [X,Y] = meshgrid(x,y);
                for i = 1:size(shape)
                    rx = shape(i).X(1:end-1);
                    ry = shape(i).Y(1:end-1);
                    obj.mask(inpolygon(X,Y,rx,ry)) = false;
                end
                if (isempty(obj.mask(obj.mask)))
                    obj.mask = [];
                end
            end
        end
        
        function applyMask(obj)
            if (~isempty(obj.image))
                if (isempty(obj.mask))
                    obj.mask = false(size(obj.image(:,:,1)));
                end
                m = roipoly(obj.showMask());
                close;
                if (~isempty(m))
                    if (isempty(obj.crop.image))
                        obj.mask = obj.mask | m;
                    else
                        obj.mask(obj.crop.yMin:obj.crop.yMax,obj.crop.xMin:obj.crop.xMax) = ...
                            obj.mask(obj.crop.yMin:obj.crop.yMax,obj.crop.xMin:obj.crop.xMax) | m;
                    end
                end
            end
        end
        
        function removeMask(obj)
            if (~isempty(obj.mask))
                poly = roipoly(obj.showMask());
                close;
                if (~isempty(poly))
                    if (isempty(obj.crop.image))
                        obj.mask(poly) = false;
                    else
                        m = false(size(obj.image(:,:,1)));
                        m(obj.crop.yMin:obj.crop.yMax,obj.crop.xMin:obj.crop.xMax) = poly;
                        obj.mask(m) = false;
                    end
                    if (isempty(obj.mask(obj.mask)))
                        obj.mask = [];
                    end
                end
            end
        end
        
        function maskFromCTR(obj, ctr, descr)
            if (nargin > 2 && ~isempty(obj.image))
                if (isempty(obj.mask))
                    obj.mask = false(size(obj.image(:,:,1)));
                end
                [x,y] = pixcenters(obj.info);
                [X,Y] = meshgrid(x,y);
                for i=1:size(ctr)
                    center = round(map2pix(obj.info.RefMatrix,ctr(i).BoundingBox(:,1),ctr(i).BoundingBox(:,2)));
                    xMin = center(1,2);
                    xMax = center(2,2);
                    yMin = center(2,1);
                    yMax = center(1,1);
                    if(ismember(ctr(i).DESCR,descr) && obj.compareBoundingBox(ctr(i).BoundingBox,obj.info.BoundingBox) && ...
                            sum(obj.image(yMin,xMin,:)~=909 & obj.image(yMax,xMax,:)~=909)==4)
                        rx = ctr(i).X(1:end-1);
                        ry = ctr(i).Y(1:end-1);
                        obj.mask = obj.mask | inpolygon(X,Y,rx,ry);
                    end
                end
            end
        end
        
        function removeMaskFromCTR(obj,ctr,descr)
            if (nargin > 2 && ~isempty(obj.mask))
                [x,y] = pixcenters(obj.info);
                [X,Y] = meshgrid(x,y);
                for i = 1:size(ctr)
                    center = round(map2pix(obj.info.RefMatrix,ctr(i).BoundingBox(:,1),ctr(i).BoundingBox(:,2)));
                    xMin = center(1,2);
                    xMax = center(2,2);
                    yMin = center(2,1);
                    yMax = center(1,1);
                    if(ismember(ctr(i).DESCR,descr) && obj.compareBoundingBox(ctr(i).BoundingBox,obj.info.BoundingBox) && ...
                            sum(obj.image(yMin,xMin,:)~=909 & obj.image(yMax,xMax,:)~=909)==4)
                        rx = ctr(i).X(1:end-1);
                        ry = ctr(i).Y(1:end-1);
                        obj.mask(inpolygon(X,Y,rx,ry)) = false;
                    end
                end
                if (isempty(obj.mask(obj.mask)))
                    obj.mask = [];
                end
            end
        end
        
        function [num, shape] = createShape(obj, Amin, Amax)
            if (~isempty(obj.water))
                obj.shape = {};
                bw = obj.water;
                if (~isempty(obj.crop.image) && ~isempty(obj.mask))
                    bw(obj.mask(obj.crop.yMin:obj.crop.yMax,...
                        obj.crop.xMin:obj.crop.xMax)) = false;
                elseif (isempty(obj.crop.image) && ~isempty(obj.mask))
                    bw(obj.mask) = false;
                end
                cc = bwconncomp(bw, 8);
                j = 1;
                temp = false(size(bw));
                if (isempty(obj.lat) || isempty(obj.lon))
                    if (~isempty(obj.crop.image))
                        [obj.lat, obj.lon] = obj.getLatLon([obj.crop.yMin,obj.crop.yMax,...
                            obj.crop.xMin,obj.crop.xMax]);
                    else
                        [obj.lat, obj.lon] = obj.getLatLon();
                    end
                end
                for i = 1:cc.NumObjects
                    temp(cc.PixelIdxList{i}) = true;
                    coord = [obj.lat(cc.PixelIdxList{i}) obj.lon(cc.PixelIdxList{i})];
                    maxLat = max(coord(:,1));
                    p1 = obj.lon(obj.lat == maxLat);
                    minLat = min(coord(:,1));
                    p2 = obj.lon(obj.lat == minLat);
                    maxLon = max(coord(:,2));
                    p3 = obj.lat(obj.lon == maxLon);
                    minLon = min(coord(:,2));
                    p4 = obj.lat(obj.lon == minLon);
                    [x, y] = projfwd(obj.info, [maxLat, p3, minLat, p4, maxLat], [p1, maxLon, p2, minLon, p1]);
                    area = polyarea(x,y);
                    checkAmin = area > 3;
                    checkAmax = true;
                    if (exist('Amin', 'var') && Amin >= 0)
                        checkAmin = area > Amin;
                    end
                    if (exist('Amax', 'var') && Amax >= Amin)
                        checkAmax = area < Amax;
                    end
                    if (size(coord,1) > 4 && area > 0 && checkAmin && checkAmax)
                        obj.shape(j).Geometry = deal('Polygon');
                        obj.shape(j).id = j;
                        obj.shape(j).X = x;
                        obj.shape(j).Y = y;
                        obj.shape(j).BoundingBox = [min(x) min(y); max(x) max(y)];
                        obj.shape(j).Description = 'vasche trovate';
                        obj.shape(j).Area = area;
                        obj.shape(j).isNew = true;
                        k = 1;
                        found = false;
                        while (k <= size(obj.ctr,1) && ~found)
                            found = obj.compareBoundingBox(obj.shape(j).BoundingBox,obj.ctr(k).BoundingBox);
                            k= k + 1;
                        end
                        if (found)
                            obj.shape(j).BoundingBox = obj.ctr(k-1).BoundingBox;
                            obj.shape(j).Area = obj.ctr(k-1).Area;
                            obj.shape(j).X = obj.ctr(k-1).X;
                            obj.shape(j).Y = obj.ctr(k-1).Y;
                            obj.shape(j).isNew = false;
                            obj.ctr(k-1).found = true;
                        end
                        j = j + 1;
                    end
                    temp(cc.PixelIdxList{i}) = false;
                end
                num = j - 1;
                shape = obj.shape;
            end
        end
        
        function [element, area] = showElement(obj,id,z)
            if (nargin>1 && ~isempty(obj.shape))
                shapeElement = obj.shape([obj.shape.id] == id);
                zoom = 100;
                if (nargin > 2)
                    zoom = z;
                end
                center = round(map2pix(obj.info.RefMatrix,shapeElement.BoundingBox(:,1),shapeElement.BoundingBox(:,2)));
                xMin = ((center(1,2)-zoom)>=0)*(center(1,2)-zoom) + ((center(1,2)-zoom)<0)*0;
                xMax = ((center(2,2)+zoom)<=obj.info.Width)*(center(2,2)+zoom) + ((center(2,2)+zoom)>obj.info.Width)*obj.info.Width;
                yMin = ((center(2,1)-zoom)>=0)*(center(2,1)-zoom) + ((center(2,1)-zoom)<0)*0;
                yMax = ((center(1,1)+zoom)<=obj.info.Height)*(center(1,1)+zoom) + ((center(1,1)+zoom)>obj.info.Height)*obj.info.Height;
                element = obj.image(yMin : yMax,xMin : xMax,1:3);
                area = shapeElement.Area;
            end
        end
        
        function [element, area] = showCTRElement(obj,id,z)
            if (nargin>1 && ~isempty(obj.ctr))
                shapeElement = obj.ctr([obj.ctr.id] == id);
                zoom = 100;
                if (nargin > 2)
                    zoom = z;
                end
                center = round(map2pix(obj.info.RefMatrix,shapeElement.BoundingBox(:,1),shapeElement.BoundingBox(:,2)));
                xMin = ((center(1,2)-zoom)>=0)*(center(1,2)-zoom) + ((center(1,2)-zoom)<0)*0;
                xMax = ((center(2,2)+zoom)<=obj.info.Width)*(center(2,2)+zoom) + ((center(2,2)+zoom)>obj.info.Width)*obj.info.Width;
                yMin = ((center(2,1)-zoom)>=0)*(center(2,1)-zoom) + ((center(2,1)-zoom)<0)*0;
                yMax = ((center(1,1)+zoom)<=obj.info.Height)*(center(1,1)+zoom) + ((center(1,1)+zoom)>obj.info.Height)*obj.info.Height;
                element = obj.image(yMin : yMax,xMin : xMax,1:3);
                area = shapeElement.Area;
            end
        end
        
        function removeElement(obj,id)
            if (nargin>1 && ~isempty(obj.shape))
                shapeElement = obj.shape([obj.shape.id] == id);
                if (isempty(obj.mask))
                    obj.mask = false(size(obj.image(:,:,1)));
                end
                if (~isempty(shapeElement))
                    center = round(map2pix(obj.info.RefMatrix,shapeElement.BoundingBox(:,1),shapeElement.BoundingBox(:,2)));
                    xMin = center(1,2);
                    xMax = center(2,2);
                    yMin = center(2,1);
                    yMax = center(1,1);
                    obj.mask(yMin : yMax,xMin : xMax) = true;
                    obj.shape([obj.shape.id] == id) = [];
                    for i = id : (obj.shape(end).id-1)
                        obj.shape(i).id = obj.shape(i).id - 1;
                    end
                end
            end
        end
        
        function applyMorphClose(obj)
            if (~isempty(obj.water))
                obj.water = bwmorph(obj.water, 'close');
            end
        end
        
        function applyMorphClean(obj)
            if (~isempty(obj.water))
                obj.water = bwmorph(obj.water, 'clean');
            end
        end
        
        function [num, ctr] = importCTR(obj,shape,descr)
            if (nargin > 2 && ~isempty(obj.info))
                obj.ctr = {};
                ctr = {};
                inCrop = 0;
                j = 1;
                z = 1;
                black = 909;
                if (~obj.histeq)
                    black = 0;
                end
                for i = 1 : size(shape,1)
                    center = round(map2pix(obj.info.RefMatrix,shape(i).BoundingBox(:,1),shape(i).BoundingBox(:,2)));
                    xMin = center(1,2);
                    xMax = center(2,2);
                    yMin = center(2,1);
                    yMax = center(1,1);
                    if (obj.compareBoundingBox(shape(i).BoundingBox,obj.info.BoundingBox) && ...
                            sum(obj.image(yMin,xMin,:)~=black & obj.image(yMax,xMax,:)~=black)==4)
                        ctr(z,1).Geometry = shape(i).Geometry;
                        ctr(z,1).BoundingBox = shape(i).BoundingBox;
                        ctr(z,1).X = shape(i).X;
                        ctr(z,1).Y = shape(i).Y;
                        ctr(z,1).LAYER = shape(i).LAYER;
                        ctr(z,1).DESCR = shape(i).DESCR;
                        ctr(z,1).id = z;
                        ctr(z,1).Area = polyarea(ctr(z,1).X(~isnan(ctr(z,1).X)),ctr(z,1).Y(~isnan(ctr(z,1).Y)));
                        ctr(z,1).found = false;
                        if (ismember(shape(i).DESCR,descr))
                            obj.ctr(j,1).Geometry = ctr(z,1).Geometry;
                            obj.ctr(j,1).BoundingBox = ctr(z,1).BoundingBox;
                            obj.ctr(j,1).X = ctr(z,1).X;
                            obj.ctr(j,1).Y = ctr(z,1).Y;
                            obj.ctr(j,1).LAYER = ctr(z,1).LAYER;
                            obj.ctr(j,1).DESCR = ctr(z,1).DESCR;
                            obj.ctr(j,1).id = j;
                            obj.ctr(j,1).Area = ctr(z,1).Area;
                            if (~isempty(obj.shape))
                                y = 1;
                                while (y<= size(obj.shape,2) && ~ctr(z,1).found)
                                    ctr(z,1).found = obj.compareBoundingBox(obj.shape(1,y).BoundingBox,ctr(z).BoundingBox);
                                    y = y + 1;
                                end
                                if (ctr(z,1).found)
                                    obj.shape(1,y-1).BoundindBox = ctr(z,1).BoundingBox;
                                    obj.shape(1,y-1).Area = ctr(z,1).Area;
                                    obj.shape(1,y-1).X = ctr(z,1).X;
                                    obj.shape(1,y-1).Y = ctr(z,1).Y;
                                    obj.shape(1,y-1).isNew = ~ctr(z,1).found;
                                end
                            end
                            obj.ctr(j,1).found = ctr(z,1).found;
                            j = j + 1;
                            if (~isempty(obj.crop.image))
                                cropBoundingBox = pix2map(obj.info.RefMatrix,[obj.crop.yMin, obj.crop.xMin]',[obj.crop.xMax, obj.crop.yMax]');
                                if (obj.compareBoundingBox(shape(i).BoundingBox,cropBoundingBox))
                                    inCrop = inCrop + 1;
                                end
                            end
                        end
                        z = z + 1;
                    end
                    if (isempty(obj.crop.image))
                        num = size(obj.ctr,1);
                    else
                        num = [size(obj.ctr,1), inCrop];
                    end
                end
            end
        end
        
        function ctr = cropCTR(obj,descr)
            ctr = {};
            if (~isempty(obj.crop.image))
                j = 1;
                for i = 1:size(obj.ctr,1)
                    cropBoundingBox = pix2map(obj.info.RefMatrix,[obj.crop.yMax; obj.crop.yMin],[obj.crop.xMin; obj.crop.xMax]);
                    check = true;
                    if (nargin>1)
                        check = ismember(string(obj.ctr(i).DESCR),string(descr));
                    end
                    if (obj.compareBoundingBox(obj.ctr(i).BoundingBox,cropBoundingBox) && check)
                        ctr(j,1).Geometry = obj.ctr(i,1).Geometry;
                        ctr(j,1).BoundingBox = obj.ctr(i,1).BoundingBox;
                        ctr(j,1).X = obj.ctr(i,1).X;
                        ctr(j,1).Y = obj.ctr(i,1).Y;
                        ctr(j,1).LAYER = obj.ctr(i,1).LAYER;
                        ctr(j,1).DESCR = obj.ctr(i,1).DESCR;
                        ctr(j,1).id = obj.ctr(i,1).id;
                        ctr(j,1).Area = obj.ctr(i,1).Area;
                        ctr(j,1).found = obj.ctr(i,1).found;
                        j = j + 1;
                    end
                end
            end
        end
        
        function addCTRElementToShape(obj,id)
            if (nargin>1)
                element = obj.ctr([obj.ctr.id] == id);
                if (~isempty(element) && ~element.found)
                    i = size(obj.shape,2) + 1;
                    obj.shape(1,i).Geometry = element.Geometry;
                    obj.shape(1,i).X = element.X;
                    obj.shape(1,i).Y = element.Y;
                    obj.shape(1,i).BoundingBox = element.BoundingBox;
                    obj.shape(1,i).Description = obj.shape(1).Description;
                    obj.shape(1,i).Area = element.Area;
                    obj.shape(1,i).isNew = false;
                    obj.shape(1,i).id = i;
                    obj.ctr([obj.ctr.id] == id).found = true;
                end
            end
        end
        
        function removeToCTR(obj,id)
            if (nargin>1)
                element = obj.ctr([obj.ctr.id] == id);
                if (~isempty(element))
                    obj.ctr([obj.ctr.id] == id) = [];
                    for i = id : (obj.ctr(end).id-1)
                        obj.ctr(i).id = obj.ctr(i).id - 1;
                    end
                end
            end
        end
        
        function removeCTR(obj)
            obj.ctr = {};
        end
        
    end
    
    methods (Access = private)
        function [lat, lon] = getLatLon(obj,rect)
            if (nargin > 0)
                yMin = 1;
                yMax = obj.info.Height;
                xMin = 1;
                xMax = obj.info.Width;
                if (nargin > 1)
                    yMin = rect(1);
                    yMax = rect(2);
                    xMin = rect(3);
                    xMax = rect(4);
                end
                [rows,cols] = meshgrid(yMin:yMax,xMin:xMax);
                [x,y] = pix2map(obj.info.RefMatrix, rows, cols);
                [lat,lon] = projinv(obj.info,x,y);
                lat = lat';
                lon = lon';
            end
        end
        
        function flag = compareBoundingBox(~,BoundingBox1,BoundingBox2)
            if (nargin > 2)
                if (BoundingBox1(:,1) >= BoundingBox2(1,1) & ...
                        BoundingBox1(:,1) <= BoundingBox2(2,1) & ...
                        BoundingBox1(:,2) >= BoundingBox2(1,2) & ...
                        BoundingBox1(:,2) <= BoundingBox2(2,2))
                    flag = true;
                else
                    flag = false;
                end
            end
        end
    end
end

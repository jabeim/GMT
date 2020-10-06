% DataUnit < handle
% Copyright (c) 2012-2020 Advanced Bionics. All rights reserved.
classdef DataUnit < handle
    properties (SetAccess = private)
        parent;
        ID;
        type;  % 'INPUT' or 'OUTPUT'
        connection = Connection.empty(0,0);
    end
    properties (Access = private)
        data = [];
    end
    
    methods (Access = {?ProcUnit, ?Strategy})
        function mod = setData(obj, data)            
            mod = ~isequal(data, obj.data);
            if mod
                obj.data = data;
            end         
        end
        
        function resetData(obj)
            obj.data = [];
        end
    end
    
    methods
        function obj = DataUnit(parent, type, ID)
            obj.parent = parent;
            obj.type = type;
            obj.ID = ID;
        end
        
        function addConnection(obj, con)
            assert(isa(con,'Connection'), 'con must be an instance of class Connection.')
            if strcmp(obj.type,'INPUT') && ~isempty(obj.connection)
                error('Attempt to create multiple connections to one input DataUnit.')
            end
            obj.connection(end+1) = con;
        end
        
        function data = getData(obj)
           data = obj.data;
        end
        
        function n = getNumberConnections(obj)
            n = length(obj.connection); 
        end
        
        function e = dataIsEmpty(obj)
            e = isempty(obj.data);
        end
    end

end

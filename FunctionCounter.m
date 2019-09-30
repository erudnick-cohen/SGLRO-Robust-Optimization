classdef FunctionCounter < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        counter
    end
    
    methods
        function obj = FunctionCounter()
            obj.counter = 0;
        end
        function out = count(this,in)
            this.counter = this.counter+1;
            out = in;
        end
    end
    
end


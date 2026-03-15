function [value,units]= KBoltzmann_ev(varargin)
% Boltzmann constant in units eV as default

units_index = 1;
format_units = {'eV','J','erg'};
if numel(varargin) == 1
    if ~isempty(varargin{1})
        units_index_temp = find(cellfun(@(x) contains(varargin{1},x,'IgnoreCase',true),format_units,'UniformOutput',true));
        if ~isempty(units_index_temp)
            units_index = units_index_temp;
        end
    end
end

units = format_units{units_index};

if strcmp(units,'eV')
    value = 8.617333262145e-5;
elseif strcmp(units,'J')
    value = 1.380649e-23;
elseif strcmp(units,'erg')
    value = 1.380649e-16;
end


end
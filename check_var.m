function res=check_var(name)

% res=check_var(name);
%
% AS, UWM, 2007

res=evalin('caller', ['exist(''' name ''',''var'') && ~isempty(' name ')']);
function value = p44_option(s,name,fallback)
%P44_OPTION Read a nonempty optional configuration field.
value=fallback;
if isfield(s,name) && ~isempty(s.(name)); value=s.(name); end
end

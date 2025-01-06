function row = dd_store_input(inp1,i,data);
    
data     = rmfield(data, {'IntlInd', 'HospInd', 'EdInd', 'B', 'C', 'x_schc', 'x_econ', 'x_elim', 'tvec', 'vly', 'vsy'});
data.NNs = data.NNs(1:45);
fields   = fieldnames(data);

row         = struct;
row.igroup  = inp1;
row.country = i;
for f = 1:length(fields)
    row.(fields{f}) = data.(fields{f})(:).';
end
row         = struct2table(row);

end
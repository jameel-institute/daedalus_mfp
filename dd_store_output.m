function dd_store_output(array,title,taskdir)

T = array2table(array);

T.Properties.VariableNames = {'country','tdur',...
                              'vlyl_1','vlyl_2','vlyl_3','vlyl_4', ...
                              'gdpl_1','gdpl_2','gdpl_3','gdpl_4','gdpl_5','gdpl_6','gdpl_7','gdpl_8','gdpl_9','gdpl_10', ...
                              'vsyl',
                              'attk_1','attk_2','attk_3','attk_4',
                              'hosp_1','hosp_2','hosp_3','hosp_4'};

writetable(T,strcat(taskdir,title,'.csv'));

end
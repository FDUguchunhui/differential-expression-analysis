library(xlsx)

# get result of all comparisons
{
sheetname <- names(res_list_all)
for(i in 1:length(res_list_all)){
  write.xlsx(x = res_list_all[[i]],
             file = 'output/result_kinetics.xlsx',
             sheet = sheetname[i], append = T)
}
}
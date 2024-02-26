library(readxl)

read_supplementary_table_S1_and_save_csv<-function(sheet, name)
{
  SheetDf<-read_excel("./data/Table_S1_defence_occurrence_in_datasets.xlsx", sheet=sheet)
  write.table(SheetDf, file =paste0("./data/",name,"_filtered.csv"),
              quote=F, sep=",", row.names = F)
}

##
read_supplementary_table_S1_and_save_csv(1,"ecoli")
read_supplementary_table_S1_and_save_csv(2,"baci")
read_supplementary_table_S1_and_save_csv(3,"burk")
read_supplementary_table_S1_and_save_csv(4,"enter")
read_supplementary_table_S1_and_save_csv(5,"pseu")

# cd /Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/02_Visualisation_Folder/Mav_Plots
# python3 Sian_MAV_plots.py

import os, sys, subprocess

gene_list=open('FstH_1Percent.txt', 'r') # temp4.txt
Rfile=open('sian_plots.R', 'w+')
Rfile.write('#Rscript sian_plots.R'+'\n'+'setwd("/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/02_Visualisation_Folder/Mav_Plots/")' +'\n'+
		'library(data.table)' +'\n'+
		'library(ggplot2)' +'\n'+
		'df<-fread("UK_dips_and_tets_MAV.table")' +'\n'+
		'df<-transform(df, tet_freq = tet_AC / tet_AN)' +'\n'+
		'df<-transform(df, dip_freq = dip_AC / dip_AN)' +'\n'+
		'df<-na.omit(df)'+'\n'+
		'df<-transform(df, AFD = abs(dip_freq - tet_freq))' +'\n'+
		'an<-fread("/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/02_Visualisation_Folder/Mav_Plots/gene_orientation_file.txt")' +'\n'+
		'genes<-read.table("FstH_1Percent.txt",h=T)' +'\n'+
		'MAVs<-read.table("british_1percent_outliers_3.5_scores.txt",h=T)'+'\n')

for count, gene in enumerate(gene_list):
	gene=gene.replace('\n', '')
	Rfile.write('ge="'+gene+'"'+'\n'+
		'an1<-subset(an,substr(gene,2,15) %'+'in%'+' substr(ge,2,15))'+'\n'+
		'df0<-subset(df,scaff %'+'in%'+' an1$scaff & pos > (an1$start - 60000)  & pos < (an1$end + 60000))'+'\n'+
		'an2<-subset(an,scaff %'+'in%'+' df0$scaff & an$start >= df0$pos[1] & an$end <= df0$pos[nrow(df0)])'+'\n'+
		'MAVs1<-subset(MAVs, scaff %'+'in%'+' df0$scaff & pos >= an1$start & pos <= an1$end)'+'\n'+
		'df2<-subset(df0, scaff %'+'in% MAVs1$scaff & pos %'+'in% MAVs1$pos) # this is the MAV SNPs in the candidate gene!'+'\n'+
		'df1<-subset(df0, !(pos %'+'in%'+' df2$pos))'+'\n'+
		'if (an1$orient %'+'in% "+") {'+'\n'+
		'    a<-ggplot() + '+'\n'+
		'    geom_point(data = df1, aes(x=as.numeric(as.character(pos)), y=as.numeric(as.character(AFD))),size = 3, alpha = 0.5) + # plots SNPs'+'\n'+
		'    try((geom_point(data = df2, aes(x=as.numeric(as.character(pos)), y=as.numeric(as.character(AFD))),size = 4, colour = "red", alpha = 0.5)), silent = TRUE) + # plots MAV SNPs'+'\n'+
		'    geom_rect(aes(xmin = an1$start, xmax = an1$end, ymin = -Inf, ymax = Inf),alpha=0.1) + '+'\n'+
		'    try((geom_segment(data = an2, aes(x = an2$start, y = 1+0.02, xend = an2$end, yend = 1+0.02), colour = "grey",size=1)), silent = TRUE) + #plots other genes'+'\n'+
		'    geom_segment(data = df1, aes(x = an1$start, y = 1+0.02, xend = an1$end, yend = 1+0.02), colour = "red",size=1.5, arrow = arrow(length = unit(0.3, "cm"),ends = "last")) + # plots gene of intrest'+'\n'+
		'    labs(x=an1$scaff, y =expression("Absolute AFD"),title = paste(ge)) + # makes labels and title'+'\n'+
		'    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "grey"),axis.text = element_text(size=18),axis.title.x = element_text(size=24),axis.title.y = element_text(size=24)) + #makes axis etc'+'\n'+
		'    expand_limits(y=1) '+'\n'+
		'} else'+'\n'+
		'{ a<-ggplot() + '+'\n'+
		'  geom_point(data = df1, aes(x=as.numeric(as.character(pos)), y=as.numeric(as.character(AFD))),size = 3, alpha = 0.5) + # plots SNPs'+'\n'+
		'  try((geom_point(data = df2, aes(x=as.numeric(as.character(pos)), y=as.numeric(as.character(AFD))),size = 4, colour = "red", alpha = 0.5)), silent = TRUE) + # plots MAV SNPs'+'\n'+
		'  geom_rect(aes(xmin = an1$start, xmax = an1$end, ymin = -Inf, ymax = Inf),alpha=0.1) + '+'\n'+
		'  try((geom_segment(data = an2, aes(x = an2$start, y = 1+0.02, xend = an2$end, yend = 1+0.02), colour = "grey",size=1)), silent = TRUE) + #plots other genes'+'\n'+
		'  geom_segment(data = df1, aes(x = an1$start, y = 1+0.02, xend = an1$end, yend = 1+0.02), colour = "red",size=1.5, arrow = arrow(length = unit(0.3, "cm"),ends = "first")) + # plots gene of intrest'+'\n'+
		'  labs(x=an1$scaff, y =expression("Absolute AFD"),title = paste(ge)) + # makes labels and title'+'\n'+
		'  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "grey"),axis.text = element_text(size=18),axis.title.x = element_text(size=24),axis.title.y = element_text(size=24)) + #makes axis etc'+'\n'+
		'  expand_limits(y=1) '+'\n'+
		'}'+'\n'+
		'pdf(paste(ge,".pdf",sep=""),height = 4,width = 8,pointsize = 12)'+'\n'+
		'print(a)'+'\n'+
		'dev.off()'+'\n')

os.system('Rscript sian_plots.R')
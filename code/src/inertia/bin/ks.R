x<-read.table(commandArgs()[4],sep="\t") # read input from "input" into table
y1=x[,1] #create 2 vectors
y2=x[,2]
ks.test(y1,y2) #do kolmogorov smirnov test to see if vector y1 and y2 are from the same distribution

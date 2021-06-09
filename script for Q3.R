
#######NNET

data.train<-read.table("entire.txt",sep = "\t",header = T,row.names = 1)
data.train$cluster=as.factor(data.train$cluster)
r=1/max(abs(data.train[,-1]))
mydata_model1=nnet(cluster~.,data=data.train,size=4,rang=r,decay=5e-4,maxit=100)

mydata_predict=predict(mydata_model1,data.train)
mydata_predict=as.data.frame(mydata_predict)
mydata_predict$pre=apply(mydata_predict, 1, function(t) colnames(mydata_predict)[which.max(t)])

table(mydata_predict$pre,data.train$cluster)

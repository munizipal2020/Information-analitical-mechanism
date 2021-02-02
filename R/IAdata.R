#Метод измерения информационной важности показателей на основе корреляции с главными компонентами
inform_value<-function(Dataset){
  B<-matrix(NA, nrow=length(Dataset[,1]),ncol=length(Dataset))
  for (i in 1:length(Dataset[1,]))
  {
    B[,i]<-scale(Dataset[,i])
  }
  myPCA <- prcomp(B)
  gamma<-myPCA$sdev/sum(myPCA$sdev)
  result<-matrix(0,nrow=length(Dataset),ncol=1)
  for (i in 1:length(Dataset))
  {
    for (j in 1:length(Dataset))
    {result[i,1]<-result[i,1]+gamma[j]*abs(cor(myPCA$x[,j],B[,i]))}

  }
  rownames(result)<-colnames(Dataset)
  return(result)
}


#Агломеративные методы
demo_agnes<-function(E,level){
  library(cluster)
  library(factoextra)
  library(grid)
  library(gridExtra)
  library(ggplot2)
  library(lattice)
  for (i in c("euclidean","manhattan")){
    for (j in c("average", "single", "complete", "ward", "weighted")){
      Res<-agnes(E, diss=FALSE, metric=i,method = j, stand=FALSE)
      plot(Res, main=paste("Agnes: ","metric=",i,", method=",j,", stand= FALSE"))
      if (Res$ac<=level){
        jpeg(paste("Agnes: ","metric=",i,", method=",j,", stand= FALSE"," AC=", Res$ac,".jpg"))
        plot(Res, main=paste("Agnes: ","metric=",i,", method=",j,", stand= FALSE"))
        dev.off()
      }
      Res<-agnes(E, diss=FALSE, metric=i,method = j, stand=TRUE)
      plot(Res, main=paste("Agnes: ","metric=",i,", method=",j,", stand= TRUE"))
      if (Res$ac<=level){
        jpeg(paste("Agnes: ","metric=",i,", method=",j,", stand= TRUE"," AC=", Res$ac,".jpg"))
        plot(Res, main=paste("Agnes: ","metric=",i,", method=",j,", stand= TRUE"))
        dev.off()
      }
    }
  }
}

#Дивизимные методы
demo_diana<-function(E,level){
  for (i in c("euclidean","manhattan")){
    Res<-diana(Examp, diss=FALSE, metric=i, stand=FALSE)
    plot(Res, main=paste("Diana: ","metric=",i, ", stand= FALSE"))
    if (Res$dc<=level){
      jpeg(paste("Diana: ","metric=",i, "stand= FALSE"," DC=", Res$dc,".jpg"))
      plot(Res, main=paste("Diana: ","metric=",i, ", stand= FALSE"))
      dev.off()
    }
    Res<-diana(Examp, diss=FALSE, metric=i,stand=TRUE)
    plot(Res, main=paste("Diana: ","metric=",i, ", stand= TRUE"))
    if (Res$dc<=level){
      jpeg(paste("Diana: ","metric=",i, "stand= TRUE"," DC=", Res$dc,".jpg"))
      plot(Res, main=paste("Diana: ","metric=",i, ", stand= TRUE"))
      dev.off()
    }
  }
}

#Метод k-means
demo_kmeans<-function(E,k){
  for (i in c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen")){
    Res<-kmeans(E, centers=k, iter.max=100, nstart=10, algorithm=i)
    g<-fviz_cluster(Res, ellipse = TRUE, data=Examp, main=paste("kmeans: method=",i))
    jpeg(paste("kmeans: method = ",i,".jpg"))
    grid.arrange(g, nrow = 1)
    dev.off()
  }
}

#Метод главных компонент
demo_hcpc<-function(E, map3D=FALSE, tofile=FALSE){
  library(FactoMineR)
  res.pca <- PCA(E, graph=FALSE)
  for (i in c("euclidean","manhattan")){
    for (j in c("average", "single", "complete", "ward")){
      Res<-HCPC(res.pca, nb.clust=-1, iter.max=100, method=j, metric = i, graph = FALSE)
      if (map3D){
        if (tofile){
          jpeg(paste("hcpc: method = ",j,"metric=",i,".jpg"))
          plot(Res, choice = "3D.map", ind.names = TRUE,sub=paste("HCPC: metric=",i," method=",j))
          dev.off()
        }
        else{
          plot(Res, choice = "3D.map", ind.names = TRUE,sub=paste("HCPC: metric=",i," method=",j))
        }}
      else{
        if (tofile){
          jpeg(paste("hcpc: method = ",j,"metric=",i,".jpg"))
          plot(Res, choice = "tree", ind.names = TRUE, main=paste("HCPC: metric=",i," method=",j))
          dev.off()
        }
        else {plot(Res, choice = "tree", ind.names = TRUE, main=paste("HCPC: metric=",i," method=",j))}
      }
    }
  }
}

#Кластеризация организаций внутри муниципальных образований
demo_cluster_orgs<-function(path){
  setwd(path)
  library(cluster)
  library(factoextra)
  library(grid)
  library(gridExtra)
  library(ggplot2)
  library(lattice)
  for (i in list.files()){
    E<-read.csv(i)
    D<-data.frame(E[,2:length(E)])
    rownames(D)<-E[,1]
    print(paste("Длина: ",length(E[,1])))
    level=as.numeric(readline("Level="))
    for (k in c("euclidean","manhattan")){
      for (j in c("average", "single", "complete", "ward", "weighted")){
        Res<-agnes(D, diss=FALSE, metric=k,method = j, stand=FALSE)
        if (Res$ac<=level){
          jpeg(paste("Agnes: ","metric=",k,", method=",j,", stand= FALSE ",i, " AC=", Res$ac,".jpg"))
          plot(Res, main=paste("Agnes: ","metric=",k,", method=",j,", stand= FALSE, \n", i))
          dev.off()
        }
        Res<-agnes(D, diss=FALSE, metric=k,method = j, stand=TRUE)
        if (Res$ac<=level){
          jpeg(paste("Agnes: ","metric=",k,", method=",j,", stand= TRUE ", i," AC=", Res$ac,".jpg"))
          plot(Res, main=paste("Agnes: ","metric=",k,", method=",j,", stand= TRUE, \n", i))
          dev.off()
        }
      }
    }
  }
}

#Комплексирование первичных показателей муниципальных образований в рамках моделей кластеризации

model_munizipal_1<-function(M){
  U<-data.frame(as.numeric(M$m1), as.numeric(M$m2), as.numeric(M$m3),
                as.numeric(M$m4), as.numeric(M$m5), as.numeric(M$m6),
                as.numeric(M$m7), as.numeric(M$m8), as.numeric(M$m9),
                as.numeric(M$m10))
  colnames(U)<-colnames(M)
  rownames(U)<-rownames(M)
  U$m5[U$m5==0]<-0.01
  U$m9[U$m9==0]<-0.01
  U$m1[U$m1==0]<-0.01
  U$m4[U$m4==0]<-0.01
  U$m8[U$m8==0]<-0.01
  U$m7[U$m7==0]<-0.01
  U$m10[U$m10==0]<-0.01
  U$m3[U$m3==0]<-0.01
  R<-data.frame(U$m3/U$m5, U$m3/U$m9, U$m6/U$m1, U$m2/U$m4,
                U$m2/U$m1, U$m8/U$m4, U$m7/U$m8, U$m10/U$m7,
                U$m5/U$m8, U$m9/U$m10, U$m6/U$m3)
  colnames(R)<-c("m3/m5", "m3/m9", "m6/m1", "m2/m4",
                 "m2/m1", "m8/m4", "m7/m8", "m10/m7",
                 "m5/m8", "m9/m10", "m6/m3")
  rownames(R)<-rownames(M)
  return(R)
}

model_munizipal_2<-function(M){
  U<-data.frame(as.numeric(M$m1), as.numeric(M$m2), as.numeric(M$m3),
                as.numeric(M$m4), as.numeric(M$m5), as.numeric(M$m6),
                as.numeric(M$m7), as.numeric(M$m8), as.numeric(M$m9),
                as.numeric(M$m10))
  colnames(U)<-colnames(M)
  rownames(U)<-rownames(M)
  U$m5[U$m5==0]<-0.01
  U$m9[U$m9==0]<-0.01
  U$m1[U$m1==0]<-0.01
  U$m4[U$m4==0]<-0.01
  U$m8[U$m8==0]<-0.01
  U$m7[U$m7==0]<-0.01
  U$m10[U$m10==0]<-0.01
  U$m3[U$m3==0]<-0.01
  R<-data.frame(U$m3/U$m5, U$m6/U$m1, U$m2/U$m4,
                U$m8/U$m4, U$m7/U$m8, U$m10/U$m7,
                U$m5/U$m8, U$m9/U$m10, U$m6/U$m3)
  colnames(R)<-c("m3/m5",  "m6/m1", "m2/m4",
                 "m8/m4", "m7/m8", "m10/m7",
                 "m5/m8", "m9/m10", "m6/m3")
  rownames(R)<-rownames(M)
  return(R)
}


#Комплексирование первичных показателей организаций в рамках моделей кластеризации

model_org_1<-function(P,city){
  R<-data.frame(as.numeric(P$p1[P$city==city]), as.numeric(P$p2[P$city==city]),
                as.numeric(P$p3[P$city==city]), as.numeric(P$p4[P$city==city]),
                as.numeric(P$p5[P$city==city]), as.numeric(P$p6[P$city==city]))
  g<-data.frame(rownames(P),P$city)
  rownames(R)<-g$rownames.P.[g$P.city==city]
  colnames(R)<-c("p1","p2","p3","p4","p5","p6")
  return(R)
}

model_org_2<-function(P,city){
  U<-data.frame(as.numeric(P$p1), as.numeric(P$p2), as.numeric(P$p3),
                as.numeric(P$p4), as.numeric(P$p5), as.numeric(P$p6), P$city)
  colnames(U)<-colnames(P)
  rownames(U)<-rownames(P)
  U$p1[U$p1==0]<-0.01
  U$p3[U$p3==0]<-0.01
  R<-data.frame(U$p2[U$city==city]/U$p1[U$city==city],
                U$p4[U$city==city]/U$p3[U$city==city],
                U$p5[U$city==city], U$p6[U$city==city])
  g<-data.frame(rownames(P),P$city)
  rownames(R)<-g$rownames.P.[g$P.city==city]
  colnames(R)<-c("p2/p1","p4/p3","p5","p6")
  return(R)
}

model_org_3<-function(P,city){
  U<-data.frame(as.numeric(P$p1), as.numeric(P$p2), as.numeric(P$p3),
                as.numeric(P$p4), as.numeric(P$p5), as.numeric(P$p6), P$city)
  colnames(U)<-colnames(P)
  rownames(U)<-rownames(P)
  U$p1[U$p1==0]<-0.01
  U$p2[U$p2==0]<-0.01
  U$p3[U$p3==0]<-0.01
  U$p4[U$p4==0]<-0.01
  R<-data.frame(U$p5[U$city==city]/U$p1[U$city==city],
                U$p5[U$city==city]/U$p2[U$city==city],
                U$p6[U$city==city]/U$p3[U$city==city],
                U$p6[U$city==city]/U$p4[U$city==city])
  g<-data.frame(rownames(P),P$city)
  rownames(R)<-g$rownames.P.[g$P.city==city]
  colnames(R)<-c("p5/p1","p5/p2","p6/p3","p6/p4")
  return(R)
}

get_models_org<-function(model,P){
  for (i in P$city[!duplicated(P$city)])
  {
    write.csv(as.function(alist(a=,b=,get(model)(a,b)))(P,i),file=paste(i,".csv"))
  }
}

#Классификация организаций НПК наукоградов по пороговым значениям

test_Q<-function(org){
  p1_q<-c(18,21,25,34,48)
  p2_q<-c(6, 8, 10, 11, 18)
  p3_q<-c(58000, 82000, 98000, 98200, 133000)
  p4_q<-c(120000, 200000, 270000, 375000, 521000)
  p5_q<-c(1, 36, 37, 42, 139)
  p6_q<-c(105000, 160000, 200000, 283000, 420000)
  r<-matrix(NA, ncol=6, nrow=length(org[,1]))
  colnames(r)<-colnames(org)
  for (i in c(1:length(org[,1]))){
    if (org$p1[i]<p1_q[1]){
      r[i,1]<-"<Q1"
    }
    else if (org$p1[i]<p1_q[2]){
      r[i,1]<-"Q1"
    }
    else if (org$p1[i]<p1_q[3]){
      r[i,1]<-"Q2"
    }
    else if (org$p1[i]<p1_q[4]){
      r[i,1]<-"Q3"
    }
    else if (org$p1[i]<p1_q[5]){
      r[i,1]<-"Q4"
    }
    else{
      r[i,1]<-">Q4"
    }

    if (org$p2[i]<p2_q[1]){
      r[i,2]<-"<Q1"
    }
    else if (org$p2[i]<p2_q[2]){
      r[i,2]<-"Q1"
    }
    else if (org$p2[i]<p2_q[3]){
      r[i,2]<-"Q2"
    }
    else if (org$p2[i]<p2_q[4]){
      r[i,2]<-"Q3"
    }
    else if (org$p2[i]<p2_q[5]){
      r[i,2]<-"Q4"
    }
    else{
      r[i,2]<-">Q4"
    }

    if (org$p3[i]<p3_q[1]){
      r[i,3]<-"<Q1"
    }
    else if (org$p3[i]<p3_q[2]){
      r[i,3]<-"Q1"
    }
    else if (org$p3[i]<p3_q[3]){
      r[i,3]<-"Q2"
    }
    else if (org$p3[i]<p3_q[4]){
      r[i,3]<-"Q3"
    }
    else if (org$p3[i]<p3_q[5]){
      r[i,3]<-"Q4"
    }
    else{
      r[i,3]<-">Q4"
    }

    if (org$p4[i]<p4_q[1]){
      r[i,4]<-"<Q1"
    }
    else if (org$p4[i]<p4_q[2]){
      r[i,4]<-"Q1"
    }
    else if (org$p4[i]<p4_q[3]){
      r[i,4]<-"Q2"
    }
    else if (org$p4[i]<p4_q[4]){
      r[i,4]<-"Q3"
    }
    else if (org$p4[i]<p4_q[5]){
      r[i,4]<-"Q4"
    }
    else{
      r[i,4]<-">Q4"
    }

    if (org$p5[i]<p5_q[1]){
      r[i,5]<-"<Q1"
    }
    else if (org$p5[i]<p5_q[2]){
      r[i,5]<-"Q1"
    }
    else if (org$p5[i]<p5_q[3]){
      r[i,5]<-"Q2"
    }
    else if (org$p5[i]<p5_q[4]){
      r[i,5]<-"Q3"
    }
    else if (org$p5[i]<p5_q[5]){
      r[i,5]<-"Q4"
    }
    else{
      r[i,5]<-">Q4"
    }

    if (org$p6[i]<p6_q[1]){
      r[i,6]<-"<Q1"
    }
    else if (org$p6[i]<p6_q[2]){
      r[i,6]<-"Q1"
    }
    else if (org$p6[i]<p6_q[3]){
      r[i,6]<-"Q2"
    }
    else if (org$p6[i]<p6_q[4]){
      r[i,6]<-"Q3"
    }
    else if (org$p6[i]<p6_q[5]){
      r[i,6]<-"Q4"
    }
    else{
      r[i,6]<-">Q4"
    }
  }
  result<-data.frame(r)
  rownames(result)<-rownames(org)
  colnames(result)<-colnames(org)
  return(result)
}


test_Q_rank<-function(org){
  p1_q<-c(18,21,25,34,48)
  p2_q<-c(6, 8, 10, 11, 18)
  p3_q<-c(58000, 82000, 98000, 98200, 133000)
  p4_q<-c(120000, 200000, 270000, 375000, 521000)
  p5_q<-c(1, 36, 37, 42, 139)
  p6_q<-c(105000, 160000, 200000, 283000, 420000)
  r<-matrix(NA, ncol=6, nrow=length(org[,1]))
  colnames(r)<-colnames(org)
  for (i in c(1:length(org[,1]))){
    if (org$p1[i]<p1_q[1]){
      r[i,1]<-0
    }
    else if (org$p1[i]<p1_q[2]){
      r[i,1]<-1
    }
    else if (org$p1[i]<p1_q[3]){
      r[i,1]<-2
    }
    else if (org$p1[i]<p1_q[4]){
      r[i,1]<-3
    }
    else if (org$p1[i]<p1_q[5]){
      r[i,1]<-4
    }
    else{
      r[i,1]<-5
    }

    if (org$p2[i]<p2_q[1]){
      r[i,2]<-0
    }
    else if (org$p2[i]<p2_q[2]){
      r[i,2]<-1
    }
    else if (org$p2[i]<p2_q[3]){
      r[i,2]<-2
    }
    else if (org$p2[i]<p2_q[4]){
      r[i,2]<-3
    }
    else if (org$p2[i]<p2_q[5]){
      r[i,2]<-4
    }
    else{
      r[i,2]<-5
    }

    if (org$p3[i]<p3_q[1]){
      r[i,3]<-0
    }
    else if (org$p3[i]<p3_q[2]){
      r[i,3]<-1
    }
    else if (org$p3[i]<p3_q[3]){
      r[i,3]<-2
    }
    else if (org$p3[i]<p3_q[4]){
      r[i,3]<-3
    }
    else if (org$p3[i]<p3_q[5]){
      r[i,3]<-4
    }
    else{
      r[i,3]<-5
    }

    if (org$p4[i]<p4_q[1]){
      r[i,4]<-0
    }
    else if (org$p4[i]<p4_q[2]){
      r[i,4]<-1
    }
    else if (org$p4[i]<p4_q[3]){
      r[i,4]<-2
    }
    else if (org$p4[i]<p4_q[4]){
      r[i,4]<-3
    }
    else if (org$p4[i]<p4_q[5]){
      r[i,4]<-4
    }
    else{
      r[i,4]<-5
    }

    if (org$p5[i]<p5_q[1]){
      r[i,5]<-0
    }
    else if (org$p5[i]<p5_q[2]){
      r[i,5]<-1
    }
    else if (org$p5[i]<p5_q[3]){
      r[i,5]<-2
    }
    else if (org$p5[i]<p5_q[4]){
      r[i,5]<-3
    }
    else if (org$p5[i]<p5_q[5]){
      r[i,5]<-4
    }
    else{
      r[i,5]<-5
    }

    if (org$p6[i]<p6_q[1]){
      r[i,6]<-0
    }
    else if (org$p6[i]<p6_q[2]){
      r[i,6]<-1
    }
    else if (org$p6[i]<p6_q[3]){
      r[i,6]<-2
    }
    else if (org$p6[i]<p6_q[4]){
      r[i,6]<-3
    }
    else if (org$p6[i]<p6_q[5]){
      r[i,6]<-4
    }
    else{
      r[i,6]<-5
    }
  }
  result<-data.frame(r)
  rownames(result)<-rownames(org)
  colnames(result)<-colnames(org)
  return(result)
}

test_model2_Q<-function(org){
  p5_p1<-c(1,2,3,5,9)
  p5_p2<-c(1,2,4,5,9)
  p6_p3<-c(0.1, 1.1, 1.4, 3, 5.4)
  p6_p4<-c(1, 15, 27, 47, 97)
  r<-matrix(NA, ncol=4, nrow=length(org[,1]))
  colnames(r)<-c("p5/p1","p5/p2", "p6/p3","p6/p4")
  p<-as.numeric(org$p5)/as.numeric(org$p1)
  p[p=="NaN"]=0; p[p=="Inf"]=max(p[p!=Inf])
  for (i in c(1:length(org[,1]))){
    if (p[i]<p5_p1[1]){
      r[i,1]<-0
    }
    else if (p[i]<p5_p1[2]){
      r[i,1]<-1
    }
    else if (p[i]<p5_p1[3]){
      r[i,1]<-2
    }
    else if (p[i]<p5_p1[4]){
      r[i,1]<-3
    }
    else if (p[i]<p5_p1[5]){
      r[i,1]<-4
    }
    else{
      r[i,1]<-5
    }
  }
  p<-as.numeric(org$p5)/as.numeric(org$p2)
  p[p=="NaN"]=0; p[p=="Inf"]=max(p[p!=Inf])
  for (i in c(1:length(org[,1]))){
    if (p[i]<p5_p2[1]){
      r[i,2]<-0
    }
    else if (p[i]<p5_p2[2]){
      r[i,2]<-1
    }
    else if (p[i]<p5_p2[3]){
      r[i,2]<-2
    }
    else if (p[i]<p5_p2[4]){
      r[i,2]<-3
    }
    else if (p[i]<p5_p2[5]){
      r[i,2]<-4
    }
    else{
      r[i,2]<-5
    }
  }

  p<-as.numeric(org$p6)/as.numeric(org$p3)
  p[p=="NaN"]=0; p[p=="Inf"]=max(p[p!=Inf])

  for (i in c(1:length(org[,1]))){
    if (p[i]<p6_p3[1]){
      r[i,3]<-0
    }
    else if (p[i]<p6_p3[2]){
      r[i,3]<-1
    }
    else if (p[i]<p6_p3[3]){
      r[i,3]<-2
    }
    else if (p[i]<p6_p3[4]){
      r[i,3]<-3
    }
    else if (p[i]<p6_p3[5]){
      r[i,3]<-4
    }
    else{
      r[i,3]<-5
    }
  }

  p<-as.numeric(org$p6)/as.numeric(org$p4)
  p[p=="NaN"]=0; p[p=="Inf"]=max(p[p!=Inf])
  for (i in c(1:length(org[,1]))){
    if (p[i]<p6_p4[1]){
      r[i,4]<-0
    }
    else if (p[i]<p6_p4[2]){
      r[i,4]<-1
    }
    else if (p[i]<p6_p4[3]){
      r[i,4]<-2
    }
    else if (p[i]<p6_p4[4]){
      r[i,4]<-3
    }
    else if (p[i]<p6_p4[5]){
      r[i,4]<-4
    }
    else{
      r[i,4]<-5
    }
  }
  result<-data.frame(r)
  rownames(result)<-rownames(org)
  colnames(result)<-c("p5/p1","p5/p2", "p6/p3","p6/p4")
  return(result)
}

test_munizipal_Q<-function(org){
  m2_m4<-c(160,310,356,396,580)
  m10_m7<-c(0.14, 0.23,0.3, 0.41, 0.49)
  m7_m8<-c(0.4, 0.46, 0.6, 0.7, 0.82)
  r<-matrix(NA, ncol=3, nrow=length(org[,1]))
  colnames(r)<-c("m2/m4","m10/m7", "m7/m8")
  p<-as.numeric(org$`m2/m4`)
  for (i in c(1:length(org[,1]))){
    if (p[i]<m2_m4[1]){
      r[i,1]<-0
    }
    else if (p[i]<m2_m4[2]){
      r[i,1]<-1
    }
    else if (p[i]<m2_m4[3]){
      r[i,1]<-2
    }
    else if (p[i]<m2_m4[4]){
      r[i,1]<-3
    }
    else if (p[i]<m2_m4[5]){
      r[i,1]<-4
    }
    else{
      r[i,1]<-5
    }
  }
  p<-as.numeric(org$`m10/m7`)
  for (i in c(1:length(org[,1]))){
    if (p[i]<m10_m7[1]){
      r[i,2]<-0
    }
    else if (p[i]<m10_m7[2]){
      r[i,2]<-1
    }
    else if (p[i]<m10_m7[3]){
      r[i,2]<-2
    }
    else if (p[i]<m10_m7[4]){
      r[i,2]<-3
    }
    else if (p[i]<m10_m7[5]){
      r[i,2]<-4
    }
    else{
      r[i,2]<-5
    }
  }
  p<-as.numeric(org$`m7/m8`)
  for (i in c(1:length(org[,1]))){
    if (p[i]<m7_m8[1]){
      r[i,3]<-0
    }
    else if (p[i]<m7_m8[2]){
      r[i,3]<-1
    }
    else if (p[i]<m7_m8[3]){
      r[i,3]<-2
    }
    else if (p[i]<m7_m8[4]){
      r[i,3]<-3
    }
    else if (p[i]<m7_m8[5]){
      r[i,3]<-4
    }
    else{
      r[i,3]<-5
    }
  }
  result<-data.frame(r)
  rownames(result)<-rownames(org)
  colnames(result)<-colnames(r)
  return(result)
}





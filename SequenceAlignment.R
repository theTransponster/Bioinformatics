#Dynamic programming for global alignment
#Mariana Sofia Flores 


#1 Inicializar la matriz de tama??o m+1 filas y n+1 columnas
#siendo m y n la longitud de cadenas de nucleotidos 


#ejemplo
n=c('A','G','C','T','A','G','C','G','A','T','C')
m=c('T','A','G','G','G','A','A','C','G')
gap=0 #Gap inicial 
w=0


matriz=matrix(0,length(m)+1,length(n)+1)
matriz[1,1]=0
matriz[1:length(m)+1,1]=gap
matriz[1,1:length(n)+1]=gap

#2 Se crea una funcion con un vector para almacenar los scores:
#length(m)+1
alineacion <- function(m, n, g, w, S1, S2){
  gap=g
  matriz=matrix(0,length(m)+1,length(n)+1)
  matriz[1,1]=0
  matriz[1:length(m)+1,1]=gap
  matriz[1,1:length(n)+1]=gap
  S=0
  for (i in 2:(length(m)+1)){
    for(j in 2:(length(n)+1)){
      if(m[i-1]==n[j-1]){
        S=S1
      }else{
        S=S2
      }
      score=matrix(0,3,1)
      score[1,1]=matriz[i-1,j-1]+S
      score[2,1]=matriz[i,j-1]+w
      score[3,1]=matriz[i-1,j]+w
      matriz[i,j]=max(score)
    }
  }
  Return(matriz)
}


#3 Traceback para determinar la alineacion con el maximo score
#comenzando desde la posicion con el maximo score
row=which(matriz == max(matriz), arr.ind = TRUE)[1]
col=which(matriz == max(matriz), arr.ind = TRUE)[2]

al=matrix(0,2,(max(length(m),length(n))))
al[2,dim(al)[2]]=n[c(col-1)]
al[1,dim(al)[2]]=m[c(row-1)]

for (i in length(m)){
  for (j in length(n)){
    if(m[i]==n[j]){
      al[1,i]<-n[i]
      al[3,i]<-m[i]
      al[2,i]<-'|'
    } else{
      al[1,i]<-n[i]
      al[2,i]<-' '
      al[3,i]<-'-'   	
    }
  }
}



#Prueba con 100 Permutaciones 
library(gtools)

per=sample(1:length(m), length(m), replace=FALSE)
mnew=m[per] 
for (k in 2:100){
  per=sample(1:9, 9, replace=FALSE)
  new=m[per] 
  mnew<-rbind(mnew,new)
}

scores<-matrix(0,100,1)


for (l in 1:100){
  seq=mnew[l,1:length(m)]
  
  matriz=matrix(0,length(m)+1,length(n)+1)
  matriz[1,1]=0
  matriz[1:length(m)+1,1]=gap
  matriz[1,1:length(n)+1]=gap
  
  for (i in 2:(length(seq)+1)){
    for(j in 2:(length(n)+1)){
      if(seq[i-1]==n[j-1]){
        S=1
      }else{
        S=0
      }
      score=matrix(0,3,1)
      score[1,1]=matriz[i-1,j-1]+S
      score[2,1]=matriz[i,j-1]+w
      score[3,1]=matriz[i-1,j]+w
      matriz[i,j]=max(score)
    }
  }
  scores[l]<-max(matriz)
  
}



#Usando la funcion:


for (l in 1:100){
  seq=mnew[l,1:length(m)]
  matriz<-alineacion(mnew,n, 0, 0, 1, 0)  
  scores[l]<-max(matriz)	
}





#Los scores son iguales a pesar de las permutaciones a menos que
#se cambie el numero establecido para el gap inicial, cuando se
#cambia a 0, el score maximo si varia. Esta variacion de scores
#se puede observar en el siguiente histograma

#install.packages("ggplot2")
library(ggplot2)
hist(scores)


#Calculamos el p value (considerando que es una distribucion normal)
p<-pnorm(-abs(scores))
mean1<-mean(p)



##----
#lo hacemos con un gap de 5 y 1000 permutaciones 

per=sample(1:length(m), length(m), replace=FALSE)
mnew=m[per] 
for (k in 2:1000){
  per=sample(1:9, 9, replace=FALSE)
  new=m[per] 
  mnew<-rbind(mnew,new)
}

scores<-matrix(0,1000,1)
for (l in 1:1000){
  seq=mnew[l,1:length(m)]
  matriz<-alineacion(mnew,n, 0, 5, 1, 0)  
  scores[l]<-max(matriz)	
}




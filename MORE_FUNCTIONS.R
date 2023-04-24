################################################################################

cpdag_list <- function(list.inc,E){    # E: end of burnIn phase
  L <- list()
  G <- list()
  
 
  nodes <- dim(list.inc[[1]][[1]])[1]
  
  
  mat.sum <- matrix(numeric(nodes*nodes),nrow=nodes)
  
 
  for (i in E:length(list.inc[[1]])){
   
    k <- cpdag(list.inc[[1]][[i]])
    dummy <- matrix(numeric(nodes*nodes),nrow=nodes)
    
    if(length(nrow(k))!=0){
      dummy[k[,1]] <- k[,5]
      L[[i]] <- dummy
    }
    if(length(nrow(k))==0 && length(k)>0){
      dummy[k[1]] <- k[5]
      L[[i]] <- dummy
    }
    mat.com <-matrix(numeric(nodes*nodes),nrow=nodes)
    mat.re <- matrix(numeric(nodes*nodes),nrow=nodes)
    com <- which(L[[i]]>0)
    re  <- which(L[[i]]<0)
    
    mat.com[com] <- 1
    mat.re[re]   <- 1
    
    mat <- mat.com + mat.re + t(mat.re)
    
    G[[i]] <- mat
    
    mat.sum <- mat.sum + mat
  }
  return(list(L,G, (mat.sum/(length(list.inc[[1]])- E+1))))
}


################################################################################

extract_cpdag_from_dag <- function(true_incidence){    
  
  L <- list()
  
  nodes <- dim(true_incidence)[1]
  
  k <- cpdag(true_incidence)
  
  
  dummy <- matrix(numeric(nodes*nodes),nrow=nodes)
  
  if(length(nrow(k))!=0){
    dummy[k[,1]] <- k[,5]
    L <- dummy
  }
  
  if(length(nrow(k))==0 && length(k)>0){
    dummy[k[1]] <- k[5]
    L <- dummy
  }
  
  mat.com <- matrix(numeric(nodes*nodes),nrow=nodes)
  mat.re  <- matrix(numeric(nodes*nodes),nrow=nodes)
  
  com <- which(L>0)
  re  <- which(L<0)
  
  mat.com[com] <- 1
  mat.re[re] <- 1
  mat <- mat.com + mat.re + t(mat.re)
  
  return(mat)
}


##########################################################################################

compute_SHD <- function(MAT_scores,TRUE_DAG){   

  n_dim   = dim(TRUE_DAG)
  n_nodes = n_dim[1]
  
  TRUE_CPDAG = cpdag(TRUE_DAG)

  DAG_REV = matrix(0,n_nodes,n_nodes)
  DAG_COM = matrix(0,n_nodes,n_nodes)
  
  ind_rev = which(TRUE_CPDAG[,5]==-1)
  ind_com = which(TRUE_CPDAG[,5]==+1)
  
  
  for (i in 1:length(ind_rev))
  {
    DAG_REV[TRUE_CPDAG[ind_rev[i],2],TRUE_CPDAG[ind_rev[i],3]] = 1
  }
  
  for (i in 1:length(ind_com))
  {
    DAG_COM[TRUE_CPDAG[ind_com[i],2],TRUE_CPDAG[ind_com[i],3]] = 1
  }
  
  
  
  DIRECTED1   = DAG_COM   
  UNDIRECTED1 = DAG_REV + t(DAG_REV)

  
  IND = which(MAT_scores >= 0.5)


  MAT_predicted      = matrix(0,n_nodes,n_nodes)
  MAT_predicted[IND] = 1;


  MAT_NEW = MAT_predicted + t(MAT_predicted)

  ind_undirected = which(MAT_NEW==2)
  

  UNDIRECTED2    = matrix(0,n_nodes,n_nodes)

  UNDIRECTED2[ind_undirected] = 1  

  DIRECTED2 = MAT_predicted 
  
  SHD = sum(abs(UNDIRECTED1 - UNDIRECTED2))/2 
  
  ind1 = which(UNDIRECTED1==1)
  ind2 = which(UNDIRECTED2==1)
  
  DIRECTED1[ind1] = 0
  DIRECTED2[ind2] = 0
  
  DIRECTED1[ind2] = 0
  DIRECTED2[ind1] = 0
  
  DIFF = abs(DIRECTED1 - DIRECTED2)
  DIFF2 = DIFF + t(DIFF)
  
  add = length(which(DIFF2[upper.tri(DIFF2)]>0))
  
  SHD = SHD + add
  
  
  return(SHD)

}


###############################################################################

compute_AUROC <- function(Aposteriori_Matrix,TRUE_DAG){

n_nodes = dim(TRUE_DAG)[1]

True_Matrix = extract_cpdag_from_dag(TRUE_DAG)   


for (i in 1:n_nodes){True_Matrix[i,i]=-1}


n_non_edges   = length(which(True_Matrix==0))
n_edges       = length(which(True_Matrix==1))

Aposterioris   = NULL

for (i in 1:n_nodes){
for (j in 1:n_nodes){
if (i!=j){Aposterioris = c(Aposterioris, Aposteriori_Matrix[i,j])}
}
}

Aposterioris = sort(Aposterioris,decreasing=TRUE)

APO_values = Aposterioris[1]

for (i in 2:length(Aposterioris)){
if (Aposterioris[i]==APO_values[length(APO_values)]){}
else
  {APO_values = c(APO_values,Aposterioris[i])}
}



ROC_x = 0
ROC_y = 0


MATRIX = matrix(0,n_nodes,n_nodes)

for (i in 1:length(APO_values)){

indicis = which(Aposteriori_Matrix>=APO_values[i])
MATRIX[indicis] = 1

TP = length(which(MATRIX==1 & True_Matrix==1))
TN = length(which(MATRIX==0 & True_Matrix==0)) 

Sensitivity = TP/n_edges;
inv_Specif  = 1 - (TN/n_non_edges)


ROC_y = c(ROC_y, Sensitivity)
ROC_x = c(ROC_x, inv_Specif)
}


ROC_x = c(ROC_x,1)
ROC_y = c(ROC_y,1)

auroc = trapz(ROC_x,ROC_y)

return(auroc)
}







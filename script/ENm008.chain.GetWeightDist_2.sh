#!/bin/bash - 
# Generate contact node weighted distance

paramlist=(1_1_1_1)
OutDirRoot=/tmp

MaxNum=1000
for par in ${paramlist[@]}; do
  par=par.$par
  SampleSize=`ls $par`
#  SampleSize=(K562.p1500_d3300)
  for Sample in ${SampleSize[@]}; do
    # 0001.pts has maximum weight
    MaxLogWeight=$(find $par/$Sample -name 0001.dst -exec grep LogWeight {} \; | sort -k3,3n | tail -n1 | cut -d" " -f3)
    
    FileListFile=$par.$Sample.ls.txt
    find $par/$Sample -name \*.dst > $FileListFile

    OutDir=$OutDirRoot/$par/$Sample
    if [ ! -d $OutDir ]; then install -d $OutDir; fi 
    R --vanilla --slave <<-RSCRIPT
MaxNum=$MaxNum
MaxLogWeight = $MaxLogWeight

FileNamesFile="$FileListFile"
OutDir = "$OutDir"

FileNames = as.vector(unlist(read.table(FileNamesFile)))

#FileNames = FileNames[1:5]

NewFileNames = NULL
NumCopyVec = NULL
for (FileName in FileNames){
  LogWeight = read.table(FileName,comment.char="",nrows=1)[,3]
  Weight = exp(LogWeight - MaxLogWeight)
  NumCopies = round(MaxNum * Weight)
  if (NumCopies != 0){ 
    NumCopyVec = c(NumCopyVec,NumCopies)
    NewFileNames = c(NewFileNames, FileName) 
  }
}

NodeNodeVec = read.table(FileNames[1])[,1:2]
iNumFile = length(NewFileNames)
DistMat = matrix(0, nrow=nrow(NodeNodeVec),ncol=iNumFile)

for (i in 1:iNumFile){
  cat(i,"/",iNumFile,"\n")
  flush.console()
  FileName = NewFileNames[i]
  DIST = read.table(FileName)[,3]
  DistMat[,i] = DIST
  
}

for( iRow in 1:nrow(NodeNodeVec)){
  FN = paste(OutDir,"/",NodeNodeVec[iRow,1],"_",NodeNodeVec[iRow,2],".txt",sep="")
  write.table(file=FN,cbind(sprintf("%.3f",DistMat[iRow,]),NumCopyVec),quote=F,row.names=F,col.names=F)
}
RSCRIPT
  done
done
